#!/usr/bin/env python3
"""
Multi-seed antibody-antigen prediction via Biohub Forge API.
No local GPU needed — compute runs on Biohub servers.

Usage:
    export BIOHUB_TOKEN="your_token"

    python predict_forge.py \\
        --heavy   VH_SEQ \\
        --light   VL_SEQ \\
        --antigen AG_SEQ \\
        --seeds   5 \\
        --out     complex.cif

Rate limit: free tier = 100 credits/day (~100 predictions).
Use --workers 1 (default) to avoid per-minute cap.
"""

import argparse
import csv
import os
import re
import sys
import concurrent.futures


def parse_chai_csv(path: str) -> list[tuple[str, int]]:
    contacts = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            for chain_col, res_col in [("chainA", "res_idxA"), ("chainB", "res_idxB")]:
                if row[chain_col].strip() == "C":
                    match = re.match(r"[A-Z](\d+)", row[res_col].strip())
                    if match:
                        contacts.append(("antigen", int(match.group(1)) - 1))
    return list(dict.fromkeys(contacts))


def _run_one_seed(args):
    token, model_name, inputs, pocket, num_loops, num_sampling_steps, seed_idx = args
    from esm.sdk.forge import SequenceStructureForgeInferenceClient
    from esm.utils.structure.input_builder import StructurePredictionInput
    from esm.sdk.api import FoldingConfig, ESMProteinError

    client = SequenceStructureForgeInferenceClient(token=token, model=model_name)
    config = FoldingConfig(num_loops=num_loops, num_sampling_steps=num_sampling_steps,
                           include_pae=True)
    result = client.fold_all_atom(
        StructurePredictionInput(sequences=inputs, pocket=pocket), config=config
    )
    if isinstance(result, ESMProteinError):
        print(f"  seed {seed_idx:>3}: ERROR — {result}", flush=True)
        return None
    iptm = result.iptm or 0.0
    ptm  = result.ptm  or 0.0
    print(f"  seed {seed_idx:>3}: ipTM={iptm:.3f}  pTM={ptm:.3f}", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(
        description="ESMFold2 multi-seed prediction via Biohub API",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--heavy",        default=None)
    parser.add_argument("--light",        default=None)
    parser.add_argument("--antigen",      default=None)
    parser.add_argument("--antigen2",     default=None,
                        help="Second antigen chain (e.g. HA2)")
    parser.add_argument("--out",          default="complex.cif")
    parser.add_argument("--seeds",        type=int, default=5)
    parser.add_argument("--loops",        type=int, default=20)
    parser.add_argument("--diff-steps",   type=int, default=100)
    parser.add_argument("--workers",      type=int, default=1,
                        help="Parallel API calls (default 1 — avoids rate limit)")
    parser.add_argument("--token",        default=None,
                        help="Biohub token (or set BIOHUB_TOKEN env var)")
    parser.add_argument("--contacts-csv", default=None)

    args = parser.parse_args()

    token = args.token or os.environ.get("BIOHUB_TOKEN", "")
    if not token:
        sys.exit("No token. Set BIOHUB_TOKEN or use --token.")

    from esm.utils.structure.input_builder import ProteinInput, PocketConditioning
    from esm.utils.constants.models import ESMFOLD2_FAST

    provided = {k: v.strip() for k, v in
                [("heavy", args.heavy), ("light", args.light),
                 ("antigen", args.antigen), ("antigen2", args.antigen2)]
                if v}
    if not provided:
        sys.exit("Provide at least one of --heavy / --light / --antigen.")

    total_aa = sum(len(s) for s in provided.values())
    print(f"Chains : " + "  ".join(f"{n}={len(s)}aa" for n, s in provided.items()))
    print(f"Total  : {total_aa} residues  (limit: 768)")
    if total_aa > 768:
        sys.exit(f"Total length {total_aa} exceeds 768 limit. Trim sequences.")

    inputs = [ProteinInput(id=name, sequence=seq) for name, seq in provided.items()]

    pocket = None
    if args.contacts_csv:
        raw = parse_chai_csv(args.contacts_csv)
        binder = "heavy" if "heavy" in provided else list(provided.keys())[0]
        pocket = PocketConditioning(binder_chain_id=binder, contacts=raw)
        print(f"PocketConditioning: {raw}")

    print(f"Model  : {ESMFOLD2_FAST}  |  loops={args.loops}  |  seeds={args.seeds}\n")

    tasks = [(token, ESMFOLD2_FAST, inputs, pocket, args.loops, args.diff_steps, i)
             for i in range(1, args.seeds + 1)]

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for r in pool.map(_run_one_seed, tasks):
            if r is not None:
                results.append(r)

    if not results:
        sys.exit("All seeds failed.")

    best = max(results, key=lambda r: r.iptm or 0.0)
    iptm_vals = sorted([round(r.iptm or 0.0, 3) for r in results], reverse=True)
    print(f"\n{len(results)}/{args.seeds} seeds succeeded.")
    print(f"iPTM scores: {iptm_vals}")
    print(f"Best — ipTM={best.iptm:.3f}  pTM={best.ptm:.3f}  "
          f"mean pLDDT={float(best.plddt.mean()):.3f}")

    with open(args.out, "w") as f:
        f.write(best.complex.to_mmcif())
    print(f"Saved → {args.out}")


if __name__ == "__main__":
    main()
