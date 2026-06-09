#!/usr/bin/env python3
"""
Local ESMFold2 inference — no API, no credits, no rate limits.
Runs entirely on your GPUs using downloaded weights.

Usage:
    python predict_local.py \\
        --heavy   EVQLVESGGGLVK... \\
        --light   DIVMTQSPDSLA... \\
        --antigen DQICIGYHANN... \\
        --out     complex.cif \\
        --seeds   25

    # Load contacts from Chai CSV for pocket conditioning:
    python predict_local.py ... --contacts-csv /path/to/constraints.csv

Weights location (set once):
    WEIGHTS_DIR env var, or pass --weights-dir
    Default: ./weights
"""

import argparse
import csv
import os
import re
import sys

WEIGHTS_DIR = os.environ.get(
    "WEIGHTS_DIR",
    os.path.join(os.path.dirname(__file__), "weights"),
)
ESMFOLD2_PATH      = os.path.join(WEIGHTS_DIR, "ESMFold2")
ESMFOLD2_FAST_PATH = os.path.join(WEIGHTS_DIR, "ESMFold2-Fast")
ESMC6B_PATH        = os.path.join(WEIGHTS_DIR, "ESMC-6B")


def parse_chai_csv(path: str) -> list[tuple[str, int]]:
    """Extract 0-indexed antigen contact residues from a Chai restraint CSV."""
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


def main():
    parser = argparse.ArgumentParser(
        description="Local ESMFold2 complex prediction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--heavy",       default=None)
    parser.add_argument("--light",       default=None)
    parser.add_argument("--antigen",     default=None)
    parser.add_argument("--out",         default=None,
                        help="Path for best CIF (default: auto-named in output/)")
    parser.add_argument("--seeds",       type=int,   default=25)
    parser.add_argument("--loops",       type=int,   default=20)
    parser.add_argument("--diff-steps",  type=int,   default=100)
    parser.add_argument("--lm-dropout",  type=float, default=0.3)
    parser.add_argument("--lm-mask-pct", type=float, default=0.0)
    parser.add_argument("--device",      default="cuda:0")
    parser.add_argument("--fast",        action="store_true", default=True,
                        help="Use ESMFold2-Fast (default: True)")
    parser.add_argument("--full",        action="store_true",
                        help="Use full ESMFold2 instead of Fast")
    parser.add_argument("--weights-dir", default=WEIGHTS_DIR)
    parser.add_argument("--out-dir",     default="output/all_seeds",
                        help="Save all seed CIFs (default: output/all_seeds)")

    pocket_grp = parser.add_mutually_exclusive_group()
    pocket_grp.add_argument("--contacts", nargs="+", metavar="CHAIN:RESIDX")
    pocket_grp.add_argument("--contacts-csv", metavar="PATH")

    args = parser.parse_args()

    fold_model_name = "ESMFold2" if args.full else "ESMFold2-Fast"
    esmfold2_path   = os.path.join(args.weights_dir, fold_model_name)
    esmc6b_path     = os.path.join(args.weights_dir, "ESMC-6B")

    for p in [esmfold2_path, esmc6b_path]:
        if not os.path.isdir(p):
            sys.exit(f"Weights not found at {p}\n"
                     f"Run the download first or set --weights-dir")

    provided = {k: v.strip() for k, v in
                [("heavy", args.heavy), ("light", args.light), ("antigen", args.antigen)]
                if v}
    if not provided:
        sys.exit("Provide at least one of --heavy / --light / --antigen.")

    # ── Pocket conditioning ────────────────────────────────────────────────
    pocket = None
    if args.contacts_csv:
        raw_contacts = parse_chai_csv(args.contacts_csv)
        print(f"Contacts from CSV: {raw_contacts}")
    elif args.contacts:
        raw_contacts = [(c.split(":")[0], int(c.split(":")[1])) for c in args.contacts]
    else:
        raw_contacts = []

    if raw_contacts:
        from esm.models.esmfold2 import StructurePredictionInput as SPI2
        # PocketConditioning lives in the same namespace
        from esm.utils.structure.input_builder import PocketConditioning
        binder = "heavy" if "heavy" in provided else list(provided.keys())[0]
        pocket = PocketConditioning(binder_chain_id=binder, contacts=raw_contacts)
        print(f"PocketConditioning: binder={binder}, contacts={raw_contacts}")

    # ── Load model ─────────────────────────────────────────────────────────
    import torch
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model
    from esm.models.esmfold2 import (
        ESMFold2InputBuilder, ProteinInput, StructurePredictionInput
    )

    print(f"\nLoading {fold_model_name} from {esmfold2_path} ...")
    model = ESMFold2Model.from_pretrained(
        esmfold2_path,
        dtype=torch.bfloat16,
        device_map=args.device,
    ).eval()
    print(f"Loading ESMC-6B backbone from {esmc6b_path} ...")
    model.load_esmc(esmc6b_path, precision="bf16")
    print(f"Model ready on {args.device}\n")

    builder = ESMFold2InputBuilder()

    # ── Build input ────────────────────────────────────────────────────────
    inputs = StructurePredictionInput(
        sequences=[ProteinInput(id=name, sequence=seq)
                   for name, seq in provided.items()],
        pocket=pocket,
    )

    total_aa = sum(len(s) for s in provided.values())
    print(f"Chains : " + "  ".join(f"{n}={len(s)}aa" for n, s in provided.items()))
    print(f"Total  : {total_aa} residues")
    print(f"Seeds  : {args.seeds}  |  loops={args.loops}  |  "
          f"diff_steps={args.diff_steps}  |  lm_dropout={args.lm_dropout}\n")

    # ── Set up run output directory ────────────────────────────────────────
    from datetime import datetime
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "output", timestamp
    )
    seeds_dir = os.path.join(run_dir, "all_seeds")
    os.makedirs(seeds_dir, exist_ok=True)

    # ── Run seeds ──────────────────────────────────────────────────────────
    results = []
    for seed in range(args.seeds):
        print(f"Seed {seed+1:>3}/{args.seeds} ...", end=" ", flush=True)
        with torch.autocast(device_type="cuda", dtype=torch.bfloat16):
            result = builder.fold(
                model,
                inputs,
                num_loops=args.loops,
                num_sampling_steps=args.diff_steps,
                lm_dropout=args.lm_dropout,
                lm_mask_pct=args.lm_mask_pct if args.lm_mask_pct != 0.0 else None,
                seed=seed,
            )
        iptm  = float(result.iptm) if result.iptm  is not None else 0.0
        ptm   = float(result.ptm)  if result.ptm   is not None else 0.0
        plddt = float(result.plddt.mean()) if result.plddt is not None else 0.0
        print(f"ipTM={iptm:.3f}  pTM={ptm:.3f}  pLDDT={plddt:.3f}")
        results.append((iptm, ptm, plddt, seed+1, result))

        seed_path = os.path.join(seeds_dir, f"seed_{seed+1:04d}.cif")
        with open(seed_path, "w") as f:
            f.write(result.complex.to_mmcif())

    # ── Rank and save ──────────────────────────────────────────────────────
    results.sort(key=lambda x: x[0], reverse=True)

    # Rename seed files to ranked names
    for rank, (iptm, ptm, plddt, seed, _) in enumerate(results, 1):
        src = os.path.join(seeds_dir, f"seed_{seed:04d}.cif")
        dst = os.path.join(seeds_dir, f"rank{rank:02d}_seed{seed:02d}_ipTM{iptm:.3f}.cif")
        os.rename(src, dst)

    best_iptm, best_ptm, best_plddt, best_seed, best_result = results[0]

    # Save best at top level of run dir
    best_path = args.out or os.path.join(
        run_dir, f"best_seed{best_seed:02d}_ipTM{best_iptm:.3f}.cif"
    )
    with open(best_path, "w") as f:
        f.write(best_result.complex.to_mmcif())

    # Write summary CSV
    with open(os.path.join(run_dir, "summary.csv"), "w") as f:
        f.write("rank,seed,ipTM,pTM,pLDDT\n")
        for rank, (iptm, ptm, plddt, seed, _) in enumerate(results, 1):
            f.write(f"{rank},{seed},{iptm:.3f},{ptm:.3f},{plddt:.3f}\n")

    iptm_vals = [round(r[0], 3) for r in results]
    print(f"\n{'='*55}")
    print(f"iPTM scores (ranked): {iptm_vals}")
    print(f"Best  seed {best_seed:2d} — ipTM={best_iptm:.3f}  pTM={best_ptm:.3f}  pLDDT={best_plddt:.3f}")
    print(f"\nRun saved to: {run_dir}/")
    print(f"  best_seed{best_seed:02d}_ipTM{best_iptm:.3f}.cif  ← open this one")
    print(f"  summary.csv")
    print(f"  all_seeds/  ({args.seeds} ranked CIFs)")


if __name__ == "__main__":
    main()
