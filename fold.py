#!/usr/bin/env python
"""Fold an antibody (H + L) + antigen complex with ESMFold2 — validated recipe by default.

The defaults ARE the recipe: ESMFold2-Fast, single-sequence (no MSA, no template, no
constraints), 20 loops, 100 sampling steps, 1 diffusion sample, lm_dropout 0.3, fp32,
seeds 0-24. So a bare `python fold.py --target <dir>` is already the right run.

Every seed is saved (CIF + PAE + pLDDT + meta), then score.py runs automatically and
writes <out_dir>/scores.csv, ranked by how many seeds agree on the epitope — never by
ipTM, which does not track epitope correctness for antibody-antigen (see README).

Input — pick one:
  --target DIR   directory holding heavy.txt / light.txt / antigen.txt
                 (+ optional reference.cif, used as the DockQ reference)
  --fasta FILE   FASTA whose record ids are H / L / A (heavy / light / antigen)
  --heavy S --light S --antigen S    raw sequences, or paths to files holding one

Usage:
  python fold.py --target targets/my_ab --gpu 0
  python fold.py --heavy EVQ... --light DIV... --antigen DQI... --tag myAb --gpu 1
"""
import argparse
import csv
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).parent.resolve()

AA = set("ACDEFGHIKLMNPQRSTVWYXBZUO")


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def _clean(seq, what):
    seq = "".join(seq.split()).upper()
    bad = sorted(set(seq) - AA)
    if not seq:
        sys.exit(f"{what}: empty sequence")
    if bad:
        sys.exit(f"{what}: not a protein sequence (bad characters {bad})")
    return seq


def read_seq(value, what):
    """A raw sequence, or a path to a file holding one (plain text or single-record FASTA)."""
    if value is None:
        return None
    p = Path(value)
    if len(value) < 4096 and p.exists() and p.is_file():
        lines = [ln for ln in p.read_text().splitlines() if not ln.startswith(">")]
        return _clean("".join(lines), f"{what} ({p})")
    return _clean(value, what)


def read_fasta(path):
    """{record_id: sequence} in file order."""
    out, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip().split()[0] if line[1:].strip() else f"seq{len(out)}"
            out[name] = ""
        elif name is not None:
            out[name] += line.strip()
    return out


HEAVY_IDS = {"h", "heavy", "vh", "hc"}
LIGHT_IDS = {"l", "light", "vl", "lc", "kappa", "lambda"}
ANTIGEN_IDS = {"a", "ag", "antigen", "target"}


def chains_from_fasta(path):
    """Map FASTA records onto heavy/light/antigen by id, else positionally (H, L, A)."""
    recs = read_fasta(path)
    if not recs:
        sys.exit(f"No records in {path}")
    picked = {}
    for name, seq in recs.items():
        key = name.strip().lower()
        if key in HEAVY_IDS:
            picked["heavy"] = seq
        elif key in LIGHT_IDS:
            picked["light"] = seq
        elif key in ANTIGEN_IDS:
            picked["antigen"] = seq
    if not picked:  # unrecognised ids -> positional: 3 records = H,L,A ; 2 records = H,A
        order = ["heavy", "light", "antigen"] if len(recs) >= 3 else ["heavy", "antigen"]
        picked = dict(zip(order, list(recs.values())))
        log("FASTA ids not recognised; mapped positionally: "
            + ", ".join(f"{k}={n}" for k, n in zip(picked, recs)))
    return {k: _clean(v, k) for k, v in picked.items()}


def load_target_dir(target):
    """A target directory -> {'heavy': .., 'light': .., 'antigen': ..} (missing keys allowed).

    Accepts either heavy.txt / light.txt / antigen.txt or a single target.fasta.
    """
    target = Path(target)
    if not target.is_dir():
        sys.exit(f"target {target} is not a directory")
    fasta = next((target / n for n in ("target.fasta", "sequences.fasta")
                  if (target / n).exists()), None)
    if fasta:
        chains = chains_from_fasta(fasta)
    else:
        chains = {name: read_seq(str(target / f"{name}.txt"), name)
                  for name in ("heavy", "light", "antigen")
                  if (target / f"{name}.txt").exists()}
    if not chains:
        sys.exit(f"target {target} has no heavy.txt/light.txt/antigen.txt and no target.fasta")
    return chains


def resolve_inputs(args):
    """-> (heavy, light, antigen, tag, out_dir, ref) ; any chain may be None except antigen."""
    heavy = light = antigen = None
    target = Path(args.target).resolve() if args.target else None

    if target:
        c = load_target_dir(target)
        heavy, light, antigen = c.get("heavy"), c.get("light"), c.get("antigen")
    elif args.fasta:
        c = chains_from_fasta(args.fasta)
        heavy, light, antigen = c.get("heavy"), c.get("light"), c.get("antigen")

    # explicit flags always win
    heavy = read_seq(args.heavy, "heavy") or heavy
    light = read_seq(args.light, "light") or light
    antigen = read_seq(args.antigen, "antigen") or antigen

    if antigen is None:
        sys.exit("No antigen sequence. Give --target/--fasta, or --antigen.")
    if heavy is None and light is None:
        sys.exit("No antibody chain. Give --heavy and/or --light (light is optional for VHH).")

    tag = args.tag or (target.name if target else "run")
    out_dir = Path(args.out_dir) if args.out_dir else ((target / "runs") if target else HERE / "runs")
    ref = args.ref
    if ref is None and target and (target / "reference.cif").exists():
        ref = str(target / "reference.cif")
    return heavy, light, antigen, tag, out_dir, ref


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    g = ap.add_argument_group("input")
    g.add_argument("--target", help="Directory with heavy.txt/light.txt/antigen.txt (+ reference.cif)")
    g.add_argument("--fasta", help="FASTA with records H / L / A")
    g.add_argument("--heavy", help="Heavy (VH) sequence, or a file holding it")
    g.add_argument("--light", help="Light (VL) sequence, or a file holding it. Omit for VHH/nanobody.")
    g.add_argument("--antigen", help="Antigen sequence, or a file holding it")

    g = ap.add_argument_group("output")
    g.add_argument("--out_dir", default=None, help="Default: <target>/runs, else ./runs")
    g.add_argument("--tag", default=None, help="Prefix for output files (default: target name)")
    g.add_argument("--no_save_all", action="store_true",
                   help="Do not save per-seed CIF/PAE/pLDDT. Disables all scoring — don't.")

    g = ap.add_argument_group("recipe (defaults = the validated recipe; change at your own risk)")
    g.add_argument("--model", default="biohub/ESMFold2-Fast",
                   help="biohub/ESMFold2-Fast (recommended for Ab-Ag) | biohub/ESMFold2 (48 layers)")
    g.add_argument("--num_loops", type=int, default=20)
    g.add_argument("--num_sampling_steps", type=int, default=100)
    g.add_argument("--num_diffusion_samples", type=int, default=1)
    g.add_argument("--num_seeds", type=int, default=25, help="Fold seeds 0..N-1 (default 25)")
    g.add_argument("--seeds", type=int, nargs="+", default=None,
                   help="Explicit seed list; overrides --num_seeds")
    g.add_argument("--lm_dropout", type=float, default=0.3)
    g.add_argument("--dtype", default="fp32", choices=["fp32", "bf16", "fp16"])
    g.add_argument("--lm_mask_pct", type=float, default=None,
                   help="Fraction of residues masked before the PLM (extra diversity knob).")
    g.add_argument("--noise_scale", type=float, default=None)
    g.add_argument("--step_scale", type=float, default=None)

    g = ap.add_argument_group("hardware")
    g.add_argument("--gpu", type=int, default=0, help="GPU index (one complex per GPU)")
    g.add_argument("--shard", action="store_true",
                   help="Shard the model across visible GPUs (device_map=auto) for large complexes")
    g.add_argument("--visible", default=None, help="CUDA_VISIBLE_DEVICES override, e.g. '2,3'")

    g = ap.add_argument_group("optional conditioning (off by default — none improved Ab-Ag accuracy)")
    g.add_argument("--heavy_msa", help="a3m MSA for the heavy chain")
    g.add_argument("--light_msa", help="a3m MSA for the light chain")
    g.add_argument("--antigen_msa", help="a3m MSA for the antigen chain")
    g.add_argument("--msa_max_depth", type=int, default=1024)
    g.add_argument("--antigen_copies", type=int, default=1,
                   help="Fold the antigen as a homo-oligomer (chains A0..An). Scoring assumes 1.")

    g = ap.add_argument_group("scoring")
    g.add_argument("--ref", default=None,
                   help="Reference structure for DockQ (default: <target>/reference.cif if present)")
    g.add_argument("--mapping", default=None,
                   help="DockQ chain mapping, e.g. HLA:HLC when scoring against a native PDB")
    g.add_argument("--no_score", action="store_true", help="Skip the automatic scoring step")

    args = ap.parse_args()
    heavy, light, antigen, tag, out_dir, ref = resolve_inputs(args)
    seeds = args.seeds if args.seeds is not None else list(range(args.num_seeds))
    save_all = not args.no_save_all

    os.environ["CUDA_VISIBLE_DEVICES"] = (
        (args.visible or "0,1,2,3") if args.shard else str(args.gpu))
    if "HF_HOME" not in os.environ:  # env.sh normally sets this; fall back to a repo-local cache
        for base in (HERE, *HERE.parents):
            if (base / "hf_cache" / "hub").is_dir():
                os.environ["HF_HOME"] = str(base / "hf_cache")
                break

    import numpy as np
    import torch
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model

    from esm.models.esmfold2 import (
        ESMFold2InputBuilder,
        ProteinInput,
        StructurePredictionInput,
    )

    dtype = {"bf16": torch.bfloat16, "fp32": torch.float32, "fp16": torch.float16}[args.dtype]
    out_dir.mkdir(parents=True, exist_ok=True)
    results_csv = out_dir / "results.csv"

    log(f"Loading {args.model} (dtype={args.dtype}) on GPU {args.gpu} ...")
    t0 = time.time()
    # The HF wrapper rejects 'sdpa' for ESMFold2Model; 'eager' satisfies it. The internal
    # esmc/esmfold2 attention independently uses flash_attn or F.scaled_dot_product_attention.
    if args.shard:
        model = ESMFold2Model.from_pretrained(
            args.model, torch_dtype=dtype, attn_implementation="eager", device_map="auto").eval()
        log(f"Model sharded: {getattr(model, 'hf_device_map', 'n/a')}")
    else:
        model = ESMFold2Model.from_pretrained(
            args.model, torch_dtype=dtype, attn_implementation="eager").cuda().eval()
    log(f"Model loaded in {time.time() - t0:.1f}s")

    builder = ESMFold2InputBuilder()

    def load_msa(path, name):
        if not path:
            return None
        from esm.models.esmfold2 import MSA
        m = MSA.from_a3m(path, remove_insertions=True)
        log(f"{name} MSA loaded: depth={m.depth} len={m.seqlen}")
        return m

    chains, lengths = [], {}
    if heavy:
        chains.append(ProteinInput(id="H", sequence=heavy, msa=load_msa(args.heavy_msa, "Heavy")))
        lengths["H"] = len(heavy)
    if light:
        chains.append(ProteinInput(id="L", sequence=light, msa=load_msa(args.light_msa, "Light")))
        lengths["L"] = len(light)
    antigen_msa = load_msa(args.antigen_msa, "Antigen")
    if args.antigen_copies > 1:  # homo-oligomer: one ProteinInput carrying N chain ids
        antigen_ids = [f"A{i}" for i in range(args.antigen_copies)]
        chains.append(ProteinInput(id=antigen_ids, sequence=antigen, msa=antigen_msa))
    else:
        antigen_ids = ["A"]
        chains.append(ProteinInput(id="A", sequence=antigen, msa=antigen_msa))
    for cid in antigen_ids:
        lengths[cid] = len(antigen)

    chain_order = [c for c in ("H", "L") if c in lengths] + antigen_ids
    total = sum(lengths.values())
    log(f"Chains: {'  '.join(f'{c}={lengths[c]}' for c in chain_order)}  (total {total} residues)")
    if total > 768:
        log("WARNING: >768 residues — trim constant domains / antigen tags if this OOMs.")

    spi = StructurePredictionInput(sequences=chains)
    ab_idx = {c: chain_order.index(c) for c in ("H", "L") if c in lengths}
    ag_idx = chain_order.index(antigen_ids[0])

    for seed in seeds:
        t0 = time.time()
        log(f"Folding seed={seed} loops={args.num_loops} steps={args.num_sampling_steps} "
            f"samples={args.num_diffusion_samples} lm_dropout={args.lm_dropout}")
        kwargs = dict(
            num_loops=args.num_loops,
            num_sampling_steps=args.num_sampling_steps,
            num_diffusion_samples=args.num_diffusion_samples,
            seed=seed,
            lm_dropout=args.lm_dropout,
            msa_max_depth=args.msa_max_depth,
        )
        for name, val in (("lm_mask_pct", args.lm_mask_pct), ("noise_scale", args.noise_scale),
                          ("step_scale", args.step_scale)):
            if val is not None:
                kwargs[name] = val

        result = builder.fold(model, spi, **kwargs)
        samples = result if isinstance(result, list) else [result]
        dt = time.time() - t0

        for si, r in enumerate(samples):
            iptm = float(r.iptm) if r.iptm is not None else float("nan")
            ptm = float(r.ptm) if r.ptm is not None else float("nan")
            plddt = float(r.plddt.mean())
            # Per-chain-pair interface ipTM (antibody chain -> antigen), parsed later by score.py.
            detail, pci = "n/a", getattr(r, "pair_chains_iptm", None)
            if pci is not None:
                try:
                    m = pci.float()
                    detail = " ".join(f"{c}-A={float(m[i, ag_idx]):.3f}" for c, i in ab_idx.items())
                except Exception as e:  # never let a metric break the run
                    detail = f"(pair_chains_iptm err {e})"
            log(f"  seed={seed} sample={si}: ipTM={iptm:.4f} pTM={ptm:.4f} "
                f"pLDDT={plddt:.2f} AbAg[{detail}]  ({dt:.1f}s)")

            row = dict(tag=tag, seed=seed, sample=si, iptm=iptm, ptm=ptm, plddt=plddt,
                       ab_ag_detail=detail, model=args.model, num_loops=args.num_loops,
                       num_sampling_steps=args.num_sampling_steps,
                       num_diffusion_samples=args.num_diffusion_samples,
                       lm_dropout=args.lm_dropout, dtype=args.dtype, seconds=round(dt, 1))
            write_header = not results_csv.exists()
            with open(results_csv, "a", newline="") as f:
                w = csv.DictWriter(f, fieldnames=list(row.keys()))
                if write_header:
                    w.writeheader()
                w.writerow(row)

            if not save_all:
                continue
            stem = f"{tag}_seed{seed}_s{si}_iptm{iptm:.3f}"
            (out_dir / f"{stem}.cif").write_text(r.complex.to_mmcif())
            # Raw PAE + per-residue pLDDT make the whole score suite (ipSAE, pDockQ,
            # pDockQ2, LIS) reproducible after the fact. Filenames follow ipsae.py's
            # Boltz convention.
            try:
                np.savez_compressed(out_dir / f"{stem}_plddt.npz",
                                    plddt=r.plddt.float().cpu().numpy())
                if getattr(r, "pae", None) is not None:
                    np.savez_compressed(out_dir / f"{stem}_pae.npz",
                                        pae=r.pae.float().cpu().numpy())
                (out_dir / f"{stem}_meta.json").write_text(json.dumps({
                    "stem": stem, "tag": tag, "seed": seed, "sample": si, "model": args.model,
                    "iptm": iptm, "ptm": ptm, "plddt_mean": plddt, "ab_ag_detail": detail,
                    "chain_order": chain_order, "chain_lengths": lengths,
                    "recipe": {k: v for k, v in kwargs.items() if k != "seed"},
                    "dtype": args.dtype,
                }, indent=2))
            except Exception as e:
                log(f"  WARN: failed to save PAE/pLDDT for {stem}: {e}")

        del result, samples
        torch.cuda.empty_cache()

    log(f"Folded {len(seeds)} seeds -> {out_dir}")

    # Score every saved prediction: ipSAE / pDockQ / pDockQ2 / LIS from the saved PAE,
    # plus DockQ vs the reference -> scores.csv, ranked by abag_ipsae. Never fatal:
    # the same numbers can be regenerated later with `python score.py <out_dir>`.
    if save_all and not args.no_score:
        try:
            import subprocess
            score_py = next((b / "score.py" for b in (HERE, *HERE.parents)
                             if (b / "score.py").exists() and (b / "tools" / "ipsae.py").exists()), None)
            if score_py is None:
                log("Auto-score skipped: score.py / tools/ipsae.py not found.")
            else:
                cmd = [sys.executable, str(score_py), str(out_dir)]
                if ref:
                    cmd += ["--ref", str(ref)]
                if args.mapping:
                    cmd += ["--mapping", args.mapping]
                log(f"Scoring -> {out_dir}/scores.csv")
                subprocess.run(cmd, check=False)
        except Exception as e:
            log(f"Auto-score failed (non-fatal, rerun `python score.py {out_dir}`): {e}")


if __name__ == "__main__":
    main()
