# esmfold2_relaxed — local usage (this machine)

This install is **ready to run** — no `setup.sh` needed. Instead of downloading its
own 26 GB of weights and rebuilding xformers, it **reuses the working runtime** from
the sibling `../../esmfold2` install, while still running **your** `esm/` source
(the relaxed MSA-subsampling processor).

## How it's wired (`activate.sh`)
- **venv**  → reuses `../../esmfold2/env` (torch 2.12+cu130, transformers 4.57, esm deps)
- **weights** → reuses `../../esmfold2/weights` (ESMFold2-Fast, ESMFold2, ESMC-6B); also symlinked as `./weights`
- **your esm** → `PYTHONPATH` puts this repo's `esm/` ahead of the installed package, so your `processor.py` changes are what actually run
- **isolation** → the pipeline's env is **not** modified; your source only takes priority inside this activated shell

## Run — antibody H/L/antigen complex
```bash
./fold.sh --heavy VH_SEQ --light VL_SEQ --antigen AG_SEQ
./fold.sh --test        # quick 1-seed smoke test
```

## Run — general structure prediction (monomer / any complex / FASTA)
```bash
# single monomer
./fold_general.sh --seq MKTAYIAKQR... --seeds 5

# explicit multi-chain complex
./fold_general.sh --chain A:MKT... --chain B:GEE... --out dimer.cif

# everything in a FASTA folded as one complex (header = chain id)
./fold_general.sh --fasta my_complex.fasta --seeds 10
```
Common flags: `--seeds N` `--loops N` `--diff-steps N` `--device cuda:1` `--full` `--out PATH`.
Ranking is by **ipTM** for complexes, **pLDDT** for monomers.

> ⚠️ **ab-ag caveat (8UME-validated):** ipTM is a *weak* arbiter of epitope correctness — a
> wrong-epitope pose scored ipTM 0.71 while correct docks scored as low as 0.31. Prefer
> **ESMFold2-Fast single-sequence** (not `--full`, no MSA, monomer antigen), scan many seeds,
> and if you have a native structure re-rank `all_seeds/` by **DockQ vs native** rather than ipTM.
> Auto loop-escalation is now **off** by default (more loops didn't help on 8UME). See the
> antibody README "Score interpretation" and `scratch/esm_fold2_8ume/` for the full writeup.

## Output
Each run writes `output/<timestamp>/`:
```
best_*.cif        ← top-ranked structure (open in PyMOL/ChimeraX)
summary.csv       ← all seeds ranked
all_seeds/        ← every seed CIF, ranked
```

## Interactive use
```bash
source activate.sh        # then call python / the scripts directly
```

## Notes
- 4× L40S available; model runs on `cuda:0` by default (`--device cuda:N` to pick another).
- Harmless startup warnings: `transformer_engine`/`flash-attn` not installed → pure-PyTorch
  fallback (correct results, slightly slower).
- `setup.sh` is still there for a clean standalone install elsewhere, but is **not** needed here.
