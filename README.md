# esmfold2_relaxed

**ESMFold2 antibody–antigen complex prediction, preconfigured with settings that work.**

A fork of [EvolutionaryScale's `esm`](https://github.com/Biohub/esm) that ships a ready-to-run
pipeline: give it a heavy chain, a light chain and an antigen, and it folds 25 seeds, scores
every one, and hands you a ranked structure. Local GPU only — no API key, no credits.

```bash
git clone git@github.com:mahdishafiei/esmfold2_relaxed.git && cd esmfold2_relaxed
bash setup.sh                                    # once per machine (~26 GB weights)

mkdir -p targets/my_ab
printf 'EVQLVESGGG...\n' > targets/my_ab/heavy.txt      # VH
printf 'DIVMTQSPDS...\n' > targets/my_ab/light.txt      # VL (omit for a VHH/nanobody)
printf 'MKTAYIAKQR...\n' > targets/my_ab/antigen.txt    # clean antigen, no tags

bash predict.sh --target targets/my_ab --gpu 0   # ~48 min for 25 seeds on one L40S
column -s, -t targets/my_ab/runs/scores.csv | less -S
```

---

## Why the defaults are what they are

1. **ESMFold2-Fast, single sequence.** For antibody–antigen the Fast model finds interfaces the
   full 48-layer model misses; MSAs, trimeric antigens and epitope constraints are off.
2. **Don't rank by ipTM.** ipTM does not track epitope correctness for Ab–Ag: a wrong-site pose
   can score higher than a correct one. The pipeline never picks a winner by ipTM.
3. **Scan seeds.** Only a fraction of seeds dock correctly, so the recipe folds 25, keeps every
   structure, and ranks them. One seed is a coin flip, and *which* seeds win changes with the
   GPU — a good seed from another machine's run is worth nothing on yours.

---

## Install

**Needs:** Linux, an NVIDIA GPU with ≥ 24 GB (one complex per GPU), ~30 GB disk, Python 3.12
(`setup.sh` fetches its own via `uv` if you don't have it).

```bash
bash setup.sh     # fold venv + DockQ venv + weights; idempotent, re-run any time
bash check.sh     # says READY, or exactly what is missing
```

| | Where | Why |
|---|---|---|
| fold venv | `$HOME/.ef2_venv` | torch + this repo's `esm` |
| DockQ venv | `$HOME/.ef2_dockq_venv` | separate on purpose — DockQ pins `numpy<2`, `esm` needs `numpy≥2` |
| weights | `$HOME/.ef2_hf_cache` | `ESMFold2-Fast` + `ESMC-6B` + the CCD dictionary, ~26 GB, once |

Override paths with `EF2_VENV_DIR`, `EF2_DOCKQ_VENV`, `EF2_HF_HOME`. On a cluster, keep weights
on shared storage and build one venv per node — a venv on network storage is slow to import.

`predict.sh` activates everything for you; for interactive work `source env.sh` once per shell.

---

## Predict

A **target** is a directory of sequences (see [`targets/README.md`](targets/README.md)). Or skip
the directory:

```bash
bash predict.sh --heavy EVQ... --light DIV... --antigen MKT... --tag my_ab --gpu 0
bash predict.sh --fasta my_ab.fasta --gpu 1        # records named H / L / A
```

One complex per GPU, so four GPUs run four targets in parallel. Use `tmux` — a 25-seed run takes
about 48 minutes.

### The parameters

| Parameter | Default | | Parameter | Default |
|---|---|---|---|---|
| `--model` | `biohub/ESMFold2-Fast` | | `--lm_dropout` | `0.3` |
| `--num_loops` | `20` | | `--dtype` | `fp32` |
| `--num_sampling_steps` | `100` | | `--num_seeds` | `25` (seeds 0–24) |
| `--num_diffusion_samples` | `1` | | MSA / template / pocket | none |

Everything else is in `python fold.py --help`. The optional knobs — `--heavy_msa` /
`--light_msa` / `--antigen_msa`, `--antigen_copies` (fold the antigen as an oligomer), `--full`
model, more loops, `--lm_mask_pct`, `--shard` across GPUs — are off because they did not improve
Ab–Ag accuracy in testing, and several cause OOM. Turn them on deliberately, not by default.

### Output

```
targets/my_ab/
├── BEST_consensus11of25_ipsae0.525_seed6_my_ab_seed6_s0_iptm0.815.cif   ← the answer
└── runs/
    ├── scores.csv          every seed, best first
    ├── results.csv         fold-time metrics and timings
    ├── *_seed<N>_*.cif     all 25 structures
    └── *_pae.npz, *_plddt.npz, *_meta.json    so any score can be recomputed
```

---

## Reading `scores.csv`

| Column | Meaning |
|---|---|
| **`consensus_n`** | how many *other* seeds put the antibody on the same epitope — **the ranking metric** |
| `epitope_n` | antigen residues within 5 Å of the antibody in this pose |
| `iptm`, `ptm`, `plddt_mean` | ESMFold2's own global confidences |
| `HA_iptm`, `LA_iptm` | interface ipTM, heavy→antigen and light→antigen |
| **`abag_ipsae`** | `max(HA_ipsae, LA_ipsae)` — the tie-break *within* a consensus cluster |
| `HA_ipsae`, `LA_ipsae`, `HL_ipsae` | ipSAE per chain pair (`HL` = the Fv's internal packing) |
| `HA_pdockq`, `HA_pdockq2`, `HA_lis` (+ `LA_`, `HL_`) | pDockQ, pDockQ2, LIS per pair |
| `HA_dockq_vs_ref`, `LA_dockq_vs_ref`, `abag_dockq_vs_ref` | DockQ vs `reference.cif`, if you supply one |
| `HA_irmsd`, `HA_lrmsd`, `HA_fnat` (+ `LA_`) | interface RMSD, ligand RMSD, native-contact fraction |

- **`consensus_n`** — a pose many independent seeds reproduce is real; a confident singleton
  usually isn't. With 25 seeds a winning cluster is typically 8–12. This needs no reference
  structure, which is why it ranks.
- **ipSAE** (0–1) — interface-restricted, stricter than ipTM. A genuine interface reads ~0.3–0.7
  and no interface ~0, but it both over- and under-calls, so it orders poses *inside* a cluster
  rather than choosing between clusters. `HL_ipsae` is ~0.8 even for a nonsense dock — that's
  just the Fv folding.
- **ipTM** (0–1) — informative, not decisive.
- **pDockQ / pDockQ2** (0–1) — interface-quality predictors; **> ~0.23** suggests a real interface.
- **LIS** (0–1) — local interaction score; higher = more confident contacts.
- **DockQ vs ref** (0–1) — against a deposited native it is ground truth (**≥ 0.23** correct,
  **≥ 0.49** medium, **≥ 0.80** high). Against your own WT prediction it measures how far a
  mutant's pose moved: ≈1 unchanged, low + large `iRMSD` = the dock shifted.

Re-score any run directory at will — it reads the saved PAE, so nothing is re-folded:

```bash
source env.sh
python score.py <run_dir> [--ref structure.cif] [--mapping HLA:HLC] [--rank ipsae] [--no_dockq]
```

### Picking the answer

1. Take the top row of `scores.csv`; that structure is copied out as `BEST_*.cif`. If ipSAE's
   favourite lost to a bigger cluster, `score.py` says so explicitly.
2. Sanity-check: a real dock sits in a large cluster with elevated `HA_iptm`/`LA_iptm`. If the
   biggest cluster is 1–2 seeds, or everything is near zero, don't trust it — fold more seeds.
3. Look at the top few poses in PyMOL before committing to one.

---

## Mutants

```bash
source env.sh
python make_mutant.py --target targets/my_ab "H:Y57F+L:Q27Y"   # verifies the WT residue first
bash targets/my_ab/mutants/mut_H-Y57F_L-Q27Y/run.sh 0          # last arg = GPU index
python best_overall.py targets/my_ab                           # rank WT + every mutant
```

Positions are 1-based on the target's own chain sequences; a WT mismatch aborts rather than
folding the wrong protein. If the parent target has a `reference.cif`, each mutant links to it,
so `abag_dockq_vs_ref` answers *"did this mutation move the dock?"*.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `no fold venv on this machine` | `bash setup.sh` on *that* machine — venvs are node-local by design |
| `weights not found` | `bash check.sh`; point `EF2_HF_HOME` at an existing cache, or re-run `setup.sh` |
| No DockQ columns | DockQ venv missing (`setup.sh`) or no `--ref` given; every other column is unaffected |
| Scoring didn't run after folding | it can never kill a fold; just run `python score.py <run_dir>` |
| CUDA OOM | one complex per GPU (~22 GB); keep total residues ≤ 768 (Fv + clean antigen); or `--shard` |
| `transformer_engine` / `flash-attn` / `xformers` not installed | harmless — pure-PyTorch fallback, correct results, slightly slower |
| Job dies when your laptop sleeps | run it under `tmux` |

---

## Layout

| Path | What |
|---|---|
| `setup.sh` / `env.sh` / `check.sh` | build the environment once · activate it · verify it |
| `predict.sh` | fold + score in one command |
| `fold.py` | fold driver; recipe as defaults, saves every seed + PAE |
| `score.py` | consensus + ipSAE / pDockQ / pDockQ2 / LIS + DockQ → `scores.csv`, `BEST_*.cif` |
| `make_mutant.py` / `best_overall.py` | scaffold a mutant · rank WT + all mutants |
| `targets/` | one directory per complex (yours; gitignored outputs) |
| `tools/ipsae.py` | vendored official Dunbrack ipSAE scorer |
| `esm/`, `pyproject.toml` | upstream EvolutionaryScale `esm`, unmodified (commit `ba4d712`) |

---

## References

- **ESMFold2 / ESMC** — `biohub/ESMFold2-Fast`, `biohub/ESMFold2`, `biohub/ESMC-6B`
- **ipSAE, pDockQ, pDockQ2, LIS** — see [`tools/README.md`](tools/README.md)
- **DockQ** — Basu & Wallner — https://github.com/bjornwallner/DockQ

```bibtex
@article{candido2025language,
  title  = {Language Modeling Materializes a World Model of Protein Biology},
  author = {Candido, Salvatore and Hayes, Thomas and Derry, Alexander and others},
  journal= {EvolutionaryScale},
  year   = {2025}
}
```
