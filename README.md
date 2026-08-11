# esmfold2_relaxed

**ESMFold2 antibody–antigen complex prediction, with the settings that actually work.**

A fork of [EvolutionaryScale's `esm`](https://github.com/Biohub/esm) whose defaults are a
*validated recipe* rather than library defaults. On PDB **8UME** it reproduces the deposited
antibody–antigen interface at **DockQ 0.936** — right epitope, right fold — from sequence
alone: no MSA, no template, no epitope constraint.

Clone it on any GPU box, run `setup.sh` once, and one command gives you 25 ranked structures
plus a full score table. Nothing here needs an API key or credits.

```bash
git clone git@github.com:mahdishafiei/esmfold2_relaxed.git && cd esmfold2_relaxed
bash setup.sh                                        # once per machine (~26 GB weights)
bash predict.sh --target targets/8ume --gpu 0        # fold + score, ~48 min on one L40S
column -s, -t targets/8ume/runs/scores.csv | less -S
```

---

## The three rules the defaults encode

Learned from a ~1500-run sweep scored against the deposited structure
([full evidence](docs/VALIDATION_8UME.md)):

1. **Fast model, single sequence.** `biohub/ESMFold2-Fast` finds the right epitope. The full
   `biohub/ESMFold2` — with or without an antigen MSA — and every trimeric-antigen setup dock
   the antibody at the *wrong* site with confident-looking scores.
2. **Never rank by ipTM.** A wrong dock hit **ipTM 0.71**; correct docks sat at **ipTM 0.31**.
   Rank by **ipSAE** (interface-restricted), and by **DockQ** against a reference when you have
   one. `fold.py` deliberately does not pick a winner by ipTM.
3. **Scan seeds.** About **5 of 25 seeds** produce a high-quality dock. Fold all 25, keep every
   structure, score each. One seed is a coin flip.

---

## Install

**Needs:** Linux, an NVIDIA GPU with ≥ 24 GB (one complex per GPU), ~30 GB disk, Python 3.12
(`setup.sh` fetches its own via `uv` if you don't have it).

```bash
bash setup.sh     # fold venv + DockQ venv + weights; idempotent, re-run any time
bash check.sh     # says READY, or exactly what is missing
```

`setup.sh` builds three things and skips whatever already exists:

| | Where | Why |
|---|---|---|
| fold venv | `$HOME/.ef2_venv` | torch 2.13 + this repo's `esm` (the validated stack) |
| DockQ venv | `$HOME/.ef2_dockq_venv` | separate on purpose — DockQ pins `numpy<2`, `esm` needs `numpy≥2` |
| weights | `$HOME/.ef2_hf_cache` | `biohub/ESMFold2-Fast` + `biohub/ESMC-6B`, ~26 GB, downloaded once |

Override any path by exporting `EF2_VENV_DIR`, `EF2_DOCKQ_VENV`, `EF2_HF_HOME` before running.
On a multi-node cluster, keep **weights on shared storage** (point `EF2_HF_HOME` at them) and
build **one venv per node** — a venv on network storage is painfully slow to import.

`predict.sh` activates everything for you. For interactive work, `source env.sh` once per shell
(it exports `$EF2_PY` and `$EF2_DOCKQ`, and sets `HF_HOME`).

---

## Predict

A **target** is a directory with `heavy.txt`, `light.txt`, `antigen.txt` (see
[`targets/README.md`](targets/README.md)):

```bash
mkdir -p targets/my_ab
printf 'EVQLVESGGG...\n' > targets/my_ab/heavy.txt      # VH only
printf 'DIVMTQSPDS...\n' > targets/my_ab/light.txt      # VL only (omit for a VHH)
printf 'DQICIGYHAN...\n' > targets/my_ab/antigen.txt    # clean ectodomain, no tags

bash predict.sh --target targets/my_ab --gpu 0
```

or without a target directory:

```bash
bash predict.sh --heavy EVQ... --light DIV... --antigen DQI... --tag my_ab --gpu 0
bash predict.sh --fasta my_ab.fasta --gpu 1             # records named H / L / A
```

Four GPUs means four complexes in parallel — one per `--gpu`. Run under `tmux` so a dropped
connection doesn't kill a 48-minute job.

### What you get

```
targets/my_ab/
├── BEST_abag_ipsae0.612_seed13_my_ab_seed13_s0_iptm0.745.cif   ← the answer
└── runs/
    ├── scores.csv          every seed, ranked by abag_ipsae (best first)
    ├── results.csv         raw fold-time metrics + timings
    ├── *_seed<N>_*.cif     all 25 structures
    └── *_pae.npz, *_plddt.npz, *_meta.json    saved so any score can be recomputed
```

### The recipe (these are the defaults — changing them is the experiment)

| Parameter | Default | | Parameter | Default |
|---|---|---|---|---|
| `--model` | `biohub/ESMFold2-Fast` | | `--lm_dropout` | `0.3` |
| `--num_loops` | `20` | | `--dtype` | `fp32` |
| `--num_sampling_steps` | `100` | | `--num_seeds` | `25` (seeds 0–24) |
| `--num_diffusion_samples` | `1` | | MSA / template / pocket | none |

`python fold.py --help` lists everything else (MSA files, `--antigen_copies`, `--shard` across
GPUs, `--lm_mask_pct`). Those knobs all made 8UME *worse* — see [dead ends](#dead-ends).

---

## Reading `scores.csv`

Every seed is scored from its saved PAE, so scoring is reproducible and re-runnable without
re-folding.

| Column | Meaning |
|---|---|
| `iptm`, `ptm`, `plddt_mean` | ESMFold2's own global confidences |
| `HA_iptm`, `LA_iptm` | interface ipTM, heavy→antigen and light→antigen |
| **`abag_ipsae`** | `max(HA_ipsae, LA_ipsae)` — **the ranking metric**, rows are sorted by it |
| `HA_ipsae`, `LA_ipsae`, `HL_ipsae` | ipSAE per chain pair (`HL` = the Fv's internal packing) |
| `HA_pdockq`, `HA_pdockq2`, `HA_lis` (+ `LA_`, `HL_`) | pDockQ, pDockQ2, LIS per pair |
| `HA_dockq_vs_ref`, `LA_dockq_vs_ref`, `abag_dockq_vs_ref` | DockQ vs `reference.cif`, if present |
| `HA_irmsd`, `HA_lrmsd`, `HA_fnat` (+ `LA_`) | interface RMSD, ligand RMSD, native-contact fraction |

**How to read them**

- **ipSAE** (0–1) — interface-restricted and stricter than ipTM. A genuine Ab–Ag interface reads
  **~0.3–0.7**; no interface reads ~0. This is the compass. Don't confuse it with `HL_ipsae`,
  which is ~0.8 even when the antibody is docked in completely the wrong place.
- **ipTM** (0–1) — informative, not decisive. Correct docks here ranged 0.31–0.85, and a wrong
  dock reached 0.71.
- **pDockQ / pDockQ2** (0–1) — interface quality predictors; **> ~0.23** suggests a real interface.
- **LIS** (0–1) — local interaction score; higher means more confident local contacts.
- **DockQ vs ref** (0–1) — against a *native* structure it is ground truth (**≥ 0.23** correct,
  **≥ 0.49** medium, **≥ 0.80** high). Against your *WT prediction* it measures how far a mutant's
  pose moved: ≈1 = unchanged, low + large `iRMSD` = the mutation shifted the dock.
- **A wrong dock's fingerprint:** `abag_ipsae ≈ 0`, `abag_dockq_vs_ref ≈ 0.007`, `iRMSD` 15–18 Å,
  while `HL_ipsae ≈ 0.8` — the Fv folds fine, it just binds the wrong place.

Re-score any run directory at will (idempotent — it reads the saved PAE):

```bash
source env.sh
python score.py <run_dir> [--ref structure.cif] [--mapping HLA:HLC] [--no_dockq]
```

### Picking the answer

1. Take the top row of `scores.csv` — it's already sorted by `abag_ipsae`, and that structure is
   copied out as `BEST_*.cif`.
2. Sanity-check it: a real dock also shows elevated `HA_iptm`/`LA_iptm` and `abag_dockq_vs_ref`
   well above the wrong-dock floor (~0.007). If the whole run tops out near zero, either the
   mutation broke binding or that target needs more seeds.
3. If a native structure exists, score against it (`--ref native.cif --mapping HLA:HLC`) and use
   **DockQ vs native** as ground truth. A novel design has no native — `abag_ipsae` plus DockQ
   against the WT prediction is the substitute.

---

## Mutants

```bash
source env.sh
python make_mutant.py --target targets/8ume "H:Y57F+L:Q27Y"   # verifies the WT residue first
bash targets/8ume/mutants/mut_H-Y57F_L-Q27Y/run.sh 0          # last arg = GPU index
python best_overall.py targets/8ume                           # rank WT + every mutant
```

Positions are 1-based on the target's own chain sequences. A WT-residue mismatch aborts, so a
typo can't silently fold the wrong protein. Each mutant folder is a self-contained target whose
`reference.cif` links back to the parent, so its `abag_dockq_vs_ref` answers *"did this mutation
move the dock?"* — DockQ tolerates the substitutions via `--allowed_mismatches`.

---

## Validate your install

```bash
bash predict.sh --target targets/8ume --gpu 0
mkdir -p validation && curl -o validation/8ume.cif https://files.rcsb.org/download/8UME.cif
python score.py targets/8ume/runs --ref validation/8ume.cif --mapping HLA:HLC
```

Expect several seeds at `HA_dockq_vs_ref ≥ 0.80` against the native structure. Full per-seed
results and the reasoning behind every default: [`docs/VALIDATION_8UME.md`](docs/VALIDATION_8UME.md).

---

## Dead ends

Measured on 8UME, all worse than the defaults — don't spend GPU hours re-discovering them:

- Full `biohub/ESMFold2`, with or without antigen MSA → DockQ ~0.03, wrong epitope (~179 seeds).
- Antigen as a trimer → ipTM ~0.72 but Ab–Ag DockQ ≤ 0.09, and OOM under 80 GB.
- More loops (40, 64), `lm_dropout 0.5`, 16 diffusion samples, per-chain or paired H+L MSA, the
  full tagged antigen construct → neutral or worse, several OOM.
- Escalating loops or seeds to chase a target ipTM → chases a metric that doesn't track
  correctness here.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `no fold venv on this machine` | `bash setup.sh` on *that* machine — venvs are node-local by design |
| `weights not found` | `bash check.sh`; point `EF2_HF_HOME` at an existing cache, or re-run `setup.sh` |
| No DockQ columns in `scores.csv` | DockQ venv missing (`setup.sh`) or no `--ref` given — ipSAE columns are unaffected |
| Scoring didn't run after folding | it's wrapped so it can never kill a fold; just run `python score.py <run_dir>` |
| CUDA OOM | one complex per GPU (~22 GB); trim to Fv + clean ectodomain (≤768 aa total); or `--shard` |
| `transformer_engine` / `flash-attn` / `xformers` not installed | harmless — pure-PyTorch fallback, correct results, a little slower |
| Job dies when your laptop sleeps | run it under `tmux` |

---

## Layout

| Path | What |
|---|---|
| `setup.sh` / `env.sh` / `check.sh` | build the environment once · activate it · verify it |
| `predict.sh` | fold + score in one command (wraps `env.sh` + `fold.py`) |
| `fold.py` | the fold driver; validated recipe as defaults, saves every seed + PAE |
| `score.py` | ipSAE / pDockQ / pDockQ2 / LIS + DockQ → `scores.csv`, `BEST_*.cif` |
| `make_mutant.py` / `best_overall.py` | scaffold a mutant from a spec · rank WT + all mutants |
| `targets/` | one directory per complex; `targets/8ume/` is the validation target |
| `tools/ipsae.py` | vendored official Dunbrack ipSAE scorer |
| `docs/VALIDATION_8UME.md` | the evidence: three rules, per-seed DockQ table, dead ends |
| `esm/`, `cookbook/`, `pyproject.toml` | upstream EvolutionaryScale `esm`, unmodified |

The `esm/` package is vendored from upstream at commit `ba4d712` — the exact snapshot the recipe
was validated against — and is otherwise untouched. (The name *relaxed* is historical: this fork
started as a patch relaxing MSA subsampling, which is upstream now.)

---

## References

- **PDB 8UME** — https://www.rcsb.org/structure/8UME
- **ESMFold2 / ESMC** — models `biohub/ESMFold2-Fast`, `biohub/ESMFold2`, `biohub/ESMC-6B`
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

Maintained by Mahdi Shafiei · Scripps Research · [mahdishafiei18@gmail.com](mailto:mahdishafiei18@gmail.com)
