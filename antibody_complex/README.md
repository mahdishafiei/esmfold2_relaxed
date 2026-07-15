# ESMFold2 — Antibody-Antigen Complex Prediction

Predict the 3D structure of antibody-antigen complexes using **ESMFold2-Fast** — state-of-the-art, fast, and accurate. Two modes: **local GPU** (no credits, no limits) or **Biohub API** (via Colab, no GPU needed).

> Built on [EvolutionaryScale/ESMFold2](https://github.com/Biohub/esm) · Paper: *Language Modeling Materializes a World Model of Protein Biology* (EvolutionaryScale, 2025)

---

## Quick start

### Option A — Local GPU (recommended for repeated use)

**Requirements:** Linux, Python 3.10+, CUDA GPU with ≥16 GB VRAM, ~27 GB disk

```bash
# 1. Clone
git clone git@github.com:mahdishafiei/esmfold2_relaxed.git
cd esmfold2_relaxed/antibody_complex

# 2. One-time setup (creates venv + downloads weights ~26 GB)
bash setup.sh

# 3. Predict
./fold.sh \
  --heavy   "EVQLVESGGGLVKPGGSLRL..." \
  --light   "DIVMTQSPDSLAVSLGERAT..." \
  --antigen "DQICIGYHANNSTEQVDTIM..."
```

Output lands in `output/YYYYMMDD_HHMMSS/`:
```
output/20260609_143022/
├── best_seed07_ipTM0.818.cif   ← open this in PyMOL / ChimeraX
├── summary.csv                  ← all 25 seeds ranked by ipTM
└── all_seeds/                   ← every seed CIF, ranked
```

---

### Option B — Google Colab (no GPU needed, requires Biohub API token)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mahdishafiei/esmfold2_relaxed/blob/main/antibody_complex/ESMFold2_Antibody_Antigen_Colab.ipynb)

1. Click the badge above
2. Enter your [Biohub API token](https://biohub.ai)
3. Paste your sequences
4. Run all cells → downloads `complex.cif`

> **Note:** The Colab version uses the Biohub Forge API (100 credits/day on free tier). Each prediction costs ~1 credit. For high-throughput use, the local GPU option is recommended.

---

## Installation details (local)

`setup.sh` does the following automatically:

| Step | What it does |
|---|---|
| Creates `env/` | Python virtualenv |
| Installs ESM package | From this repo |
| Installs xformers | Required for numerical accuracy |
| Downloads `biohub/ESMFold2-Fast` | ~900 MB — folding model |
| Downloads `biohub/ESMC-6B` | ~25 GB — language model backbone |

---

## Usage

```bash
# Minimal — any chain is optional
./fold.sh --heavy VH_SEQ --antigen AG_SEQ

# With light chain
./fold.sh --heavy VH --light VL --antigen AG

# With epitope constraints (from a Chai CSV file)
./fold.sh --heavy VH --light VL --antigen AG \
          --contacts-csv /path/to/constraints.csv

# Use full ESMFold2 model instead of Fast
./fold.sh --heavy VH --light VL --antigen AG --full

# Custom seeds and output
./fold.sh --heavy VH --light VL --antigen AG --seeds 50
```

### Sequence length limit

Total residues (heavy + light + antigen) must be ≤ 768.  
Trim antibody constant regions and antigen fusion tags if needed.

---

## Default settings

| Parameter | Value | Notes |
|---|---|---|
| Model | ESMFold2-Fast | 24-layer model — same as Biohub Colab default. **Use Fast single-sequence first** — it beat full ESMFold2 + MSA + trimer on the 8UME benchmark (see below). |
| Seeds | 25 | Every seed is saved to `all_seeds/`; **highly seed-dependent for ab-ag** (only ~1/2–1/3 of seeds find the right epitope). |
| Loops | 20 | Paper recommendation for ab-ag. More loops (40/64) did **not** help on 8UME — don't auto-escalate loops chasing ipTM. |
| Diffusion steps | 100 | |
| lm_dropout | 0.3 | Drives conformation diversity across seeds |
| MSA / trimer | off | Single-sequence, monomer antigen. Antigen MSA and HA-trimer both mis-docked 8UME (high ipTM, wrong epitope). |

---

## Score interpretation — read this before trusting ipTM

> **⚠️ ipTM is NOT a reliable arbiter of epitope correctness for antibody–antigen complexes.**
> On the 8UME benchmark (validated against the deposited structure by DockQ), a pose with
> **ipTM 0.71 docked at the completely wrong epitope (DockQ 0.008)**, while **correct docks had
> ipTM as low as 0.31**. Full-model + antigen-MSA and HA-trimer runs reached ipTM 0.65–0.9 but
> DockQ ≤ 0.09 (wrong site). **The single-sequence Fast model at ipTM 0.75 gave DockQ 0.936 (correct).**

**What to actually do:**

1. **If you have a native/reference structure** (benchmarking, close homolog): rank candidates by
   **DockQ vs native**, *not* ipTM. Always run many seeds with every CIF saved, and DockQ each.
   See `../../../../../scratch/esm_fold2_8ume/dockq_eval.py` for a working scorer (needs a
   separate DockQ venv — numpy2/scipy≥1.14, conflicts with the folding env).
2. **If you have NO native structure** (true prospective prediction): you lose the DockQ selector.
   Fold many seeds, treat ipTM only as a weak prior, and **visually inspect the top poses** for
   sane CDR-mediated contacts and epitope consistency across seeds (a pose that recurs across
   independent seeds is more trustworthy than a single high-ipTM outlier).

| ipTM | Loose reading (weak signal only) |
|---|---|
| > 0.8 | *Usually* a correct dock — but not guaranteed (see the 8UME false positive) |
| 0.5–0.8 | Ambiguous — inspect carefully; a correct pose can sit here or lower |
| < 0.5 | Low ipTM does **not** mean wrong — some correct 8UME docks were ~0.31 |

**pTM** reflects individual chain fold confidence. **pLDDT** > 0.7 indicates well-folded regions.

---

## Lessons from the 8UME benchmark (validated)

8UME = an anti-influenza-HA antibody Fv + HA ectodomain. ESMFold2-Fast reproduces the deposited
interface at **DockQ 0.936** (right epitope + right fold). Full writeup and replication guide:
`/home/jovyan/work/scratch/esm_fold2_8ume/` (`README.md`, `SOLUTION_8UME_REPLICATION.md`).

**Winning recipe:** `ESMFold2-Fast`, single-sequence, Fv (VH+VL) + clean HA ectodomain (monomer,
no tag), `--loops 20 --diff-steps 100 --lm-dropout 0.3`, scan 25 seeds, save all, select by DockQ.

**Dead ends — do NOT repeat (all mis-dock this antibody or waste compute):**
- Full `ESMFold2` (`--full`), single-seq or + antigen MSA → DockQ ~0.03, wrong epitope (~179 seeds).
- HA **trimer** (3 antigen copies) → overall ipTM ~0.72 but Ab–Ag DockQ ≤ 0.09; also OOMs <80 GB.
- More loops (40/64), higher lm_dropout, 16 diffusion samples, per-chain/paired H+L MSA,
  full-547 tagged antigen → all neutral/worse, several OOM.
- **Selecting by ipTM instead of DockQ → picks confidently-wrong poses.**

---

## Epitope constraints (Chai CSV format)

If you know the binding epitope, provide constraints to guide docking:

```csv
restraint_id,chainA,res_idxA,chainB,res_idxB,connection_type,confidence,...
restraint0,B,Y32,C,R53,contact,1.0,...
restraint1,A,I103,C,N273,contact,1.0,...
```

- Chain A = heavy, B = light, C = antigen
- Residue indices follow antibody Kabat/IMGT numbering for antibody chains
- Antigen residue indices are sequential (1-indexed)

Pass with `--contacts-csv yourfile.csv`.

---

## Hardware

Tested on 4 × NVIDIA L40S (48 GB). The model fits on a single GPU with ≥16 GB VRAM.  
Use `--device cuda:1` etc. to select a specific GPU.

---

## Citation

```bibtex
@article{candido2025language,
  title   = {Language Modeling Materializes a World Model of Protein Biology},
  author  = {Candido, Salvatore and Hayes, Thomas and Derry, Alexander and others},
  journal = {EvolutionaryScale},
  year    = {2025}
}
```

---

## Contact

Mahdi Shafiei · [mahdishafiei18@gmail.com](mailto:mahdishafiei18@gmail.com)  
Scripps Research Institute
