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
| Model | ESMFold2-Fast | 24-layer model — same as Biohub Colab default |
| Seeds | 25 | Best by ipTM is saved |
| Loops | 20 | Paper recommendation for ab-ag |
| Diffusion steps | 100 | |
| lm_dropout | 0.3 | Drives conformation diversity across seeds |

---

## Score interpretation

| ipTM | Meaning |
|---|---|
| > 0.8 | Confident — the predicted binding pose is reliable |
| 0.5–0.8 | Plausible — inspect carefully, check CDR contacts |
| < 0.5 | Low confidence — run more seeds or verify epitope constraints |

**ipTM** (interface predicted TM-score) is the key metric for complex quality.  
**pTM** reflects individual chain fold confidence.  
**pLDDT** > 0.7 indicates well-folded regions.

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
