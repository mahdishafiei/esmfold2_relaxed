# ESMFold2 — Antibody-Antigen Complex Prediction

Fork of [Biohub/esm](https://github.com/Biohub/esm) specialized for predicting antibody-antigen complex structures using **ESMFold2-Fast**. Runs locally on GPU (no API credits) or via Biohub Forge API on Colab.

> Full documentation → [antibody_complex/README.md](antibody_complex/README.md)

---

## Install

```bash
git clone git@github.com:mahdishafiei/esmfold2_relaxed.git
cd esmfold2_relaxed/antibody_complex
bash setup.sh   # creates venv + downloads weights (~26 GB), run once
```

## Run

```bash
./fold.sh \
  --heavy   "VH_SEQUENCE" \
  --light   "VL_SEQUENCE" \
  --antigen "ANTIGEN_SEQUENCE"
```

Output → `output/YYYYMMDD_HHMMSS/best_seedXX_ipTMX.XXX.cif`

## Or run on Colab (no GPU needed)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mahdishafiei/esmfold2_relaxed/blob/main/antibody_complex/ESMFold2_Antibody_Antigen_Colab.ipynb)

Requires a [Biohub API token](https://biohub.ai).

---

## Key settings

| Parameter | Value |
|---|---|
| Model | ESMFold2-Fast (24-layer) |
| Seeds | 25 — best by ipTM kept |
| Loops | 20 |
| Diffusion steps | 100 |
| Max sequence length | 768 residues total |

## Score guide

| ipTM | |
|---|---|
| > 0.8 | Confident |
| 0.5–0.8 | Plausible, inspect |
| < 0.5 | Low confidence |

---

*Mahdi Shafiei · Scripps Research Institute · [mahdishafiei18@gmail.com](mailto:mahdishafiei18@gmail.com)*
