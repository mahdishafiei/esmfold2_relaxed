# ESMFold2 — Antibody-Antigen Complex Prediction

Fork of [Biohub/esm](https://github.com/Biohub/esm) specialized for predicting antibody-antigen complex structures. Uses **ESMFold2-Fast** — the same model used in the Biohub Colab, optimized for antibody-antigen interactions.

Two modes: **local GPU** (no credits, unlimited) or **Biohub API via Colab** (no GPU needed).

---

## What this does

Takes heavy chain (VH), light chain (VL), and antigen sequences → predicts the 3D complex structure → saves ranked `.cif` files you can open in PyMOL or ChimeraX.

Runs 25 independent seeds and picks the best by **ipTM** (interface confidence score). Each seed gets a different random diffusion initialization, so you get a genuine ensemble.

---

## Install (local GPU)

**Requirements:** Linux · **Python 3.12 exactly** · CUDA GPU ≥16 GB VRAM · ~27 GB disk

```bash
git clone git@github.com:mahdishafiei/esmfold2_relaxed.git
cd esmfold2_relaxed/antibody_complex
bash setup.sh
```

Verify setup works before running real jobs:
```bash
./fold.sh --test
```

`setup.sh` creates a virtualenv, installs the package, installs `xformers` (required for numerical accuracy), and downloads the model weights (~26 GB) from HuggingFace automatically.

---

## Run (local)

```bash
cd antibody_complex

./fold.sh \
  --heavy   "EVQLVESGGGLVKPGGSLRL..." \
  --light   "DIVMTQSPDSLAVSLGERAT..." \
  --antigen "DQICIGYHANNSTEQVDTIM..."
```

Any chain is optional — run with just `--heavy` + `--antigen` if you have no light chain.

**Output** is saved to a timestamped folder:
```
output/20260609_143022/
├── best_seed07_ipTM0.818.cif   ← open this one
├── summary.csv                  ← all 25 seeds ranked by ipTM
└── all_seeds/                   ← every seed CIF, ranked
```

---

## Run on Colab (no GPU, requires Biohub API token)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/mahdishafiei/esmfold2_relaxed/blob/main/antibody_complex/ESMFold2_Antibody_Antigen_Colab.ipynb)

Get a token at [biohub.ai](https://biohub.ai). Free tier is 100 credits/day (~100 single predictions). For high-throughput use, the local GPU option is better.

---

## All options

```
--heavy        VH amino acid sequence
--light        VL amino acid sequence (optional)
--antigen      Antigen amino acid sequence
--seeds        Number of seeds to run (default: 25)
--loops        Recurrent folding loops (default: 20)
--diff-steps   Diffusion sampling steps (default: 100)
--lm-dropout   LM dropout for seed diversity (default: 0.3)
--device       GPU to use, e.g. cuda:1 (default: cuda:0)
--full         Use full ESMFold2 (48 layers) instead of Fast (24 layers)
--contacts-csv Path to a Chai-format epitope constraint CSV
--out          Override output CIF path (default: auto-named)
```

---

## Settings used (and why)

| Setting | Value | Why |
|---|---|---|
| Model | ESMFold2-Fast | Same as Biohub Colab — produces correct structures for ab-ag |
| Seeds | 25 | Paper shows pass rate rises from 49% (1 seed) → ~62% (25 seeds) |
| Loops | 20 | Paper recommendation for antibody-antigen — ~5% better than 10 |
| xformers | installed | Required for numerical equivalence to Biohub servers |
| lm_dropout | 0.3 | Drives conformation diversity; paper inference default |
| lm_mask_pct | 0.0 | Default for full + fast ESMFold2 model |

---

## Sequence length limit

Total residues (VH + VL + antigen) must be **≤ 768**. If over the limit:
- Trim antibody constant regions (keep VH/VL only)
- Remove C-terminal fusion tags from the antigen

---

## Citation

```bibtex
@article{candido2025language,
  title   = {Language Modeling Materializes a World Model of Protein Biology},
  author  = {Candido, Salvatore and Hayes, Thomas and Derry, Alexander and others},
  year    = {2025}
}
```

