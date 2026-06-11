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

`setup.sh` creates a virtualenv, installs the package, builds `xformers` from source to match your torch/CUDA version (~15 min — PyPI wheels don't cover newer CUDA), and downloads the model weights (~26 GB) from HuggingFace automatically.

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
--heavy           VH amino acid sequence
--light           VL amino acid sequence (optional)
--antigen         Antigen sequence
--antigen2        Second antigen chain (e.g. HA2 for a cleaved HA1+HA2 heterodimer)
--seeds           Seeds to run (default: 25)
--loops           Recurrent folding loops (default: 20)
--diff-steps      Diffusion steps (default: 100)
--lm-dropout      LM dropout for seed diversity (default: 0.3)
--device          GPU to use, e.g. cuda:1 (default: cuda:0)
--full            Use full ESMFold2 (48 layers) instead of Fast (24 layers)
--contacts-csv    Path to a Chai-format epitope constraint CSV
--out             Override output CIF path (default: auto-named)
--iptm-target     Auto-escalate loops if best ipTM < this (default: 0.80, set 0 to disable)
--escalate-loops  Loop counts to retry when target unmet (default: 48 64)
--escalate-seeds  Seeds per escalation rung (default: 10)
```

### Auto-escalation

If the best ipTM after the main 25 seeds is below `--iptm-target` (default 0.80), the script automatically retries at higher loop counts (48 → 64) with 10 seeds each, stopping as soon as the target is met. The output `summary.csv` gains a `loops` column so you can see which rung produced the best result.

Disable with `--iptm-target 0` if you know the target is a moderate-confidence binder and want to skip the extra compute.

---

## Settings used (and why)

| Setting | Value | Why |
|---|---|---|
| Model | ESMFold2-Fast | Same as Biohub Colab — best accuracy for ab-ag without MSA |
| Seeds | 25 | Paper: pass rate rises from 49% (1 seed) → ~62% (25 seeds) |
| Loops | 20 | Paper recommendation for ab-ag — ~5% better than 10 loops |
| xformers | built from source | Must match your torch/CUDA version — `setup.sh` builds it automatically. Falls back to PyTorch attention if build fails (results still correct, slightly slower) |
| lm_dropout | 0.3 | Drives conformation diversity across seeds; paper default |
| lm_mask_pct | 0.0 | Default for ESMFold2 and ESMFold2-Fast |
| ipTM target | 0.80 | Auto-escalates to 48 → 64 loops if unmet |

---

## Other ESMFold2 capabilities for antibody-antigen prediction

**Models:**

| Model | Loops | MSA support | Notes |
|---|---|---|---|
| `ESMFold2-Fast` | 24 | No | Used here — best ab-ag accuracy without MSA |
| `ESMFold2` | 48 | Yes (up to 16k seqs) | Better with MSA; use `--full` flag |

**Additional inference options (available via the [Biohub API](https://github.com/Biohub/esm)):**
- **MSA per chain** — provide homologous sequences to improve interface accuracy (ab-ag: 50% → 53% DockQ pass rate)
- **Pocket conditioning** — specify known epitope residues to guide antibody docking (`--contacts-csv`)
- **Distogram conditioning** — provide pairwise distance priors between residues
- **More seeds** — pass rate keeps rising with seed count: 49% at 1 seed → ~62% at 25 → ~65% at 1000
- **More loops** — `--loops 20` vs default 10 gives ~5% better ab-ag pass rate
- **Multi-sample diffusion** — `num_diffusion_samples > 1` returns multiple structures per forward pass

**Benchmarks (DockQ pass rate, from paper):**

| Mode | Pass rate vs AlphaFold3 (47% with MSA) |
|---|---|
| ESMFold2-Fast, no MSA | 50% |
| ESMFold2, no MSA | 50% |
| ESMFold2, with MSA | 53% |
| ESMFold2, 20 loops + MSA | 55% |

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

