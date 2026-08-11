# Validation — 8UME, and where the defaults come from

Every default in this repo comes from a ~1500-run sweep on PDB **8UME** (an anti-influenza
antibody Fv + HA ectodomain), scored against the deposited structure with DockQ. This is the
evidence, including the settings that looked promising and were wrong.

**Result: ESMFold2-Fast reproduces the deposited interface at DockQ 0.936** — correct epitope,
correct fold — from sequence alone, with no MSA, template, or epitope constraint. That
structure is `targets/8ume/reference.cif` (seed 13).

---

## 1. The three rules

**1. Use the Fast model, single-sequence.**
The full `biohub/ESMFold2` (with or without an antigen MSA) and every HA-trimer setup dock this
antibody at the **wrong epitope** — high ipTM, DockQ 0.03–0.09, zero native contacts, over ~179
seeds. Only `biohub/ESMFold2-Fast` finds the right site.

**2. Do not rank by ipTM.**
For this target ipTM is actively misleading. A *wrong* dock reached **ipTM 0.71** (seed 5,
DockQ 0.008) while *correct* docks sat at **ipTM ~0.31** (seeds 1, 7). Rank by structure-aware
scores instead: ipSAE, or DockQ against a reference when you have one.

**3. It is seed-dependent — scan many seeds.**
Of 14 DockQ-scored WT seeds, 7 docked correctly and 5 were high quality (DockQ ≥ 0.80): about
**5 high-quality docks per 25 seeds**. One seed is a coin flip; 25 seeds is the recipe.

---

## 2. The exact winning recipe

| Parameter | Value | | Parameter | Value |
|---|---|---|---|---|
| model | `biohub/ESMFold2-Fast` | | `lm_dropout` | `0.3` |
| `num_loops` | `20` | | dtype | `fp32` |
| `num_sampling_steps` | `100` | | seeds | `0 … 24` |
| `num_diffusion_samples` | `1` | | save all seeds | required (PAE ⇒ ipSAE) |
| antigen copies | 1 (monomer, untagged) | | MSA / template / pocket | none |

Loaded as `ESMFold2Model.from_pretrained("biohub/ESMFold2-Fast", attn_implementation="eager")`.
These are the defaults of `fold.py`, so `python fold.py --target targets/8ume` *is* this recipe.

Validated stack: **torch 2.13.0+cu130, transformers 4.57.6** (the Biohub fork pulled in by
`esm`), numpy 2.5.1, Python 3.12, on NVIDIA L40S (46 GB) and A100 (80 GB). ~115 s/seed on an
L40S ⇒ ~48 min for 25 seeds. One complex per GPU (~22 GB).

---

## 3. Per-seed results (WT, Fast model, DockQ vs deposited 8UME)

| seed | ipTM | heavy–Ag DockQ | fnat | verdict |
|---|---|---|---|---|
| 13 | 0.75 | **0.936** | 1.00 | HIGH ✅ (the reference structure) |
| 19 | 0.81 | 0.915 | 1.00 | HIGH ✅ |
| 2 | 0.80 | 0.904 | 1.00 | HIGH ✅ |
| 3 | 0.76 | 0.876 | 1.00 | HIGH ✅ |
| 6 | 0.75 | 0.840 | 0.93 | HIGH ✅ |
| 1 | 0.36 | 0.732 | 0.93 | medium ✅ |
| 7 | 0.31 | 0.413 | 0.29 | acceptable ✅ |
| 5 | 0.71 | 0.008 | 0.00 | **WRONG** ❌ — the ipTM false positive |
| 0, 4, 8, 9, 10, 24 | 0.26–0.29 | ≤ 0.08 | 0.00 | wrong ❌ |

Scope: 25 seeds folded, these 14 DockQ-scored (triage stopped once enough correct docks were
found). High ipTM *usually* means a correct dock — seed 5 shows it is not a guarantee, and
seeds 1 and 7 show correct docks can look unconfident.

---

## 4. Dead ends — do not repeat these

| Tried | Outcome |
|---|---|
| Full `biohub/ESMFold2`, single-sequence or + antigen MSA | DockQ ~0.03, wrong epitope (~179 seeds) |
| HA **trimer** (3 antigen copies) | overall ipTM ~0.72 but Ab–Ag DockQ ≤ 0.09; OOM under 80 GB |
| More loops (40, 64) | neutral to worse — do not escalate loops chasing ipTM |
| `lm_dropout 0.5` | worse |
| 16 diffusion samples per seed | no better than 1 sample × more seeds; several OOM |
| Per-chain or paired H+L MSA | neutral/worse |
| Full 547-aa tagged antigen construct | worse than the clean 502-aa ectodomain |
| Selecting the winner by ipTM | picks confidently-wrong poses (see seed 5) |

MSA and multimer conditioning are *generally* useful in the ESMFold2 paper (ab–ag DockQ pass
rate 50% → 53% with MSA). They were not useful **here**. The flags remain in `fold.py`; the
defaults reflect what was measured on this target.

---

## 5. Reproducing the validation

```bash
bash predict.sh --target targets/8ume --gpu 0          # 25 seeds, ~48 min on an L40S

mkdir -p validation
curl -o validation/8ume.cif https://files.rcsb.org/download/8UME.cif
python score.py targets/8ume/runs --ref validation/8ume.cif --mapping HLA:HLC
```

`--mapping HLA:HLC` maps the model's H, L, A onto the native's H, L, C. Expect several seeds
at `HA_dockq_vs_ref ≥ 0.80`. Diffusion is seeded, but exact per-seed numbers can drift a little
with GPU model, driver, and dtype — the *distribution* (≈5 high-quality docks per 25 seeds) is
the thing that reproduces, not a specific seed's third decimal.
