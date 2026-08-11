# 8UME — the validation target

Anti-influenza antibody Fv + hemagglutinin (HA) ectodomain, from PDB
[8UME](https://www.rcsb.org/structure/8UME). This is the complex the recipe in this repo was
tuned and validated on: ESMFold2-Fast reproduces the deposited interface at **DockQ 0.936**
(right epitope, right fold) from sequence alone.

Use it to confirm a fresh install works before trusting the pipeline on your own antibody.

| File | What |
|---|---|
| `heavy.txt` | VH, 120 aa — verified against RCSB |
| `light.txt` | VL, 114 aa |
| `antigen.txt` | HA ectodomain, 502 aa — one protomer, **no trimerization tag** |
| `reference.cif` | the DockQ 0.936 prediction (seed 13). Default DockQ reference for mutants. |

Each Fv binds a single HA protomer (~6 Å epitope), so folding the monomer is enough — the
trimer costs 3× the residues, mis-docks, and OOMs.

## Reproduce it

```bash
bash predict.sh --target targets/8ume --gpu 0        # 25 seeds, ~48 min on one L40S
column -s, -t targets/8ume/runs/scores.csv | less -S
```

Expect several seeds with `abag_dockq_vs_ref ≈ 0.9–1.0` (i.e. matching `reference.cif`) and
`abag_ipsae` around 0.3–0.7 on those seeds. Roughly 5 of 25 seeds land a high-quality dock —
that is the expected hit rate, not a failure.

## Score against the deposited structure instead

```bash
mkdir -p validation && curl -o validation/8ume.cif https://files.rcsb.org/download/8UME.cif
python score.py targets/8ume/runs --ref validation/8ume.cif --mapping HLA:HLC
```

Native 8UME names its chains H, L, C — hence `--mapping HLA:HLC`. DockQ ≥ 0.23 is a correct
epitope, ≥ 0.49 medium, ≥ 0.80 high. Per-seed results from the original sweep are in
[`../../docs/VALIDATION_8UME.md`](../../docs/VALIDATION_8UME.md).
