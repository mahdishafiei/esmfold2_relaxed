# Targets

A **target** is a directory describing one antibody–antigen complex. Point `fold.py` at it
and everything else (output location, DockQ reference, mutants) follows from it.

```
targets/my_ab/
├── heavy.txt        VH sequence, one line, plain text        (required*)
├── light.txt        VL sequence                              (optional — omit for VHH/nanobody)
├── antigen.txt      antigen sequence                         (required)
├── reference.cif    DockQ reference, if you have one         (optional)
├── runs/            created by fold.py: per-seed CIF/PAE/meta + scores.csv
└── mutants/         created by make_mutant.py
```
\* at least one of `heavy.txt` / `light.txt`.

A single `target.fasta` with records `H`, `L`, `A` works instead of the three `.txt` files.

```bash
mkdir -p targets/my_ab
printf 'EVQLVESGGG...\n' > targets/my_ab/heavy.txt
printf 'DIVMTQSPDS...\n' > targets/my_ab/light.txt
printf 'DQICIGYHAN...\n' > targets/my_ab/antigen.txt

bash predict.sh --target targets/my_ab --gpu 0
```

## Sequence rules

- **Fv only** — VH + VL. Constant domains add residues without helping the interface.
- **Clean antigen** — the biological ectodomain/monomer. Strip expression and trimerization
  tags: on 8UME the tagged 547-aa construct and the HA trimer both mis-docked.
- **≤ 768 residues total** across all chains, or you risk OOM on a 48 GB card.
- Mutation positions used by `make_mutant.py` are 1-based on exactly these sequences.

## reference.cif

Optional, and only used for the DockQ columns of `scores.csv`:

- **A deposited native structure** → DockQ is ground truth (≥0.23 correct, ≥0.49 medium,
  ≥0.80 high). Chain ids rarely match, so pass e.g. `--mapping HLA:HLC`.
- **Your WT prediction** (as in `targets/8ume`) → DockQ measures how far a mutant's dock moved
  from the parent: ≈1 means unchanged, low means the mutation shifted the pose.
- **Nothing** → DockQ columns are skipped; rank by `abag_ipsae` alone.
