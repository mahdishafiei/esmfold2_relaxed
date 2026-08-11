# tools/

Vendored third-party scoring code.

## `ipsae.py`

The official ipSAE script by Roland Dunbrack (Fox Chase Cancer Center), version 4, MIT
licensed — copied here unmodified so scoring is reproducible without a network fetch.

It computes, per chain pair, from the saved PAE:

| Metric | Reference |
|---|---|
| **ipSAE** | Dunbrack, bioRxiv 2025 — https://www.biorxiv.org/content/10.1101/2025.02.10.637595v2 |
| **pDockQ** | Bryant, Pozzati & Elofsson, Nat. Commun. 2022 — https://www.nature.com/articles/s41467-022-28865-w |
| **pDockQ2** | Zhu, Shenoy, Kundrotas & Elofsson, Bioinformatics 2023 — https://academic.oup.com/bioinformatics/article/39/7/btad424/7219714 |
| **LIS** | Kim et al., bioRxiv 2024 — https://www.biorxiv.org/content/10.1101/2024.02.19.580970v1 |

Upstream: https://github.com/DunbrackLab/IPSAE — `score.py` calls it as
`python tools/ipsae.py <pae.npz> <model.cif> 10 15` (the Boltz-style npz input path)
and parses the `max` row of each chain pair.

> The rest of EvolutionaryScale's original `tools/` moved to Apps under
> [Forge](https://forge.evolutionaryscale.ai/).
