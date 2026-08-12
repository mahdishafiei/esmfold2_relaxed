#!/usr/bin/env python
"""Rank the best prediction of every run under a directory (WT + all mutants).

Scans <dir>/**/runs/scores.csv, takes each run's top row by abag_ipsae, ranks them, and
writes <dir>/BEST_overall.csv with a best_cif column pointing at each run's BEST_*.cif.

Usage:
  python best_overall.py targets/my_ab          # the WT run and every mutant of it
  python best_overall.py targets/my_ab/mutants  # mutants only
"""
import csv
import sys
from pathlib import Path

COLS = ["run", "seed", "iptm", "abag_ipsae", "HA_ipsae", "LA_ipsae",
        "HA_pdockq2", "HA_lis", "abag_dockq_vs_ref", "HA_irmsd", "best_cif"]


def _f(row, key):
    try:
        return float(row.get(key) or 0)
    except (TypeError, ValueError):
        return 0.0


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
    rows = []
    for scores in sorted(root.glob("**/runs/scores.csv")):
        run_dir = scores.parent.parent
        with open(scores) as fh:
            data = list(csv.DictReader(fh))
        if not data:
            continue
        top = max(data, key=lambda r: _f(r, "abag_ipsae"))
        best_cif = next(iter(sorted(run_dir.glob("BEST_*.cif"))), "")
        rows.append({
            "run": run_dir.name if run_dir != root else ".",
            "seed": top.get("seed"), "iptm": top.get("iptm"),
            "abag_ipsae": top.get("abag_ipsae"), "HA_ipsae": top.get("HA_ipsae"),
            "LA_ipsae": top.get("LA_ipsae"), "HA_pdockq2": top.get("HA_pdockq2"),
            "HA_lis": top.get("HA_lis"),
            "abag_dockq_vs_ref": top.get("abag_dockq_vs_ref"),
            "HA_irmsd": top.get("HA_irmsd"),
            "best_cif": str(best_cif),
        })
    if not rows:
        sys.exit(f"No */runs/scores.csv found under {root}")
    rows.sort(key=lambda r: _f(r, "abag_ipsae"), reverse=True)

    out = root / "BEST_overall.csv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)

    width = max(len(r["run"]) for r in rows) + 2
    print(f"{'run':{width}} {'seed':>4} {'ipTM':>6} {'abag_ipsae':>10} {'DockQ_vs_ref':>12}")
    for r in rows:
        print(f"{r['run']:{width}} {str(r['seed']):>4} {_f(r, 'iptm'):6.3f} "
              f"{_f(r, 'abag_ipsae'):10.3f} {_f(r, 'abag_dockq_vs_ref'):12.3f}")
    best = rows[0]
    print(f"\nBEST: {best['run']} seed {best['seed']} (abag_ipsae {_f(best, 'abag_ipsae'):.3f})")
    print(f"  -> {best['best_cif']}")
    print(f"[best_overall] wrote {out}")


if __name__ == "__main__":
    main()
