#!/usr/bin/env python
"""Score every prediction in a run directory and rank them — by seed agreement, not by ipTM.

For each prediction (one ``*_meta.json`` per seed, written by ``fold.py``) this collects
one row of ``scores.csv``:

  Epitope agreement, computed from the structures themselves
    epitope_n, consensus_n                       <- consensus_n is the ranking metric

  ESMFold2 confidences (from ``*_meta.json``)
    iptm, ptm, plddt_mean, HA_iptm, LA_iptm      (HA/LA = per-chain-pair interface ipTM)

  PAE-derived, via the official Dunbrack ``tools/ipsae.py``, for pairs H-A, L-A, H-L
    ipSAE, pDockQ, pDockQ2, LIS
    abag_ipsae = max(HA_ipsae, LA_ipsae)         <- the tie-break within a cluster

  DockQ vs a reference structure (optional, needs --ref and the DockQ binary)
    HA_dockq_vs_ref / _irmsd / _lrmsd / _fnat, same for LA, and abag_dockq_vs_ref
    Reference = the WT prediction -> "did the mutation move the dock?"
    Reference = a deposited native -> ground-truth correctness (>=0.23 correct,
    >=0.49 medium, >=0.80 high); pass --mapping HLA:HLC etc. to match native chain ids.

The winner is copied out as ``BEST_consensus<n>of<total>_ipsae<val>_seed<N>_*.cif``.

Usage:
  python score.py <run_dir> [--ref ref.cif] [--mapping HLA:HLA] [--rank ipsae] [--no_dockq]
"""
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
IPSAE = HERE / "tools" / "ipsae.py"


def default_dockq():
    """DockQ lives in its own venv (it pins numpy<2, which conflicts with esm)."""
    for cand in (os.environ.get("EF2_DOCKQ"),
                 str(Path.home() / ".ef2_dockq_venv" / "bin" / "DockQ"),
                 shutil.which("DockQ")):
        if cand and Path(cand).exists():
            return cand
    return ""


def _pair_key(a, b):
    return frozenset((a, b))


def run_ipsae(pae_npz, cif, pae_cutoff, dist_cutoff, quiet=True):
    """Run the official ipsae.py and parse its .txt output.

    Returns {frozenset({c1, c2}): {ipsae, iptm_d0chn, pdockq, pdockq2, lis}} from the
    'max' (symmetric) row of each chain pair, or {} on failure.
    """
    ps = str(int(pae_cutoff)).zfill(2) if pae_cutoff < 10 else str(int(pae_cutoff))
    ds = str(int(dist_cutoff)).zfill(2) if dist_cutoff < 10 else str(int(dist_cutoff))
    out_txt = Path(str(cif)[:-4] + f"_{ps}_{ds}.txt")
    if not out_txt.exists():
        r = subprocess.run([sys.executable, str(IPSAE), str(pae_npz), str(cif),
                            str(pae_cutoff), str(dist_cutoff)], capture_output=True, text=True)
        if r.returncode != 0 and not quiet:
            sys.stderr.write(f"ipsae failed for {cif}: {r.stderr[-400:]}\n")
    if not out_txt.exists():
        return {}
    out = {}
    with open(out_txt) as fh:
        for line in fh:
            p = line.split()
            if len(p) < 13 or p[0] == "Chn1" or p[4] != "max":
                continue
            try:
                out[_pair_key(p[0], p[1])] = dict(
                    ipsae=float(p[5]), iptm_d0chn=float(p[9]),
                    pdockq=float(p[10]), pdockq2=float(p[11]), lis=float(p[12]))
            except (ValueError, IndexError):
                continue
    return out


def run_dockq(cif, ref, dockq_bin, mapping, allowed_mismatches, quiet=True):
    """DockQ model vs reference -> {'HA': {...}, 'LA': {...}, ...}, or {} on failure."""
    tmp = Path(str(cif)[:-4] + "_dockq_vs_ref.json")
    r = subprocess.run([dockq_bin, str(cif), str(ref), "--mapping", mapping,
                        "--allowed_mismatches", str(allowed_mismatches), "--json", str(tmp)],
                       capture_output=True, text=True)
    if not tmp.exists():
        if not quiet:
            sys.stderr.write(f"DockQ failed for {cif}: {r.stderr[-400:]}\n")
        return {}
    best = json.load(open(tmp)).get("best_result", {})
    return {k: dict(dockq=v.get("DockQ", 0.0), irmsd=v.get("iRMSD", 0.0),
                    lrmsd=v.get("LRMSD", 0.0), fnat=v.get("fnat", 0.0))
            for k, v in best.items()}


def _iptm_from_detail(detail):
    """'H-A=0.587 L-A=0.584' -> {'H': 0.587, 'L': 0.584}."""
    out = {}
    for tok in (detail or "").split():
        if "-A=" in tok:
            chain, val = tok.split("-A=")
            try:
                out[chain] = float(val)
            except ValueError:
                pass
    return out


def score_prediction(meta_path, ref, pae_cutoff, dist_cutoff, dockq_bin, mapping,
                     allowed_mismatches, do_dockq=True, quiet=True):
    meta = json.load(open(meta_path))
    d = Path(meta_path).parent
    stem = meta["stem"]
    cif, pae = d / f"{stem}.cif", d / f"{stem}_pae.npz"
    iptm_pairs = _iptm_from_detail(meta.get("ab_ag_detail", ""))

    row = dict(stem=stem, seed=meta.get("seed"), sample=meta.get("sample"),
               model=meta.get("model"), iptm=meta.get("iptm"), ptm=meta.get("ptm"),
               plddt_mean=round(meta.get("plddt_mean", float("nan")), 4),
               HA_iptm=iptm_pairs.get("H", ""), LA_iptm=iptm_pairs.get("L", ""))

    ip = run_ipsae(pae, cif, pae_cutoff, dist_cutoff, quiet) if pae.exists() else {}
    for tag, key in (("HA", _pair_key("H", "A")), ("LA", _pair_key("L", "A")),
                     ("HL", _pair_key("H", "L"))):
        s = ip.get(key, {})
        row[f"{tag}_ipsae"] = round(s.get("ipsae", 0.0), 4)
        row[f"{tag}_pdockq"] = round(s.get("pdockq", 0.0), 4)
        row[f"{tag}_pdockq2"] = round(s.get("pdockq2", 0.0), 4)
        row[f"{tag}_lis"] = round(s.get("lis", 0.0), 4)
    row["abag_ipsae"] = round(max(row["HA_ipsae"], row["LA_ipsae"]), 4)

    if do_dockq:
        dq = run_dockq(cif, ref, dockq_bin, mapping, allowed_mismatches, quiet)
        for tag in ("HA", "LA"):
            s = dq.get(tag, {})
            row[f"{tag}_dockq_vs_ref"] = round(s.get("dockq", 0.0), 4)
            row[f"{tag}_irmsd"] = round(s.get("irmsd", 0.0), 3)
            row[f"{tag}_lrmsd"] = round(s.get("lrmsd", 0.0), 3)
            row[f"{tag}_fnat"] = round(s.get("fnat", 0.0), 3)
        row["abag_dockq_vs_ref"] = round(max(row["HA_dockq_vs_ref"], row["LA_dockq_vs_ref"]), 4)
    return row


def epitope_of(cif, ab_chains=("H", "L"), ag_chain="A", cutoff=5.0):
    """Antigen residues within `cutoff` Å of any antibody atom — the predicted epitope."""
    import biotite.structure.io.pdbx as pdbx
    import numpy as np
    from scipy.spatial import cKDTree

    arr = pdbx.get_structure(pdbx.CIFFile.read(str(cif)), model=1)
    ab = arr[np.isin(arr.chain_id, list(ab_chains))]
    ag = arr[arr.chain_id == ag_chain]
    if len(ab) == 0 or len(ag) == 0:
        return frozenset()
    hits = cKDTree(ab.coord).query_ball_point(ag.coord, r=cutoff)
    return frozenset(int(r) for r, h in zip(ag.res_id, hits) if h)


def add_consensus(runs, rows, min_jaccard=0.5, quiet=True):
    """How many *other* seeds put the antibody on the same epitope.

    No single confidence score is a reliable arbiter here: on a 25-seed 8UME run the
    top seed by both ipTM (0.825) and ipSAE (0.651) was docked 17.7 Å away at the wrong
    site, while three correct docks scored ipSAE 0.000. What does separate them is
    agreement: a pose 11 independent seeds reproduce is real, a confident singleton is
    not. This needs no reference structure, so it works on novel designs too.
    """
    ag = "A"
    meta = next(iter(sorted(Path(runs).glob("*_meta.json"))), None)
    if meta:
        order = json.load(open(meta)).get("chain_order") or []
        ag = next((c for c in order if c.startswith("A")), "A")
    try:
        epitopes = {r["stem"]: epitope_of(Path(runs) / f"{r['stem']}.cif", ag_chain=ag)
                    for r in rows}
    except Exception as e:
        sys.stderr.write(f"[score] epitope consensus unavailable ({e}); ranking by ipSAE only.\n")
        return False
    for r in rows:
        a = epitopes[r["stem"]]
        r["epitope_n"] = len(a)
        r["consensus_n"] = sum(
            1 for other, b in epitopes.items()
            if other != r["stem"] and a and b and len(a & b) / len(a | b) >= min_jaccard)
    if not quiet:
        top = max(rows, key=lambda r: r["consensus_n"])
        print(f"[score] epitope consensus: largest cluster = {top['consensus_n'] + 1}"
              f"/{len(rows)} seeds")
    return True


def default_mapping(runs):
    """Model chains actually folded -> DockQ mapping (e.g. HLA:HLA, or HA:HA for a VHH)."""
    meta = next(iter(sorted(Path(runs).glob("*_meta.json"))), None)
    if meta:
        order = json.load(open(meta)).get("chain_order") or []
        chains = "".join(c for c in order if c in ("H", "L", "A"))
        if chains:
            return f"{chains}:{chains}"
    return "HLA:HLA"


def write_best(runs, rows):
    """Copy the winner (rows[0] — they arrive already ranked) to BEST_*.cif next to the
    run dir, so it's obvious without opening scores.csv. Stale copies are removed."""
    if not rows:
        return None
    top = rows[0]
    src = Path(runs) / f"{top['stem']}.cif"
    if not src.exists():
        return None
    dest_dir = Path(runs).parent if Path(runs).name == "runs" else Path(runs)
    for old in dest_dir.glob("BEST_*.cif"):
        try:
            old.unlink()
        except OSError:
            pass
    val = float(top.get("abag_ipsae", 0) or 0)
    label = (f"consensus{top['consensus_n']}of{len(rows)}_ipsae{val:.3f}"
             if top.get("consensus_n") not in (None, "") else f"abag_ipsae{val:.3f}")
    dest = dest_dir / f"BEST_{label}_seed{top.get('seed', '?')}_{top['stem']}.cif"
    shutil.copy2(src, dest)
    return dest


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs_dir")
    ap.add_argument("--ref", default=None,
                    help="DockQ reference (default: reference.cif beside the run dir)")
    ap.add_argument("--mapping", default=None,
                    help="DockQ chain mapping model:ref (default: HLA:HLA / HA:HA as folded)")
    ap.add_argument("--allowed_mismatches", type=int, default=20,
                    help="Sequence mismatches DockQ tolerates (point mutants vs the WT ref)")
    ap.add_argument("--rank", choices=["auto", "consensus", "ipsae"], default="auto",
                    help="auto: epitope consensus (ipSAE tie-break) once there are >=5 seeds; "
                         "ipsae: rank by abag_ipsae alone")
    ap.add_argument("--pae_cutoff", type=float, default=10.0)
    ap.add_argument("--dist_cutoff", type=float, default=15.0)
    ap.add_argument("--dockq_bin", default=default_dockq())
    ap.add_argument("--no_dockq", action="store_true")
    ap.add_argument("--out", default=None, help="Output CSV (default: <runs_dir>/scores.csv)")
    ap.add_argument("--best-only", action="store_true",
                    help="Don't rescore; just rebuild BEST_*.cif from an existing scores.csv")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    runs = Path(args.runs_dir)
    csv_path = Path(args.out) if args.out else runs / "scores.csv"

    if args.best_only:
        if not csv_path.exists():
            sys.exit(f"--best-only needs an existing {csv_path}")
        with open(csv_path) as fh:
            best = write_best(runs, list(csv.DictReader(fh)))
        print(f"[score] BEST -> {best}" if best else "[score] no BEST written")
        return

    metas = sorted(runs.glob("*_meta.json"))
    if not metas:
        sys.exit(f"No *_meta.json in {runs} (fold.py writes them unless --no_save_all).")

    ref = args.ref
    if ref is None:  # a target dir keeps its DockQ reference next to runs/
        for cand in (runs.parent / "reference.cif", runs / "reference.cif"):
            if cand.exists():
                ref = str(cand)
                break
    mapping = args.mapping or default_mapping(runs)

    do_dockq = not args.no_dockq
    if do_dockq and not ref:
        sys.stderr.write("[score] no --ref given; skipping DockQ columns.\n")
        do_dockq = False
    elif do_dockq and not Path(ref).exists():
        sys.stderr.write(f"[score] reference not found ({ref}); skipping DockQ.\n")
        do_dockq = False
    if do_dockq and not (args.dockq_bin and Path(args.dockq_bin).exists()):
        sys.stderr.write("[score] DockQ binary not found (build it with setup.sh); skipping DockQ.\n")
        do_dockq = False
    if do_dockq and not args.quiet:
        print(f"[score] DockQ vs {ref}  (mapping {mapping})")

    rows = []
    for m in metas:
        try:
            row = score_prediction(m, ref, args.pae_cutoff, args.dist_cutoff, args.dockq_bin,
                                   mapping, args.allowed_mismatches, do_dockq, args.quiet)
            rows.append(row)
            if not args.quiet:
                print(f"  {row['stem']}: iptm={float(row['iptm']):.3f} "
                      f"abag_ipsae={row['abag_ipsae']:.3f} HA_ipsae={row['HA_ipsae']:.3f}"
                      + (f" abag_dockq_vs_ref={row.get('abag_dockq_vs_ref', 0):.3f}"
                         if do_dockq else ""))
        except Exception as e:
            sys.stderr.write(f"[score] failed on {m}: {e}\n")

    if not rows:
        sys.exit("No rows scored.")

    # Rank by how many seeds agree on the epitope, ipSAE as tie-break — a lone confident
    # pose loses to one that eleven independent seeds reproduce. Needs enough seeds to
    # mean anything; --rank ipsae restores pure ipSAE ordering.
    have_consensus = add_consensus(runs, rows, quiet=args.quiet) if len(rows) > 1 else False
    use_consensus = have_consensus and (args.rank == "consensus"
                                        or (args.rank == "auto" and len(rows) >= 5))
    by_ipsae = sorted(rows, key=lambda r: r.get("abag_ipsae", 0.0), reverse=True)
    if use_consensus:
        rows.sort(key=lambda r: (r["consensus_n"], r.get("abag_ipsae", 0.0)), reverse=True)
        if rows[0]["stem"] != by_ipsae[0]["stem"]:
            print(f"[score] NOTE: top by ipSAE is seed {by_ipsae[0]['seed']} but only "
                  f"{by_ipsae[0]['consensus_n']}/{len(rows) - 1} other seeds agree with its "
                  f"epitope — picked seed {rows[0]['seed']} instead "
                  f"({rows[0]['consensus_n']}/{len(rows) - 1} agree).")
    else:
        rows[:] = by_ipsae
    cols = list(rows[0].keys())
    for r in rows:  # union of keys (DockQ columns may be absent)
        cols += [k for k in r if k not in cols]
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    basis = "epitope consensus, ipSAE tie-break" if use_consensus else "abag_ipsae"
    print(f"[score] wrote {len(rows)} rows -> {csv_path}  (ranked by {basis})")
    best = write_best(runs, rows)
    if best:
        print(f"[score] BEST structure -> {best}")


if __name__ == "__main__":
    main()
