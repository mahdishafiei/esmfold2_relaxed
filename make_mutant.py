#!/usr/bin/env python
"""Scaffold a mutant of a target from a mutation spec, ready to fold with the same recipe.

A spec is chain-tagged point mutations joined by '+' or ',':
    H:Y57F+L:Q27Y        H:Y57F,L:Q27Y        A:N145D
Each mutation is <WT><pos><MUT>, 1-based on that chain's sequence in the parent target.
The WT residue is verified before substituting, so a typo aborts instead of silently
folding the wrong sequence.

Creates <target>/mutants/mut_<...>/ containing:
    heavy.txt / light.txt / antigen.txt   mutated chains (unchanged chains copied as-is)
    reference.cif -> ../../reference.cif  the parent's structure, so DockQ vs WT works
    run.sh                                the validated recipe, one command

Usage:
    python make_mutant.py --target targets/my_ab "H:Y57F+L:Q27Y"
    python make_mutant.py --target targets/my_ab "H:Y57F" --print-only
    bash targets/my_ab/mutants/mut_H-Y57F/run.sh 0        # fold it on GPU 0
"""
import argparse
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).parent.resolve()
sys.path.insert(0, str(HERE))
from fold import load_target_dir  # noqa: E402  (stdlib-only import path — no torch)

CHAINS = {"H": "heavy", "L": "light", "A": "antigen"}
MUT_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

RUN_SH = """#!/bin/bash
# {desc} — validated recipe (ESMFold2-Fast, single-sequence, fp32, 25 seeds),
# then automatic scoring (ipSAE / pDockQ / pDockQ2 / LIS + DockQ vs the parent structure).
# Usage: bash run.sh [GPU]     (default GPU 0)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
EF2_ROOT="${{EF2_ROOT:-{repo}}}"
# shellcheck disable=SC1091
source "$EF2_ROOT/env.sh"
exec python "$EF2_ROOT/fold.py" --target "$HERE" --tag {tag} --gpu "${{1:-0}}"
"""


def parse_spec(spec):
    """'H:Y57F+L:Q27Y' -> {'H': [('Y', 57, 'F')], 'L': [('Q', 27, 'Y')]}"""
    by_chain = {}
    for token in re.split(r"[+,]", spec.strip()):
        token = token.strip()
        if not token:
            continue
        if ":" not in token:
            sys.exit(f"Bad mutation '{token}': expected CHAIN:WT<pos>MUT, e.g. H:Y57F")
        chain, mut = token.split(":", 1)
        chain = chain.strip().upper()
        if chain not in CHAINS:
            sys.exit(f"Unknown chain '{chain}' (allowed: H, L, A)")
        m = MUT_RE.match(mut.strip().upper())
        if not m:
            sys.exit(f"Bad mutation '{mut}': expected e.g. Y57F")
        by_chain.setdefault(chain, []).append((m.group(1), int(m.group(2)), m.group(3)))
    if not by_chain:
        sys.exit("No mutations parsed from the spec.")
    return by_chain


def apply_mutations(seq, muts, chain):
    out = list(seq)
    for wt, pos, new in muts:
        i = pos - 1
        if i < 0 or i >= len(out):
            sys.exit(f"{chain}:{wt}{pos}{new} out of range (chain {chain} is {len(out)} aa)")
        if out[i] != wt:
            sys.exit(f"{chain}:{wt}{pos}{new} WT mismatch: position {pos} is '{out[i]}', not '{wt}'")
        out[i] = new
    return "".join(out)


def folder_and_tag(by_chain):
    parts, tagparts = [], []
    for chain in ("H", "L", "A"):
        if chain in by_chain:
            muts = sorted(by_chain[chain], key=lambda m: m[1])
            parts.append(f"{chain}-" + "+".join(f"{w}{p}{n}" for w, p, n in muts))
            tagparts += [f"{chain}{w}{p}{n}" for w, p, n in muts]
    return "mut_" + "_".join(parts), "MUT_" + "_".join(tagparts)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("spec", help="e.g. 'H:Y57F+L:Q27Y'")
    ap.add_argument("--target", required=True, help="Parent target directory (the WT)")
    ap.add_argument("--out", default=None,
                    help="Where to create the mutant folder (default: <target>/mutants)")
    ap.add_argument("--print-only", action="store_true",
                    help="Print the mutated sequences; create nothing")
    ap.add_argument("--force", action="store_true", help="Overwrite an existing mutant folder")
    args = ap.parse_args()

    target = Path(args.target).resolve()
    parent = load_target_dir(target)
    by_chain = parse_spec(args.spec)

    seqs = {}
    for chain, key in CHAINS.items():
        if key not in parent:
            if chain in by_chain:
                sys.exit(f"Spec mutates chain {chain} but {target} has no {key} sequence.")
            continue
        seqs[key] = apply_mutations(parent[key], by_chain.get(chain, []), chain)

    folder, tag = folder_and_tag(by_chain)
    desc = " + ".join(f"{c}:{','.join(f'{w}{p}{n}' for w, p, n in by_chain[c])}"
                      for c in ("H", "L", "A") if c in by_chain)

    print(f"parent : {target}")
    print(f"spec   : {desc}")
    for chain, key in CHAINS.items():
        if key in seqs:
            state = "mutated" if chain in by_chain else "WT     "
            print(f"{key:<8} ({state}): {len(seqs[key])} aa"
                  + (f"  {seqs[key]}" if chain in by_chain else ""))
    if args.print_only:
        return

    mutdir = (Path(args.out).resolve() if args.out else target / "mutants") / folder
    if mutdir.exists() and not args.force:
        sys.exit(f"{mutdir} already exists (use --force to overwrite).")
    mutdir.mkdir(parents=True, exist_ok=True)

    for key, seq in seqs.items():
        (mutdir / f"{key}.txt").write_text(seq + "\n")

    # DockQ reference = the parent structure, so scores.csv answers "did the mutation
    # move the dock?". A relative symlink survives moving the whole tree.
    ref = target / "reference.cif"
    link = mutdir / "reference.cif"
    if ref.exists() and not link.exists():
        try:
            link.symlink_to(os.path.relpath(ref, mutdir))
        except OSError:
            link.write_bytes(ref.read_bytes())

    (mutdir / "run.sh").write_text(RUN_SH.format(desc=f"{target.name} mutant {desc}",
                                                 repo=HERE, tag=tag))
    (mutdir / "run.sh").chmod(0o755)

    rel = os.path.relpath(mutdir, Path.cwd())
    print(f"\ncreated {rel}")
    print(f"fold it: bash {rel}/run.sh 0        # last arg = GPU index")


if __name__ == "__main__":
    main()
