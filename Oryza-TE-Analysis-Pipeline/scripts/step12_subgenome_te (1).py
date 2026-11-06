#!/usr/bin/env python3
"""
step12_subgenome_te.py

TE content per SUBGENOME (the leaves of the step8 tree) and a test of TE expansion/contraction
between subgenomes of the same tetraploid.

Inputs
  te_density_windows.tsv  (step11)       TE fractions in 1 Mb windows
  subgenome_segments.tsv  (step9)        runs of subgenome labels along chromosomes
Method
  * each window is given the label of the step9 segment that contains its midpoint
    (windows outside any segment, e.g. on chromosomes with < 200 assigned genes, are ignored)
  * leaf value = mean TE fraction over the windows of that species + label (weighted by window length);
    officinalis_C = all windows of officinalis (diploid CC, no subgenome split)
  * test: for each tetraploid, units = (chromosome, label) with >= --min-win windows;
    chromosome-level means are compared between the two labels (BBCC: B vs C; CCDD: C vs D)
    by a label-permutation test across units (--perms, two-sided). Using chromosomes, not windows,
    as units avoids counting autocorrelated neighbouring windows as independent.
    Caveat: units of one tetraploid are not independent of the other species (shared ancestry);
    results are per species and not pooled.
Outputs in --outdir
  leaf_te_table.tsv       leaf, species, subgenome, n_windows, bp, pct_all, pct_LTR, pct_DNA, pct_intact, intact_share
  subgenome_te_tests.tsv  species, labels, category, n_units_1, n_units_2, mean_1, mean_2, diff, p_perm
  subgenome_te_report.txt readable summary

Usage:
    python3 step12_subgenome_te.py \
        --windows dual_students_te_oryza_nku/data/te_atlas_step11/te_density_windows.tsv \
        --segments dual_students_te_oryza_nku/data/te_atlas_step9/subgenome_segments.tsv \
        --outdir dual_students_te_oryza_nku/data/te_atlas_step12
Stdlib only.
"""

import argparse
import csv
import os
import random
from collections import defaultdict

CATS = ["all", "LTR", "DNA", "intact"]
PAIRS = {"minuta": ("B", "C"), "malampuzhaensis": ("B", "C"),
         "alta": ("C", "D"), "latifolia": ("C", "D"), "grandiglumis": ("C", "D")}


def wmean(rows, key):
    tw = sum(r["w"] for r in rows)
    return sum(r["w"] * r[key] for r in rows) / tw if tw else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", required=True)
    ap.add_argument("--segments", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--min-win", type=int, default=5)
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)

    segs = defaultdict(list)
    with open(args.segments) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            segs[(r["species"], r["seq"])].append((int(r["start"]), int(r["end"]), r["label"]))

    wins = defaultdict(list)       # (species, label) -> window rows
    units = defaultdict(list)      # (species, chrom, label) -> window rows
    with open(args.windows) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ws, we = int(r["win_start"]), int(r["win_end"])
            row = {"w": we - ws, "all": float(r["frac_all"]), "LTR": float(r["frac_LTR"]),
                   "DNA": float(r["frac_DNA"]), "intact": float(r["frac_intact"])}
            sp = r["species"]
            if sp == "officinalis":
                wins[(sp, "C")].append(row)
                continue
            mid = (ws + we) / 2
            lab = next((l for s, e, l in segs.get((sp, r["chrom"]), []) if s <= mid <= e), None)
            if lab is None:
                continue
            wins[(sp, lab)].append(row)
            units[(sp, r["chrom"], lab)].append(row)

    leaf_rows = []
    for (sp, lab), rows in sorted(wins.items()):
        bp = sum(r["w"] for r in rows)
        v = {c: 100 * wmean(rows, c) for c in CATS}
        share = v["intact"] / v["all"] if v["all"] else float("nan")
        leaf_rows.append([f"{sp}_{lab}", sp, lab, len(rows), bp] + [round(v[c], 2) for c in CATS] + [round(share, 3)])
    with open(os.path.join(args.outdir, "leaf_te_table.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["leaf", "species", "subgenome", "n_windows", "bp", "pct_all", "pct_LTR", "pct_DNA",
                    "pct_intact", "intact_share"])
        w.writerows(leaf_rows)

    test_rows, lines = [], ["TE % per subgenome leaf (all / LTR / DNA / intact):"]
    for r in leaf_rows:
        lines.append(f"  {r[0]:22s} windows={r[3]:4d}  {r[5]:6.2f} {r[6]:6.2f} {r[7]:6.2f} {r[8]:6.2f}")
    lines.append("\nSubgenome comparison (chromosome-level units, permutation test):")
    for sp, (l1, l2) in PAIRS.items():
        u1 = [rows for (s, ch, lab), rows in units.items() if s == sp and lab == l1 and len(rows) >= args.min_win]
        u2 = [rows for (s, ch, lab), rows in units.items() if s == sp and lab == l2 and len(rows) >= args.min_win]
        if len(u1) < 3 or len(u2) < 3:
            lines.append(f"  {sp}: too few units ({len(u1)} {l1}, {len(u2)} {l2})")
            continue
        for c in CATS:
            m1 = [wmean(u, c) for u in u1]
            m2 = [wmean(u, c) for u in u2]
            a, b = sum(m1) / len(m1), sum(m2) / len(m2)
            obs = a - b
            pool = m1 + m2
            n1, ge = len(m1), 0
            for _ in range(args.perms):
                rng.shuffle(pool)
                d = sum(pool[:n1]) / n1 - sum(pool[n1:]) / (len(pool) - n1)
                if abs(d) >= abs(obs) - 1e-12:
                    ge += 1
            p = (1 + ge) / (args.perms + 1)
            test_rows.append([sp, f"{l1}-{l2}", c, len(m1), len(m2), round(100 * a, 2), round(100 * b, 2),
                              round(100 * obs, 2), round(p, 4)])
            lines.append(f"  {sp:16s} {l1} vs {l2} {c:7s} n={len(m1)}/{len(m2)}  "
                         f"{100 * a:6.2f}% vs {100 * b:6.2f}%  diff={100 * obs:+6.2f} pp  p={p:.4f}")
    with open(os.path.join(args.outdir, "subgenome_te_tests.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "labels", "category", "n_units_1", "n_units_2", "mean_1_pct", "mean_2_pct",
                    "diff_pp", "p_perm"])
        w.writerows(test_rows)
    rep = "\n".join(lines)
    open(os.path.join(args.outdir, "subgenome_te_report.txt"), "w").write(rep + "\n")
    print(rep)
    print(f"\nFiles in {args.outdir}: leaf_te_table.tsv, subgenome_te_tests.tsv, subgenome_te_report.txt")


if __name__ == "__main__":
    main()
