#!/usr/bin/env python3
"""
step5_synteny_blocks.py

From the BUSCO ortholog tables (step4) call synteny blocks and rearrangement
breakpoints between O. officinalis (CC reference) and each allotetraploid.
No genome alignment: the only evidence is the order and strand of shared
single-copy genes.

Method
  1. Keep chromosome-scale sequences only (default >= 10 Mb) in ref and query.
  2. For every query chromosome find its dominant reference partner
     (the reference chromosome holding most of its orthologs) and the
     "purity" (share of its orthologs on that partner). A tetraploid should
     show two query chromosomes per reference chromosome (the two
     homeologous subgenomes); purity well below ~0.9 flags a translocation.
  3. Along each query chromosome label every gene (ref_chr, relative strand),
     smooth with a majority vote over --window neighbours, and cut the chromosome
     into runs of identical label.
  4. A boundary between two runs is a breakpoint:
        inversion     same ref chromosome, strand flips
        translocation ref chromosome changes
     Runs shorter than --min-genes are absorbed by their neighbours.
  5. breakpoints.tsv also carries left_ref_pos / right_ref_pos: the O. officinalis
     coordinate (gene midpoint) of the two genes flanking each breakpoint. step7
     uses them to recognise the SAME rearrangement seen in several species or on
     both homeologous chromosomes of one species.

RESOLUTION CAVEAT: BUSCO genes are ~1 per 90-160 kb, so a breakpoint is an
interval between two flanking genes (reported as bp_width), not a base-pair
position, and inversions spanning only a few genes are invisible. This finds
large rearrangements only. Say so in the thesis methods.

Run from /xdisk/dcopetti/mdalam/rte:
    python3 step5_synteny_blocks.py \
        --step4-dir dual_students_te_oryza_nku/data/te_atlas_step4 \
        --outdir    dual_students_te_oryza_nku/data/te_atlas_step5
Stdlib only.
"""

import argparse
import csv
import glob
import os
import re
from collections import Counter, defaultdict


def read_lengths(path):
    lens = defaultdict(dict)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lens[r["genome"]][r["seqid"]] = int(r["length"])
    return lens


def genome_of(lens, short_name):
    for g in lens:
        if f"_Oryza_{short_name}_" in g:
            return g
    return None


def majority_smooth(labels, window):
    half = window // 2
    out = []
    for i, lab in enumerate(labels):
        win = labels[max(0, i - half): i + half + 1]
        c = Counter(win).most_common()
        if len(c) > 1 and c[0][1] == c[1][1]:
            out.append(lab)           # tie: keep original
        else:
            out.append(c[0][0])
    return out


def runs_of(labels):
    runs, start = [], 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            runs.append((start, i - 1, labels[start]))
            start = i
    return runs


def absorb_short(runs, genes_labels, min_genes):
    """Merge runs shorter than min_genes into the larger neighbour."""
    changed = True
    while changed and len(runs) > 1:
        changed = False
        for k, (a, b, lab) in enumerate(runs):
            if b - a + 1 < min_genes:
                left = runs[k - 1] if k > 0 else None
                right = runs[k + 1] if k < len(runs) - 1 else None
                tgt = None
                if left and right:
                    tgt = left if (left[1] - left[0]) >= (right[1] - right[0]) else right
                else:
                    tgt = left or right
                for i in range(a, b + 1):
                    genes_labels[i] = tgt[2]
                runs = runs_of(genes_labels)
                changed = True
                break
    return runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--step4-dir", default="dual_students_te_oryza_nku/data/te_atlas_step4")
    ap.add_argument("--outdir", default="dual_students_te_oryza_nku/data/te_atlas_step5")
    ap.add_argument("--min-chr-mb", type=float, default=10.0)
    ap.add_argument("--window", type=int, default=5)
    ap.add_argument("--min-genes", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    lens = read_lengths(os.path.join(args.step4_dir, "genome_seq_lengths.tsv"))
    min_len = int(args.min_chr_mb * 1e6)

    partner_rows, block_rows, bp_rows = [], [], []

    for path in sorted(glob.glob(os.path.join(args.step4_dir, "busco_orthologs_*_vs_officinalis.tsv"))):
        qshort = re.search(r"busco_orthologs_(.+)_vs_officinalis", os.path.basename(path)).group(1)
        qgen = genome_of(lens, qshort)
        rgen = genome_of(lens, "officinalis")
        qchr = {s for s, l in lens[qgen].items() if l >= min_len}
        rchr = {s for s, l in lens[rgen].items() if l >= min_len}
        print(f"\n[{qshort}] chromosome-scale seqs: query={len(qchr)} ref={len(rchr)} "
              f"(expect ~24 and ~12)")

        by_q = defaultdict(list)
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if r["q_seq"] in qchr and r["ref_seq"] in rchr:
                    by_q[r["q_seq"]].append(r)

        for qs in sorted(by_q, key=lambda s: -lens[qgen][s]):
            rows = sorted(by_q[qs], key=lambda r: int(r["q_start"]))
            partner, n_partner = Counter(r["ref_seq"] for r in rows).most_common(1)[0]
            purity = n_partner / len(rows)
            partner_rows.append([qshort, qs, lens[qgen][qs], len(rows), partner,
                                 round(purity, 3)])

            labels = [(r["ref_seq"], r["rel_orientation"]) for r in rows]
            sm = majority_smooth(labels, args.window)
            runs = absorb_short(runs_of(sm), sm, args.min_genes)

            for (a, b, lab) in runs:
                seg = rows[a:b + 1]
                block_rows.append([
                    qshort, qs, min(int(r["q_start"]) for r in seg),
                    max(int(r["q_end"]) for r in seg), lab[0],
                    min(int(r["ref_start"]) for r in seg),
                    max(int(r["ref_end"]) for r in seg), lab[1], len(seg)])
            for (a1, b1, l1), (a2, b2, l2) in zip(runs, runs[1:]):
                left, right = rows[b1], rows[a2]
                kind = "translocation" if l1[0] != l2[0] else "inversion"
                bs, be = int(left["q_end"]), int(right["q_start"])
                lpos = (int(left["ref_start"]) + int(left["ref_end"])) // 2
                rpos = (int(right["ref_start"]) + int(right["ref_end"])) // 2
                bp_rows.append([qshort, qs, kind, bs, be, max(be - bs, 0),
                                l1[0], l1[1], l2[0], l2[1], b1 - a1 + 1, b2 - a2 + 1,
                                lpos, rpos])

    def w(name, header, rows):
        with open(os.path.join(args.outdir, name), "w", newline="") as fh:
            wr = csv.writer(fh, delimiter="\t")
            wr.writerow(header)
            wr.writerows(rows)

    w("chromosome_partners.tsv",
      ["query", "q_chr", "q_len", "n_orthologs", "ref_partner", "purity"], partner_rows)
    w("synteny_blocks.tsv",
      ["query", "q_chr", "q_start", "q_end", "ref_chr", "ref_start", "ref_end",
       "rel_orientation", "n_genes"], block_rows)
    w("breakpoints.tsv",
      ["query", "q_chr", "type", "bp_start", "bp_end", "bp_width", "left_ref_chr",
       "left_orient", "right_ref_chr", "right_orient", "n_genes_left", "n_genes_right",
       "left_ref_pos", "right_ref_pos"],
      bp_rows)

    print("\n=== summary ===")
    for q in sorted({r[0] for r in partner_rows}):
        pr = [r for r in partner_rows if r[0] == q]
        low = [r for r in pr if r[5] < 0.9]
        n_inv = sum(1 for r in bp_rows if r[0] == q and r[2] == "inversion")
        n_tr = sum(1 for r in bp_rows if r[0] == q and r[2] == "translocation")
        per_ref = Counter(r[4] for r in pr)
        print(f"{q:16s} q_chrs={len(pr):2d}  ref chrs hit={len(per_ref):2d}  "
              f"chrs with purity<0.9: {len(low)}  "
              f"inversion breakpoints={n_inv}  translocation breakpoints={n_tr}")
    print("\nPositive control: bioRxiv/Nat Genet report a reciprocal translocation "
          "between Chr1BB and Chr1CC in BOTH O. minuta and O. malampuzhaensis. "
          "If the method works, those two species should each show low-purity "
          "chromosomes / translocation breakpoints in the same place.")


if __name__ == "__main__":
    main()
