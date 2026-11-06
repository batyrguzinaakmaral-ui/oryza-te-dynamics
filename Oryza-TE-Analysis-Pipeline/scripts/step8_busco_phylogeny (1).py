#!/usr/bin/env python3
"""
step8_busco_phylogeny.py

Species tree of the CC-clade from the BUSCO (poales_odb12) protein sequences that are
already stored in cm/*.poales_odb12.tsv.gz (column protein_sequence). Stdlib only:
no aligner / tree program is needed on the HPC.

Method (alignment-free, appropriate for closely related Oryza proteins):
  1. Orthologs = BUSCO ids that are single-copy in both diploid anchors
       officinalis (CC)  GCA_008326285.1
       punctata    (BB)  GCA_000573905.2
  2. Homeolog-aware leaves for tetraploids (two copies per BUSCO = two subgenomes):
       BBCC (minuta, malampuzhaensis): the copy closer to officinalis is called C, the other B
                                       (best of the two assignments, by summed k-mer distance)
       CCDD (alta, latifolia, grandiglumis): the copy closer to officinalis is C, the other D
     Genes where the two assignments are nearly tied (--margin) are skipped for that species.
     Leaves: officinalis_C, punctata_B, minuta_B, minuta_C, ..., alta_C, alta_D, ...
     (Assigning C by similarity to officinalis is partly circular for the C side; the
      non-circular check is whether minuta_B and malampuzhaensis_B end up closer to each
      other and to punctata_B than to any C leaf.)
  3. Distance between two protein copies: Mash-like distance on amino-acid k-mers
       d = -ln(2J/(1+J))/k,  J = Jaccard index of the k-mer sets  (k = --k, default 5)
     Leaf-to-leaf distance = mean over genes where both leaves exist.
  4. Neighbor-joining tree (unrooted) written as Newick.

Outputs in --outdir:
  tree.nwk                      Newick (view in R: ape::read.tree + plot, or iTOL / FigTree)
  leaf_distance_matrix.tsv      mean distances
  leaf_pair_gene_counts.tsv     how many genes support each distance
  homeolog_assignment.tsv       busco_id, species, subgenome, seq, start, end, strand
                                (genomic positions of every assigned copy: input for a
                                 subgenome-aware synteny / positive-control analysis)
  report.txt                    counts and sanity checks
  fasta/<busco_id>.faa          (only with --write-fasta) per-gene FASTA with leaf names,
                                for MAFFT + IQ-TREE if those exist on the server

Usage (from /xdisk/dcopetti/mdalam/rte):
    python3 dual_students_te_oryza_nku/data/step8_busco_phylogeny.py \
        --cm-dir cm --outdir dual_students_te_oryza_nku/data/te_atlas_step8
"""

import argparse
import csv
import glob
import gzip
import itertools
import math
import os
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)

TAXA = {
    "officinalis": ("GCA_008326285.1", "C"),
    "punctata": ("GCA_000573905.2", "B"),
    "minuta": ("GCA_048166525.1", "BBCC"),
    "malampuzhaensis": ("GCA_048564985.1", "BBCC"),
    "alta": ("GCA_047899615.1", "CCDD"),
    "latifolia": ("GCA_048174585.1", "CCDD"),
    "grandiglumis": ("GCA_048188845.1", "CCDD"),
}


def find_cm(cm_dir, acc):
    hits = sorted(glob.glob(os.path.join(cm_dir, f"{acc}*.poales_odb12.tsv.gz")))
    return hits[0] if hits else None


def load(path):
    out = defaultdict(list)
    with gzip.open(path, "rt", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r.get("Status") not in ("Single", "Duplicated"):
                continue
            seq = (r.get("protein_sequence") or "").replace("*", "").upper()
            if len(seq) < 50:
                continue
            try:
                s, e = int(r["Gene Start"]), int(r["Gene End"])
            except (ValueError, KeyError):
                continue
            out[r["Gene"]].append({"seq": seq, "chrom": r["Sequence"], "start": min(s, e),
                                   "end": max(s, e), "strand": r["Strand"], "status": r["Status"]})
    return out


def neighbor_joining(names, D):
    dist = dict(D)

    def g(a, b):
        return dist[(a, b)] if (a, b) in dist else dist[(b, a)]

    active = list(names)
    tree = {n: n for n in names}
    cnt = 0
    while len(active) > 3:
        n = len(active)
        r = {a: sum(g(a, b) for b in active if b != a) for a in active}
        best = None
        for a, b in itertools.combinations(active, 2):
            q = (n - 2) * g(a, b) - r[a] - r[b]
            if best is None or q < best[0]:
                best = (q, a, b)
        _, a, b = best
        dab = g(a, b)
        la = max(0.0, 0.5 * dab + (r[a] - r[b]) / (2 * (n - 2)))
        lb = max(0.0, dab - la)
        cnt += 1
        u = f"_n{cnt}"
        tree[u] = f"({tree[a]}:{la:.5f},{tree[b]}:{lb:.5f})"
        for c in active:
            if c not in (a, b):
                dist[(u, c)] = 0.5 * (g(a, c) + g(b, c) - dab)
        active = [c for c in active if c not in (a, b)] + [u]
    a, b, c = active
    dab, dac, dbc = g(a, b), g(a, c), g(b, c)
    la = max(0.0, 0.5 * (dab + dac - dbc))
    lb = max(0.0, 0.5 * (dab + dbc - dac))
    lc = max(0.0, 0.5 * (dac + dbc - dab))
    return f"({tree[a]}:{la:.5f},{tree[b]}:{lb:.5f},{tree[c]}:{lc:.5f});"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cm-dir", default="cm")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--margin", type=float, default=0.002,
                    help="skip a tetraploid gene when the two homeolog assignments differ by less than this")
    ap.add_argument("--min-genes", type=int, default=30,
                    help="drop a leaf pair supported by fewer genes (the leaf is dropped if it has no pair left)")
    ap.add_argument("--write-fasta", action="store_true")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    K = args.k

    def kmers(s):
        return {s[i:i + K] for i in range(len(s) - K + 1)}

    def dist(a, b):
        inter = len(a & b)
        union = len(a) + len(b) - inter
        if not union:
            return 3.0
        j = inter / union
        return 3.0 if j == 0 else -math.log(2 * j / (1 + j)) / K

    data = {}
    for sp, (acc, kind) in TAXA.items():
        p = find_cm(args.cm_dir, acc)
        if not p:
            print(f"[WARN] no cm table for {sp} ({acc})")
            continue
        data[sp] = load(p)
        print(f"[LOAD] {sp:16s} {acc}  BUSCOs with protein: {len(data[sp])}")
    for need in ("officinalis", "punctata"):
        if need not in data:
            sys.exit(f"[ERROR] anchor {need} missing; both officinalis (CC) and punctata (BB) are needed")

    off, pun = data["officinalis"], data["punctata"]
    genes = [g for g in off if len(off[g]) == 1 and off[g][0]["status"] == "Single"
             and g in pun and len(pun[g]) == 1 and pun[g][0]["status"] == "Single"]
    print(f"[INFO] BUSCOs single-copy in both anchors: {len(genes)}")

    sums = defaultdict(float)
    cnts = defaultdict(int)
    assign_rows = []
    stats = defaultdict(lambda: {"used": 0, "tied": 0, "not2": 0})
    if args.write_fasta:
        os.makedirs(os.path.join(args.outdir, "fasta"), exist_ok=True)

    for g in genes:
        leaves = {"officinalis_C": off[g][0], "punctata_B": pun[g][0]}
        ks = {name: kmers(c["seq"]) for name, c in leaves.items()}
        for sp, (acc, kind) in TAXA.items():
            if kind not in ("BBCC", "CCDD") or sp not in data:
                continue
            copies = data[sp].get(g, [])
            if len(copies) != 2:
                stats[sp]["not2"] += 1
                continue
            k1, k2 = kmers(copies[0]["seq"]), kmers(copies[1]["seq"])
            if kind == "BBCC":
                o1 = dist(k1, ks["officinalis_C"]) + dist(k2, ks["punctata_B"])   # copy1 = C, copy2 = B
                o2 = dist(k2, ks["officinalis_C"]) + dist(k1, ks["punctata_B"])   # copy2 = C, copy1 = B
                if abs(o1 - o2) < args.margin:
                    stats[sp]["tied"] += 1
                    continue
                c_idx = 0 if o1 < o2 else 1
                other = "B"
            else:
                d1, d2 = dist(k1, ks["officinalis_C"]), dist(k2, ks["officinalis_C"])
                if abs(d1 - d2) < args.margin:
                    stats[sp]["tied"] += 1
                    continue
                c_idx = 0 if d1 < d2 else 1
                other = "D"
            stats[sp]["used"] += 1
            for idx, sub in ((c_idx, "C"), (1 - c_idx, other)):
                leaves[f"{sp}_{sub}"] = copies[idx]
                ks[f"{sp}_{sub}"] = k1 if idx == 0 else k2
                c = copies[idx]
                assign_rows.append([g, sp, sub, c["chrom"], c["start"], c["end"], c["strand"]])
        for sp in ("officinalis", "punctata"):
            c = leaves["officinalis_C" if sp == "officinalis" else "punctata_B"]
            assign_rows.append([g, sp, "C" if sp == "officinalis" else "B",
                                c["chrom"], c["start"], c["end"], c["strand"]])
        names = sorted(leaves)
        for a, b in itertools.combinations(names, 2):
            sums[(a, b)] += dist(ks[a], ks[b])
            cnts[(a, b)] += 1
        if args.write_fasta:
            with open(os.path.join(args.outdir, "fasta", f"{g}.faa"), "w") as fh:
                for n in names:
                    fh.write(f">{n}\n{leaves[n]['seq']}\n")

    all_leaves = sorted({x for pair in cnts for x in pair})
    D, N = {}, {}
    for (a, b), c in cnts.items():
        if c >= args.min_genes:
            D[(a, b)] = sums[(a, b)] / c
            N[(a, b)] = c
    leaves_ok = [x for x in all_leaves
                 if all(((x, y) in D or (y, x) in D) for y in all_leaves if y != x)]
    dropped = [x for x in all_leaves if x not in leaves_ok]

    def get(a, b):
        return D[(a, b)] if (a, b) in D else D[(b, a)]

    with open(os.path.join(args.outdir, "leaf_distance_matrix.tsv"), "w") as fh:
        fh.write("leaf\t" + "\t".join(leaves_ok) + "\n")
        for a in leaves_ok:
            fh.write(a + "\t" + "\t".join("0" if a == b else f"{get(a, b):.5f}" for b in leaves_ok) + "\n")
    with open(os.path.join(args.outdir, "leaf_pair_gene_counts.tsv"), "w") as fh:
        fh.write("leaf_a\tleaf_b\tn_genes\n")
        for (a, b), c in sorted(N.items()):
            fh.write(f"{a}\t{b}\t{c}\n")
    with open(os.path.join(args.outdir, "homeolog_assignment.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["busco_id", "species", "subgenome", "seq", "start", "end", "strand"])
        w.writerows(assign_rows)

    lines = [f"BUSCOs single-copy in both anchors: {len(genes)}",
             "tetraploid genes used / skipped (tied or not exactly 2 copies):"]
    for sp, s in stats.items():
        lines.append(f"  {sp:16s} used={s['used']:5d} tied={s['tied']:5d} not_2_copies={s['not2']:5d}")
    if dropped:
        lines.append(f"leaves dropped for lack of data: {dropped}")
    nwk = ""
    if len(leaves_ok) >= 4:
        nwk = neighbor_joining(leaves_ok, {(a, b): get(a, b) for a, b in itertools.combinations(leaves_ok, 2)})
        open(os.path.join(args.outdir, "tree.nwk"), "w").write(nwk + "\n")
        lines.append("\nNewick:\n" + nwk)
    else:
        lines.append("fewer than 4 usable leaves: no tree")
    # sanity check: B leaves of the two BBCC species should be closer to each other / to punctata_B than to C leaves
    chk = [("minuta_B", "malampuzhaensis_B"), ("minuta_B", "punctata_B"), ("minuta_B", "minuta_C"),
           ("minuta_B", "officinalis_C"), ("minuta_C", "officinalis_C"), ("minuta_C", "malampuzhaensis_C")]
    lines.append("\nsanity distances (mean aa-kmer distance):")
    for a, b in chk:
        if a in leaves_ok and b in leaves_ok:
            lines.append(f"  {a:20s} {b:22s} {get(a, b):.4f}")
    rep = "\n".join(lines)
    open(os.path.join(args.outdir, "report.txt"), "w").write(rep + "\n")
    print(rep)
    print(f"\nFiles in {args.outdir}")


if __name__ == "__main__":
    main()
