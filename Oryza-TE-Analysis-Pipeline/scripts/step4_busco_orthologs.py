#!/usr/bin/env python3
"""
step4_busco_orthologs.py

Builds ortholog tables for macro-synteny (dot plots, inversion/breakpoint
detection) from the compleasm/BUSCO per-gene tables already sitting in cm/.
No genome alignment needed.

Reference  : O. officinalis (CC diploid) -> uses Status == Single only
Queries    : 5 allotetraploids (BBCC / CCDD) -> uses Single AND Duplicated rows
             (in a tetraploid most BUSCOs are expected to be Duplicated: one
             copy per subgenome, so keeping only Single would discard most of
             the data). Homeologs on the non-C subgenome will appear as a
             second, weaker diagonal -- that is real biology, not noise.

Confirmed cm/*.poales_odb12.tsv.gz columns (2026-10-14):
  Gene, Status, Sequence, Gene Start, Gene End, Strand, Score, Length,
  Identity, Fraction, Frameshift events, Best gene, Codons, protein_sequence

Also writes sequence-length summaries (N50, #sequences) from the .fai files,
because O. officinalis looks contig-level (BDMV01xxxxxx.1): if it is highly
fragmented, only small-scale rearrangements will be visible against it.

Run from /xdisk/dcopetti/mdalam/rte :
    python3 step4_busco_orthologs.py \
        --outdir dual_students_te_oryza_nku/data/te_atlas_step4

Stdlib only.
"""

import argparse
import csv
import gzip
import os
import sys
from collections import Counter, defaultdict

csv.field_size_limit(sys.maxsize)  # protein_sequence column is huge

REF = "GCA_008326285.1_Oryza_officinalis_CC"
QUERIES = [
    "GCA_048166525.1_Oryza_minuta_BBCC",
    "GCA_048564985.1_Oryza_malampuzhaensis_BBCC",
    "GCA_047899615.1_Oryza_alta_CCDD",
    "GCA_048174585.1_Oryza_latifolia_CCDD",
    "GCA_048188845.1_Oryza_grandiglumis_CCDD",
]


def short(genome):
    # GCA_..._Oryza_alta_CCDD -> alta
    parts = genome.split("_Oryza_")
    return parts[1].rsplit("_", 1)[0] if len(parts) == 2 else genome


def find_fai(root, genome):
    for cand in (f"gfai/{genome}.fa.gz.fai", f"gfa/{genome}.fa.fai"):
        p = os.path.join(root, cand)
        if os.path.isfile(p):
            return p
    return None


def read_fai(path):
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                out.append((f[0], int(f[1])))
    return out


def n50(lengths):
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2
    acc = 0
    for x in s:
        acc += x
        if acc >= half:
            return x
    return 0


def load_busco(path, keep_status):
    """Returns (dict busco_id -> list of loci, Counter of all statuses)."""
    loci = defaultdict(list)
    status_counts = Counter()
    with gzip.open(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            st = row.get("Status", "")
            status_counts[st] += 1
            if st not in keep_status:
                continue
            try:
                s, e = int(row["Gene Start"]), int(row["Gene End"])
            except (ValueError, KeyError):
                continue
            loci[row["Gene"]].append({
                "seq": row["Sequence"], "start": min(s, e), "end": max(s, e),
                "strand": row["Strand"], "status": st,
            })
    return loci, status_counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rte-dir", default=".")
    ap.add_argument("--outdir",
                    default="dual_students_te_oryza_nku/data/te_atlas_step4")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # ---- sequence length summaries -------------------------------------
    len_rows = []
    summary_rows = []
    for g in [REF] + QUERIES:
        fai = find_fai(args.rte_dir, g)
        if not fai:
            print(f"[WARN] no .fai found for {g}")
            continue
        seqs = read_fai(fai)
        lens = [l for _, l in seqs]
        total = sum(lens)
        # how many sequences are needed to cover 90% of the assembly
        acc, n90 = 0, 0
        for l in sorted(lens, reverse=True):
            acc += l
            n90 += 1
            if acc >= 0.9 * total:
                break
        summary_rows.append([g, len(seqs), total, n50(lens), n90])
        for name, l in seqs:
            len_rows.append([g, name, l])
        print(f"[SEQ] {short(g):15s} seqs={len(seqs):6d} total={total/1e6:8.1f} Mb "
              f"N50={n50(lens)/1e6:7.2f} Mb  seqs_for_90pct={n90}")

    with open(os.path.join(args.outdir, "genome_seq_lengths.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "seqid", "length"])
        w.writerows(len_rows)
    with open(os.path.join(args.outdir, "genome_assembly_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["genome", "n_seqs", "total_bp", "n50", "n_seqs_for_90pct"])
        w.writerows(summary_rows)

    # ---- BUSCO tables ----------------------------------------------------
    ref_path = os.path.join(args.rte_dir, "cm", f"{REF}.poales_odb12.tsv.gz")
    ref_loci, ref_status = load_busco(ref_path, keep_status={"Single"})
    print(f"\n[REF] {short(REF)}: statuses={dict(ref_status)}; "
          f"usable single-copy BUSCOs={len(ref_loci)}")

    for q in QUERIES:
        qpath = os.path.join(args.rte_dir, "cm", f"{q}.poales_odb12.tsv.gz")
        if not os.path.isfile(qpath):
            print(f"[WARN] missing {qpath}")
            continue
        q_loci, q_status = load_busco(qpath, keep_status={"Single", "Duplicated"})
        shared = sorted(set(ref_loci) & set(q_loci))
        print(f"[QRY] {short(q):15s} statuses={dict(q_status)}; "
              f"shared BUSCO ids with ref={len(shared)}")

        out = os.path.join(args.outdir, f"busco_orthologs_{short(q)}_vs_officinalis.tsv")
        n_rows = 0
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["busco_id", "ref_seq", "ref_start", "ref_end", "ref_strand",
                        "q_seq", "q_start", "q_end", "q_strand", "q_status",
                        "n_q_copies", "rel_orientation"])
            for bid in shared:
                r = ref_loci[bid][0]
                copies = q_loci[bid]
                for c in copies:
                    rel = "+" if r["strand"] == c["strand"] else "-"
                    w.writerow([bid, r["seq"], r["start"], r["end"], r["strand"],
                                c["seq"], c["start"], c["end"], c["strand"],
                                c["status"], len(copies), rel])
                    n_rows += 1
        print(f"      wrote {out} ({n_rows} rows)")

    print("\nDone. Paste the [SEQ], [REF] and [QRY] lines back.")
    print("Key question: is officinalis contiguous enough (N50, seqs_for_90pct)?")


if __name__ == "__main__":
    main()
