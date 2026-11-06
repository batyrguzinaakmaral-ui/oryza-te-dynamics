#!/usr/bin/env python3
"""
step3_full_genome_te_composition.py

Whole-genome TE composition (NOT intact-only) for the 6 CC-clade genomes,
parsed directly from each genome's full TEanno.gff3.

Why this script exists: edta/testats/ (used in step2) was built only from
*.intact.gff3 files -- structurally intact / recently active TE copies.
That undercounts total TE content badly (testats showed ~2% per genome vs.
an expected ~35-45% typical for rice-relative genomes). This script parses
the FULL *.fa.mod.EDTA.TEanno.gff3 (homology + structural calls combined)
to get genome-wide % TE per classification.

Confirmed GFF3 format (from `grep -v '^#' <file> | head -5` on all 6 files,
2026-10-08):
    col1 seqid, col2 source(EDTA), col3 feature_type, col4 start, col5 end,
    col6 score, col7 strand, col8 phase(.), col9 attributes
    attributes include: ID=...;Name=...;classification=<Class/Subclass>;
    sequence_ontology=...;identity=...;method=homology|structural

CAVEAT (important, read before trusting absolute %): EDTA's full TEanno
combines structural-search and homology-search calls, and these CAN
overlap on the same locus (e.g. a structural LTR call overlapping a
homology call for the same copy). Counting both naively double-counts bp.
To reduce this, this script merges overlapping intervals PER
(seqid, classification) before summing bp -- i.e. it collapses overlaps
within the same TE class on the same contig, so a region is not counted
twice for the same classification. It does NOT attempt to resolve
overlaps BETWEEN different classifications (e.g. a DTA call overlapping a
Helitron call) -- that would need the original EDTA conflict-resolution
logic, which this script does not reimplement. So totals are a good-faith
non-redundant-within-class estimate, not a certified masking-free total.
Flag this caveat in the thesis methods section.

Usage (run from /xdisk/dcopetti/mdalam/rte):
    python3 step3_full_genome_te_composition.py \
        --edta-dir edta \
        --genome-annot edta/testats/genome_annot.tsv \
        --outdir dual_students_te_oryza_nku/data/te_atlas_step3

Stdlib only -- no pandas, no pip needed.
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

CC_CLADE_GENOMES = {
    "GCA_008326285.1_Oryza_officinalis_CC":        ("CC",   "Oryza_officinalis"),
    "GCA_048166525.1_Oryza_minuta_BBCC":           ("BBCC", "Oryza_minuta"),
    "GCA_048564985.1_Oryza_malampuzhaensis_BBCC":  ("BBCC", "Oryza_malampuzhaensis"),
    "GCA_047899615.1_Oryza_alta_CCDD":             ("CCDD", "Oryza_alta"),
    "GCA_048174585.1_Oryza_latifolia_CCDD":        ("CCDD", "Oryza_latifolia"),
    "GCA_048188845.1_Oryza_grandiglumis_CCDD":     ("CCDD", "Oryza_grandiglumis"),
}

# genome_annot.tsv (testats/) has NO header: genome, species, genome_type,
# mating_system, n_intact_elements, total_bp_intact, genome_size
GENOME_ANNOT_COLS = ["genome", "species", "genome_type", "mating_system",
                     "n_intact_elements", "total_bp_intact", "genome_size"]


def load_genome_sizes(path):
    sizes = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or len(row) < len(GENOME_ANNOT_COLS):
                continue
            rec = dict(zip(GENOME_ANNOT_COLS, row))
            try:
                sizes[rec["genome"]] = int(float(rec["genome_size"]))
            except ValueError:
                continue
    return sizes


def parse_attributes(attr_field):
    out = {}
    for part in attr_field.strip().split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k] = v
    return out


def merge_intervals(intervals):
    """intervals: list of (start, end) 1-based inclusive. Returns merged bp sum."""
    if not intervals:
        return 0
    intervals.sort()
    merged_bp = 0
    cur_start, cur_end = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_end + 1:  # overlap or adjacent
            cur_end = max(cur_end, e)
        else:
            merged_bp += cur_end - cur_start + 1
            cur_start, cur_end = s, e
    merged_bp += cur_end - cur_start + 1
    return merged_bp


def parse_gff3(path):
    """
    Returns:
      n_elements_raw: dict classification -> raw line count (no merging)
      feature_type_of: dict classification -> feature_type (as seen in col3)
      intervals: dict (seqid, classification) -> list of (start, end)
    """
    n_elements_raw = defaultdict(int)
    feature_type_of = {}
    intervals = defaultdict(list)

    with open(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            seqid, source, ftype, start, end = cols[0], cols[1], cols[2], cols[3], cols[4]
            attrs = parse_attributes(cols[8])
            classification = attrs.get("classification")
            if not classification:
                continue
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            if e < s:
                s, e = e, s
            n_elements_raw[classification] += 1
            feature_type_of.setdefault(classification, ftype)
            intervals[(seqid, classification)].append((s, e))

    return n_elements_raw, feature_type_of, intervals


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--edta-dir", default="edta",
                     help="Path to the edta/ directory (default: edta)")
    ap.add_argument("--genome-annot", default="edta/testats/genome_annot.tsv",
                     help="Path to genome_annot.tsv for genome_size lookup")
    ap.add_argument("--outdir", default="dual_students_te_oryza_nku/data/te_atlas_step3",
                     help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"[INFO] Loading genome sizes from {args.genome_annot} ...")
    genome_sizes = load_genome_sizes(args.genome_annot)
    missing_sizes = [g for g in CC_CLADE_GENOMES if g not in genome_sizes]
    if missing_sizes:
        print(f"[WARN] No genome_size found for: {missing_sizes}")

    long_rows = []  # genome, species, genome_type, classification, feature_type,
                     # n_elements_raw, total_bp_merged, genome_size, pct_of_genome

    for genome_id, (genome_type, species) in CC_CLADE_GENOMES.items():
        gff3_path = os.path.join(args.edta_dir, genome_id,
                                  f"{genome_id}.fa.mod.EDTA.TEanno.gff3")
        print(f"\n[INFO] Parsing {genome_id} ...")
        if not os.path.isfile(gff3_path):
            print(f"[ERROR] Not found: {gff3_path} -- skipping")
            continue

        n_elements_raw, feature_type_of, intervals = parse_gff3(gff3_path)

        genome_size = genome_sizes.get(genome_id)
        if not genome_size:
            print(f"[ERROR] No genome_size for {genome_id} -- skipping pct calc")
            continue

        # Merge overlaps within (seqid, classification), then sum across seqids
        # per classification.
        bp_by_class = defaultdict(int)
        for (seqid, classification), ivs in intervals.items():
            bp_by_class[classification] += merge_intervals(ivs)

        total_classes = len(bp_by_class)
        total_bp_all = sum(bp_by_class.values())
        print(f"  classes found: {total_classes}, "
              f"raw elements: {sum(n_elements_raw.values())}, "
              f"merged TE bp: {total_bp_all:,} / genome_size {genome_size:,} "
              f"= {100*total_bp_all/genome_size:.2f}%")

        for classification, bp in sorted(bp_by_class.items()):
            long_rows.append({
                "genome": genome_id,
                "species": species,
                "genome_type": genome_type,
                "classification": classification,
                "feature_type": feature_type_of.get(classification, ""),
                "n_elements_raw": n_elements_raw.get(classification, 0),
                "total_bp_merged": bp,
                "genome_size": genome_size,
                "pct_of_genome": round(100 * bp / genome_size, 6),
            })

    if not long_rows:
        print("[ERROR] No data parsed -- check paths above.")
        sys.exit(1)

    # ---- Write long-format table ----
    long_path = os.path.join(args.outdir, "full_genome_per_class.tsv")
    fieldnames = ["genome", "species", "genome_type", "classification",
                  "feature_type", "n_elements_raw", "total_bp_merged",
                  "genome_size", "pct_of_genome"]
    with open(long_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(long_rows)
    print(f"\n[INFO] Wrote {long_path} ({len(long_rows)} rows)")

    # ---- Write wide pct table (classification x species) ----
    classes = sorted(set(r["classification"] for r in long_rows))
    species_order = [CC_CLADE_GENOMES[g][1] for g in CC_CLADE_GENOMES]
    pct_lookup = {(r["classification"], r["species"]): r["pct_of_genome"]
                  for r in long_rows}
    ftype_lookup = {r["classification"]: r["feature_type"] for r in long_rows}

    wide_path = os.path.join(args.outdir, "full_genome_class_wide_pct.tsv")
    with open(wide_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["classification", "feature_type"] + species_order)
        for c in classes:
            row = [c, ftype_lookup.get(c, "")]
            for sp in species_order:
                row.append(pct_lookup.get((c, sp), 0.0))
            writer.writerow(row)
    print(f"[INFO] Wrote {wide_path}")

    # ---- Report ----
    report_path = os.path.join(args.outdir, "full_genome_report.txt")
    with open(report_path, "w") as fh:
        fh.write("Full-genome (non-intact-only) TE composition -- CC-clade\n")
        fh.write("=" * 60 + "\n\n")
        fh.write("CAVEAT: overlaps merged WITHIN each classification only, not\n")
        fh.write("across classifications. See script docstring.\n\n")
        for genome_id, (genome_type, species) in CC_CLADE_GENOMES.items():
            rows = [r for r in long_rows if r["genome"] == genome_id]
            if not rows:
                fh.write(f"{genome_id}: NOT PARSED (see errors above)\n\n")
                continue
            total_pct = sum(r["pct_of_genome"] for r in rows)
            fh.write(f"{genome_id} ({species}, {genome_type})\n")
            fh.write(f"  Total TE: {total_pct:.2f}% of genome "
                      f"(sum across {len(rows)} classes, within-class "
                      f"overlap-merged)\n")
            for r in sorted(rows, key=lambda x: -x["pct_of_genome"])[:5]:
                fh.write(f"    {r['classification']:15s} {r['pct_of_genome']:6.2f}%  "
                          f"({r['feature_type']})\n")
            fh.write("\n")
    print(f"[INFO] Wrote {report_path}")
    print("\nDone. Compare full_genome_class_wide_pct.tsv (this script) against")
    print("cc_clade_class_wide_pct.tsv (step2, intact-only) -- the gap between")
    print("them IS the 'decayed / non-intact' TE fraction per class.")


if __name__ == "__main__":
    main()
