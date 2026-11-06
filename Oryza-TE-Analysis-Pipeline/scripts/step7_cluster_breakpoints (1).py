#!/usr/bin/env python3
"""
step7_cluster_breakpoints.py

Collapses breakpoint records into rearrangement EVENTS, because one rearrangement
is reported many times:
  * in every species that carries it (e.g. a junction that is really specific to the
    O. officinalis reference appears in all five tetraploids),
  * on both homeologous chromosomes of one tetraploid,
  * at both ends of one inversion (two breakpoints, same event).
Treating those records as independent inflates n and makes p-values far too small.

Two records belong to the same event when they have the same type and the same
reference chromosome(s), and the O. officinalis positions of their flanking genes agree
within --tol-mb (default 1 Mb):
  inversion      (ref_chr, lo, hi)   lo/hi = the two flanking-gene positions, sorted
  translocation  ((chrA,posA),(chrB,posB)) sorted, so reading direction does not matter
Single-linkage clustering (union-find).

Event class (relative to the number of species in the file):
  shared_all     seen in every species     -> most likely a property of the reference
                                              (officinalis lineage or assembly), NOT a
                                              rearrangement of the tetraploids
  shared_subset  seen in 2..n-1 species    -> phylogenetically informative
  private        seen in one species only

Output weight = 1 / (number of records in the event), so each event counts once in step6.
This is conservative: homologous loci in different species could be called semi-independent
observations of TE environment, but the rearrangement itself happened once.

Needs breakpoints.tsv from the PATCHED step5 (columns left_ref_pos, right_ref_pos).

Usage:
    python3 step7_cluster_breakpoints.py \
        --breakpoints dual_students_te_oryza_nku/data/te_atlas_step5_mg20/breakpoints.tsv \
        --outdir dual_students_te_oryza_nku/data/te_atlas_step7_mg20
Stdlib only.
"""

import argparse
import csv
import os
import statistics
import sys
from collections import Counter, defaultdict

GROUPS = {"minuta": "BBCC", "malampuzhaensis": "BBCC",
          "alta": "CCDD", "grandiglumis": "CCDD", "latifolia": "CCDD"}


def signature(r):
    lp, rp = int(r["left_ref_pos"]), int(r["right_ref_pos"])
    if r["type"] == "inversion":
        lo, hi = sorted((lp, rp))
        return ("inversion", r["left_ref_chr"], lo, r["left_ref_chr"], hi)
    a, b = sorted([(r["left_ref_chr"], lp), (r["right_ref_chr"], rp)])
    return ("translocation", a[0], a[1], b[0], b[1])


def similar(s, t, tol):
    return (s[0] == t[0] and s[1] == t[1] and s[3] == t[3]
            and abs(s[2] - t[2]) <= tol and abs(s[4] - t[4]) <= tol)


def short_chr(c):
    return c.replace("BDMV010000", "R").replace(".1", "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--breakpoints", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tol-mb", type=float, default=1.0)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    tol = int(args.tol_mb * 1e6)

    with open(args.breakpoints) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rows, fields = list(rd), rd.fieldnames
    if not {"left_ref_pos", "right_ref_pos"} <= set(fields):
        sys.exit("breakpoints.tsv has no left_ref_pos/right_ref_pos: re-run the PATCHED step5 first")

    sigs = [signature(r) for r in rows]
    parent = list(range(len(rows)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(len(rows)):
        for j in range(i + 1, len(rows)):
            if similar(sigs[i], sigs[j], tol):
                parent[find(i)] = find(j)

    clusters = defaultdict(list)
    for i in range(len(rows)):
        clusters[find(i)].append(i)

    species_all = sorted({r["query"] for r in rows})
    n_sp = len(species_all)

    ev_rows, ev_of = [], {}
    order = sorted(clusters.values(), key=lambda idx: (sigs[idx[0]][0], sigs[idx[0]][1], sigs[idx[0]][2]))
    for k, idx in enumerate(order, 1):
        eid = f"E{k:03d}"
        sp = sorted({rows[i]["query"] for i in idx})
        cls = "shared_all" if len(sp) == n_sp else ("shared_subset" if len(sp) > 1 else "private")
        typ = sigs[idx[0]][0]
        s = [sigs[i] for i in idx]
        pos1 = int(statistics.median(x[2] for x in s))
        pos2 = int(statistics.median(x[4] for x in s))
        groups = sorted({GROUPS.get(x, "?") for x in sp})
        ev_rows.append([eid, typ, short_chr(s[0][1]), pos1, short_chr(s[0][3]), pos2,
                        len(idx), len(sp), "+".join(sp), "/".join(groups), cls])
        for i in idx:
            ev_of[i] = (eid, cls, len(sp), 1.0 / len(idx))

    with open(os.path.join(args.outdir, "breakpoint_events.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["event_id", "type", "ref_chr_1", "ref_pos_1", "ref_chr_2", "ref_pos_2",
                    "n_records", "n_species", "species", "genome_groups", "event_class"])
        w.writerows(ev_rows)

    with open(os.path.join(args.outdir, "breakpoints_with_events.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(fields + ["event_id", "event_class", "n_species_in_event", "weight"])
        for i, r in enumerate(rows):
            e = ev_of[i]
            w.writerow([r[f] for f in fields] + [e[0], e[1], e[2], round(e[3], 5)])

    # ---------------- report ----------------
    lines = []
    lines.append(f"records={len(rows)}  events={len(ev_rows)}  species={n_sp}  tol={args.tol_mb} Mb")
    for typ in ("inversion", "translocation"):
        n_rec = sum(1 for r in rows if r["type"] == typ)
        n_ev = sum(1 for e in ev_rows if e[1] == typ)
        lines.append(f"  {typ:14s} records={n_rec:4d} -> independent events={n_ev:4d}")
    lines.append("\nevents by class (records in brackets):")
    for cls in ("shared_all", "shared_subset", "private"):
        for typ in ("inversion", "translocation"):
            ev = [e for e in ev_rows if e[10] == cls and e[1] == typ]
            lines.append(f"  {cls:14s} {typ:14s} events={len(ev):4d}  records={sum(e[6] for e in ev):4d}")
    lines.append("\nrecords per species, and how many are private to it:")
    for sp in species_all:
        n = sum(1 for r in rows if r["query"] == sp)
        pr = sum(1 for i, r in enumerate(rows) if r["query"] == sp and ev_of[i][1] == "private")
        lines.append(f"  {sp:16s} records={n:3d}  private={pr:3d}")
    lines.append("\nTRANSLOCATION events (the ones worth reading by hand):")
    lines.append(f"  {'id':5s}{'ref_1':>16s}{'ref_2':>16s}{'rec':>5s}{'sp':>4s}  species / groups / class")
    for e in ev_rows:
        if e[1] == "translocation":
            lines.append(f"  {e[0]:5s}{e[2] + ':' + str(round(e[3] / 1e6, 1)) + 'Mb':>16s}"
                         f"{e[4] + ':' + str(round(e[5] / 1e6, 1)) + 'Mb':>16s}{e[6]:>5d}{e[7]:>4d}  "
                         f"{e[8]} / {e[9]} / {e[10]}")
    lines.append("\nEvents seen in >=2 species that are restricted to ONE genome group:")
    for grp in ("BBCC", "CCDD"):
        ev = [e for e in ev_rows if e[7] >= 2 and e[9] == grp]
        lines.append(f"  {grp}-only: {len(ev)} events " +
                     ", ".join(f"{e[0]}({e[1][:3]} {e[2]}:{round(e[3]/1e6,1)}Mb)" for e in ev[:12]))
    lines.append("  -> A BBCC-only translocation between the B and C subgenomes would be the "
                 "positive control (Chr1BB/Chr1CC in O. minuta and O. malampuzhaensis).")
    rep = "\n".join(lines)
    open(os.path.join(args.outdir, "events_report.txt"), "w").write(rep + "\n")
    print(rep)
    print(f"\nFiles in {args.outdir}: breakpoint_events.tsv, breakpoints_with_events.tsv, events_report.txt")


if __name__ == "__main__":
    main()
