#!/usr/bin/env python3
"""
step6_te_at_breakpoints.py

Tests H1/H2/H3: is TE content higher in rearrangement-breakpoint intervals
than in comparable control regions?

Input
  breakpoints.tsv (step5)   intervals on the QUERY genome (gene-to-gene gap)
  edta/<genome>/<genome>.fa.mod.EDTA.TEanno.gff3   all TE calls
  edta/<genome>/<genome>.fa.mod.EDTA.intact.gff3   structurally intact copies
  genome_seq_lengths.tsv (step4)

Categories (union coverage, so overlapping calls are NOT double counted):
  all     every feature with a classification=
  LTR     classification starts with LTR            (Class I)
  DNA     classification starts with DNA or MITE    (Class II, incl. Helitron)
  intact  every feature in the intact.gff3          (recently active proxy)
  NOTE on intact: sub-features (LTRs, TSDs) nest inside or sit next to the element,
  so the union is the element footprint plus a few bp of TSD. Negligible at 1 kb bins.

Observed value = mean TE density (bp covered / bp) over the breakpoint intervals
of a scope. Two null models, each resampling ONE region per breakpoint, of the
same width, --perms times:
  chrom   uniform on the same chromosome
  local   uniform within +-local-mb of the breakpoint (excluding it) -- controls
          for the centromere-to-arm TE gradient, the stricter and preferred test.
  gap     (needs --orthologs-dir) resamples with replacement the OTHER inter-ortholog
          gaps of the same chromosome within +-local-mb whose width is within
          0.5x-2x of the breakpoint's. Breakpoint intervals lie BETWEEN two BUSCO
          genes, i.e. by construction they exclude the flanking gene bodies, which
          carry almost no TE; random windows (chrom/local) can overlap genes and so
          look TE-poorer -> spurious "enrichment" of ~ the gene-body fraction (~5-15%).
          The gap null removes this bias and is the one to trust.
Reported: ratio observed / null mean and one-sided empirical p (enrichment, depletion).

Pseudo-replication: the same rearrangement can be reported in several species and on both
homeologous chromosomes of one species (e.g. a junction that is really specific to the
O. officinalis reference). step7 clusters such records into events; --use-weights gives each
event total weight 1 and --drop-class shared_all removes events seen in every species.

Statistical caveats (put these in the thesis):
  * The breakpoint is an interval between two flanking BUSCO genes (~100-400 kb),
    not a position; a real TE hotspot of a few kb is diluted by the interval width.
  * Breakpoints within a species are not fully independent (nearby, nested events).
  * Few breakpoints per scope -> low power; non-significant != no effect.
  * Breakpoints shared by ALL five species at the same officinalis position are
    probably reference-lineage events; they are counted once per species here.

Usage (run from /xdisk/dcopetti/mdalam/rte):
    python3 step6_te_at_breakpoints.py \
        --breakpoints dual_students_te_oryza_nku/data/te_atlas_step5_mg40/breakpoints.tsv \
        --lens dual_students_te_oryza_nku/data/te_atlas_step4/genome_seq_lengths.tsv \
        --orthologs-dir dual_students_te_oryza_nku/data/te_atlas_step4 \
        --edta-dir edta --outdir dual_students_te_oryza_nku/data/te_atlas_step6_mg40
Stdlib only.
"""

import argparse
import csv
import glob
import os
import random
import re
from array import array
from collections import defaultdict

BIN = 1000
CLS = re.compile(r"classification=([^;\s]+)")


def load_intervals(path, intact):
    out = defaultdict(lambda: defaultdict(list))   # cat -> chrom -> [(s, e)]
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 9:
                continue
            m = CLS.search(c[8])
            if not m:
                continue
            try:
                s, e = int(c[3]), int(c[4])
            except ValueError:
                continue
            if e < s:
                s, e = e, s
            cls = m.group(1)
            if intact:
                out["intact"][c[0]].append((s, e))
            else:
                out["all"][c[0]].append((s, e))
                if cls.startswith("LTR"):
                    out["LTR"][c[0]].append((s, e))
                elif cls.startswith(("DNA", "MITE")):
                    out["DNA"][c[0]].append((s, e))
    return out


def build_prefix(intervals, nbins):
    merged = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1] + 1:
            if e > merged[-1][1]:
                merged[-1][1] = e
        else:
            merged.append([s, e])
    bins = array("d", [0.0]) * nbins
    for s, e in merged:
        s0 = s - 1
        b0, b1 = s0 // BIN, min((e - 1) // BIN, nbins - 1)
        for b in range(b0, b1 + 1):
            lo, hi = max(s0, b * BIN), min(e, (b + 1) * BIN)
            if hi > lo:
                bins[b] += hi - lo
    pref = array("d", [0.0]) * (nbins + 1)
    acc = 0.0
    for i in range(nbins):
        acc += bins[i]
        pref[i + 1] = acc
    return pref


def dens(pref, lo, w):
    return (pref[lo + w] - pref[lo]) / (w * BIN)


def sample_local(lo, w, nbins, span, n, rng):
    a, b = max(0, lo - span), min(nbins - w, lo + span)
    loose = (b - a) < 2 * w
    res = []
    while len(res) < n:
        x = rng.randint(a, b)
        if loose or abs(x - lo) >= w:
            res.append(x)
    return res


def load_gaps(orth_dir, species, min_gap=2000):
    """chrom -> sorted list of (start_bp, end_bp) gaps between consecutive orthologs."""
    path = os.path.join(orth_dir, f"busco_orthologs_{species}_vs_officinalis.tsv")
    genes = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            genes[r["q_seq"]].append((int(r["q_start"]), int(r["q_end"])))
    gaps = {}
    for ch, g in genes.items():
        g.sort()
        out, maxend = [], g[0][1]
        for s, e in g[1:]:
            if s - maxend >= min_gap:
                out.append((maxend, s))
            maxend = max(maxend, e)
        gaps[ch] = out
    return gaps


def mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--breakpoints", required=True)
    ap.add_argument("--lens", required=True)
    ap.add_argument("--edta-dir", default="edta")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--orthologs-dir", default=None,
                    help="step4 dir with busco_orthologs_*_vs_officinalis.tsv (enables the gap null)")
    ap.add_argument("--use-weights", action="store_true",
                    help="weight each breakpoint by the 'weight' column (step7 output: 1/size of its "
                         "event cluster) so that one event counts once however many species/chromosomes show it")
    ap.add_argument("--drop-class", nargs="*", default=[],
                    help="drop breakpoints whose 'event_class' (step7) is in this list, e.g. shared_all")
    ap.add_argument("--perms", type=int, default=5000)
    ap.add_argument("--local-mb", type=float, default=3.0)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = random.Random(args.seed)
    span = int(args.local_mb * 1e6 / BIN)

    lens = defaultdict(dict)
    with open(args.lens) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            lens[r["genome"]][r["seqid"]] = int(r["length"])

    bps = defaultdict(list)
    n_dropped = 0
    with open(args.breakpoints) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if args.drop_class:
                if "event_class" not in r:
                    raise SystemExit("--drop-class needs the step7 output (breakpoints_with_events.tsv)")
                if r["event_class"] in args.drop_class:
                    n_dropped += 1
                    continue
            if args.use_weights and "weight" not in r:
                raise SystemExit("--use-weights needs the step7 output (breakpoints_with_events.tsv)")
            bps[r["query"]].append(r)
    if args.drop_class:
        print(f"[INFO] dropped {n_dropped} breakpoints of class {args.drop_class}")

    cats = ["all", "LTR", "DNA", "intact"]
    records = []

    for species, rows in sorted(bps.items()):
        genome = next((g for g in lens if f"_Oryza_{species}_" in g), None)
        if genome is None:
            print(f"[WARN] no genome for {species}")
            continue
        base = os.path.join(args.edta_dir, genome, f"{genome}.fa.mod.EDTA.")
        gff_path, intact_path = base + "TEanno.gff3", base + "intact.gff3"
        if not os.path.isfile(gff_path):
            # fallback: look for the species name anywhere under --edta-dir
            hits = [f for f in glob.glob(os.path.join(args.edta_dir, "**", "*.TEanno.gff3"), recursive=True)
                    if f"_Oryza_{species}_" in f]
            if not hits:
                raise SystemExit(f"[ERROR] TEanno.gff3 for {species} not found under '{args.edta_dir}' "
                                 f"(current dir: {os.getcwd()}). Run from rte/ or pass --edta-dir with the right path, "
                                 f"e.g. ../../edta when running from dual_students_te_oryza_nku/data")
            gff_path = hits[0]
            intact_path = gff_path.replace(".TEanno.gff3", ".intact.gff3")
        print(f"[{species}] parsing TE annotation from {gff_path} ({len(rows)} breakpoints)...")
        ivs = load_intervals(gff_path, intact=False)
        if os.path.isfile(intact_path):
            ivs.update(load_intervals(intact_path, intact=True))
        else:
            print(f"[WARN] {species}: no intact.gff3 at {intact_path}; 'intact' category will be empty")

        pref = {}
        for ch in {r["q_chr"] for r in rows}:
            if ch not in lens[genome]:
                continue
            nb = lens[genome][ch] // BIN + 1
            for cat in cats:
                pref[(cat, ch)] = build_prefix(ivs[cat].get(ch, []), nb)
        if not any(k[0] == "all" for k in pref):
            print(f"[WARN] {species}: no chromosome ids matched TE annotation seqids")
        gaps = load_gaps(args.orthologs_dir, species) if args.orthologs_dir else {}
        n_used = 0
        for r in rows:
            ch = r["q_chr"]
            if ("all", ch) not in pref:
                continue
            nbins = len(pref[("all", ch)]) - 1
            lo = int(r["bp_start"]) // BIN
            hi = -(-int(r["bp_end"]) // BIN)
            w = max(1, hi - lo)
            if lo + w > nbins or nbins - w < 1:
                continue
            rec = {"species": species, "q_chr": ch, "type": r["type"],
                   "start": int(r["bp_start"]), "end": int(r["bp_end"]),
                   "width": int(r["bp_end"]) - int(r["bp_start"]),
                   "obs": {}, "chrom": {}, "local": {}, "gap": {},
                   "wt": float(r["weight"]) if args.use_weights else 1.0}
            bs, be = int(r["bp_start"]), int(r["bp_end"])
            cand = []
            if gaps:
                mid, wbp = (bs + be) / 2, max(be - bs, 1)
                for mult in (1, 2, 4):   # widen the neighbourhood if too few candidates
                    cand = [(g0 // BIN, max(1, -(-g1 // BIN) - g0 // BIN))
                            for g0, g1 in gaps.get(ch, [])
                            if abs((g0 + g1) / 2 - mid) <= span * BIN * mult
                            and (g1 <= bs or g0 >= be)
                            and 0.5 * wbp <= (g1 - g0) <= 2 * wbp]
                    if len(cand) >= 5:
                        break
                if len(cand) < 5:
                    cand = []
            pos_chrom = [rng.randint(0, nbins - w) for _ in range(args.perms)]
            pos_local = sample_local(lo, w, nbins, span, args.perms, rng)
            for cat in cats:
                p = pref[(cat, ch)]
                rec["obs"][cat] = dens(p, lo, w)
                rec["chrom"][cat] = [dens(p, x, w) for x in pos_chrom]
                rec["local"][cat] = [dens(p, x, w) for x in pos_local]
                if cand:
                    cd = [dens(p, lo2, w2) for lo2, w2 in cand if lo2 + w2 <= nbins]
                    if cd:
                        rec["gap"][cat] = [rng.choice(cd) for _ in range(args.perms)]
            records.append(rec)
            n_used += 1
        print(f"    used {n_used}/{len(rows)} breakpoints")
        del ivs, pref

    # ---- per-breakpoint table ----------------------------------------------------
    with open(os.path.join(args.outdir, "te_at_breakpoints_per_bp.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["species", "q_chr", "type", "bp_start", "bp_end", "width"]
                   + [f"obs_{c}" for c in cats] + [f"local_null_{c}" for c in cats])
        for r in records:
            w.writerow([r["species"], r["q_chr"], r["type"], r["start"], r["end"], r["width"]]
                       + [round(r["obs"][c], 5) for c in cats]
                       + [round(mean(r["local"][c]), 5) for c in cats])

    # ---- summaries -------------------------------------------------------------------
    scopes = sorted({r["species"] for r in records}) + ["POOLED"]
    types = ["all", "inversion", "translocation"]
    out_rows = []
    B = args.perms
    for scope in scopes:
        for typ in types:
            sel = [r for r in records
                   if (scope == "POOLED" or r["species"] == scope)
                   and (typ == "all" or r["type"] == typ)]
            if not sel:
                continue
            for cat in cats:
                row = [scope, typ, cat, len(sel)]
                for mode in ("chrom", "local", "gap"):
                    sm = [r for r in sel if r[mode]]
                    if not sm:
                        row += ["NA"] * 6
                        continue
                    ws = [r["wt"] for r in sm]
                    sw = sum(ws)
                    obs = sum(wt * r["obs"][cat] for wt, r in zip(ws, sm)) / sw
                    cols = zip(*[r[mode][cat] for r in sm])
                    nm = [sum(wt * x for wt, x in zip(ws, c)) / sw for c in cols]
                    avg = mean(nm)
                    pe = (1 + sum(m >= obs for m in nm)) / (B + 1)
                    pd = (1 + sum(m <= obs for m in nm)) / (B + 1)
                    row += [len(sm), round(obs, 5), round(avg, 5),
                            round(obs / avg, 3) if avg else "NA", round(pe, 4), round(pd, 4)]
                row.append(round(sum(r["wt"] for r in sel), 1))   # weighted n = number of events when pooled
                out_rows.append(row)
    hdr = ["scope", "bp_type", "category", "n_breakpoints"]
    for m in ("chrom", "local", "gap"):
        hdr += [f"{m}_n", f"{m}_obs", f"{m}_null", f"{m}_ratio", f"{m}_p_enrich", f"{m}_p_deplete"]
    hdr += ["n_weighted"]
    with open(os.path.join(args.outdir, "te_at_breakpoints_summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr)
        w.writerows(out_rows)

    print("\n=== POOLED across species ===  (gap null = the one to trust; n_gap may be < n)")
    print(f"{'type':14s}{'cat':8s}{'n':>4s}{'n_w':>7s}{'local_ratio':>13s}{'local_p':>9s}"
          f"{'n_gap':>7s}{'gap_ratio':>11s}{'gap_p_enr':>11s}{'gap_p_dep':>11s}")
    for r in out_rows:
        if r[0] == "POOLED":
            print(f"{r[1]:14s}{r[2]:8s}{r[3]:>4d}{r[22]:>7.1f}{r[13]!s:>13s}{r[14]!s:>9s}"
                  f"{r[16]!s:>7s}{r[19]!s:>11s}{r[20]!s:>11s}{r[21]!s:>11s}")
    if args.use_weights:
        print("(weighted: each rearrangement event counts once; n_w = number of independent events)")
    print(f"\nFiles in {args.outdir}: te_at_breakpoints_summary.tsv, te_at_breakpoints_per_bp.tsv")


if __name__ == "__main__":
    main()
