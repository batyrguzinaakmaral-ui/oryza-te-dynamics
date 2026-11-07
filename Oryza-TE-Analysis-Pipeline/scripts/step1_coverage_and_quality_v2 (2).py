#!/usr/bin/env python3
"""
Step 1.1 + 1.2 (v2) — Oryza TE Annotation Dataset & Atlas
===========================================================
Версия под РЕАЛЬНУЮ схему файлов (проверено по head -n6 на HPC):

  oryza_genomes.tsv содержит НАПРЯМУЮ:
    accession, genome_type, species, ...
    poales_single, poales_dup, poales_frag, poales_incomplete,
    poales_missing, poales_total, poales_complete_pct   <- completeness готов,
                                                             парсить cm/*.txt не нужно
    rv2_redundant, rv2_label, rv2_note                  <- мост к te/-покрытию

  rv2_crossref.tsv содержит:
    rv2_label, rv2_n_seqs, rv2_total_length,
    rv2_has_te_annotation, rv2_subgenome,
    rte_accession, rte_organism, rte_assembly_name,
    match_type, note

Ключ джойна между таблицами — rv2_label (текстовое совпадение), а не
accession, т.к. rv2_label уже присутствует в обеих таблицах напрямую.

Запуск на HPC:
    cd /xdisk/dcopetti/mdalam/rte/dual_students_te_oryza_nku
    python3 step1_coverage_and_quality_v2.py \
        --base /xdisk/dcopetti/mdalam/rte \
        --outdir ./data/te_atlas_step1

Выход (в --outdir):
    coverage_full.tsv          - все 198 геномов + has_te_annotation + completeness
    coverage_by_type.tsv       - агрегат по genome_type
    candidates_per_species.tsv - лучший геном на (genome_type, species)
    representative_set.tsv     - финальный шорт-лист
    step1_report.txt           - текстовая сводка
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def to_bool(series: pd.Series) -> pd.Series:
    """Приводит колонку вида yes/no/true/false/1/0/NaN к bool."""
    truthy = {"yes", "true", "1", "y", "t"}
    return series.astype(str).str.strip().str.lower().isin(truthy)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True,
                     help="Путь к /xdisk/dcopetti/mdalam/rte")
    ap.add_argument("--outdir", default="./te_atlas_step1")
    ap.add_argument("--min-completeness", type=float, default=90.0,
                     help="Минимальный poales_complete_pct для passes_quality "
                          "(default: 90.0)")
    args = ap.parse_args()

    base = Path(args.base)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genomes_path = base / "oryza_genomes.tsv"
    crossref_path = base / "rv2_crossref.tsv"

    for p in [genomes_path, crossref_path]:
        if not p.exists():
            sys.exit(f"[ОШИБКА] Не найден файл: {p}")

    # ------------------------------------------------------------------
    # 1) Загрузка
    # ------------------------------------------------------------------
    genomes = pd.read_csv(genomes_path, sep="\t", dtype=str, low_memory=False)
    crossref = pd.read_csv(crossref_path, sep="\t", dtype=str, low_memory=False)

    print(f"[OK] oryza_genomes.tsv: {genomes.shape[0]} строк, "
          f"{genomes.shape[1]} колонок")
    print(f"[OK] rv2_crossref.tsv: {crossref.shape[0]} строк, "
          f"{crossref.shape[1]} колонок")

    required_genome_cols = ["accession", "genome_type", "species",
                             "poales_complete_pct", "rv2_label"]
    missing = [c for c in required_genome_cols if c not in genomes.columns]
    if missing:
        sys.exit(f"[ОШИБКА] В oryza_genomes.tsv нет колонок: {missing}\n"
                  f"Реальные колонки: {list(genomes.columns)}")

    required_cross_cols = ["rv2_label", "rv2_has_te_annotation", "rte_accession"]
    missing = [c for c in required_cross_cols if c not in crossref.columns]
    if missing:
        sys.exit(f"[ОШИБКА] В rv2_crossref.tsv нет колонок: {missing}\n"
                  f"Реальные колонки: {list(crossref.columns)}")

    # ------------------------------------------------------------------
    # 2) Приведение типов
    # ------------------------------------------------------------------
    genomes["poales_complete_pct"] = pd.to_numeric(
        genomes["poales_complete_pct"], errors="coerce")

    crossref["rv2_has_te_annotation"] = to_bool(crossref["rv2_has_te_annotation"])

    # ------------------------------------------------------------------
    # 3) Join по rv2_label: только строки crossref с rv2_has_te_annotation=True
    # ------------------------------------------------------------------
    te_labels = set(
        crossref.loc[crossref["rv2_has_te_annotation"], "rv2_label"].dropna()
    )
    print(f"[OK] rv2_label с rv2_has_te_annotation=True: {len(te_labels)}")

    genomes["has_te_annotation"] = (
        genomes["rv2_label"].notna()
        & (genomes["rv2_label"].str.strip() != "")
        & genomes["rv2_label"].isin(te_labels)
    )

    n_covered = int(genomes["has_te_annotation"].sum())
    print(f"[OK] Геномов с TE-аннотацией (te/): {n_covered} / {len(genomes)}")

    # Сверка: сколько rv2_label в genomes вообще не нашли соответствия
    # в crossref (опечатка/рассинхрон меток) — полезно для отладки
    genome_labels = set(genomes.loc[genomes["rv2_label"].notna()
                                     & (genomes["rv2_label"].str.strip() != ""),
                                     "rv2_label"])
    unmatched_labels = genome_labels - set(crossref["rv2_label"].dropna())
    if unmatched_labels:
        print(f"[ПРЕДУПРЕЖДЕНИЕ] {len(unmatched_labels)} rv2_label из "
              f"oryza_genomes.tsv не найдены в rv2_crossref.tsv: "
              f"{sorted(unmatched_labels)[:10]}{'...' if len(unmatched_labels) > 10 else ''}")

    # ------------------------------------------------------------------
    # 4) Фильтр качества (Step 1.2) — используем готовый poales_complete_pct
    # ------------------------------------------------------------------
    genomes["passes_quality"] = (
        genomes["poales_complete_pct"] >= args.min_completeness
    )

    # ------------------------------------------------------------------
    # 5) Агрегация по genome_type
    # ------------------------------------------------------------------
    by_type = genomes.groupby("genome_type").agg(
        n_total=("accession", "count"),
        n_with_te=("has_te_annotation", "sum"),
        n_pass_quality=("passes_quality", "sum"),
        n_species=("species", "nunique"),
        mean_completeness=("poales_complete_pct", "mean"),
    ).reset_index().sort_values("n_total", ascending=False)

    # ------------------------------------------------------------------
    # 6) Кандидаты: лучший геном на (genome_type, species)
    #    приоритет has_te_annotation=True, затем max poales_complete_pct
    # ------------------------------------------------------------------
    genomes_sorted = genomes.sort_values(
        by=["has_te_annotation", "poales_complete_pct"],
        ascending=[False, False]
    )
    candidates = genomes_sorted.groupby(["genome_type", "species"],
                                         as_index=False).first()

    keep_cols = ["accession", "species", "genome_type",
                 "has_te_annotation", "poales_complete_pct",
                 "passes_quality", "rv2_label"]
    candidates_out = candidates[keep_cols].sort_values(
        ["genome_type", "has_te_annotation", "poales_complete_pct"],
        ascending=[True, False, False]
    )

    # ------------------------------------------------------------------
    # 7) representative_set: 1 лучший геном на genome_type
    # ------------------------------------------------------------------
    rep_rows = []
    for gtype, grp in candidates_out.groupby("genome_type"):
        covered = grp[grp["has_te_annotation"]]
        if len(covered) > 0:
            best = covered.sort_values("poales_complete_pct",
                                        ascending=False).iloc[0]
            status = "READY_te_dir"
        else:
            best = grp.sort_values("poales_complete_pct",
                                    ascending=False).iloc[0]
            status = "REQUIRES_NEW_panEDTA_TEsorter"
        row = best.to_dict()
        row["status"] = status
        row["n_species_in_type"] = int(grp["species"].nunique())
        rep_rows.append(row)

    representative_set = pd.DataFrame(rep_rows).sort_values(
        ["status", "genome_type"]
    )

    # ------------------------------------------------------------------
    # 8) Запись файлов
    # ------------------------------------------------------------------
    genomes.to_csv(outdir / "coverage_full.tsv", sep="\t", index=False)
    by_type.to_csv(outdir / "coverage_by_type.tsv", sep="\t", index=False)
    candidates_out.to_csv(outdir / "candidates_per_species.tsv", sep="\t",
                           index=False)
    representative_set.to_csv(outdir / "representative_set.tsv", sep="\t",
                               index=False)

    # ------------------------------------------------------------------
    # 9) Текстовый отчёт
    # ------------------------------------------------------------------
    lines = []
    lines.append("=" * 70)
    lines.append("STEP 1.1 + 1.2 (v2) — Coverage & Quality Report")
    lines.append("=" * 70)
    lines.append(f"Всего геномов: {len(genomes)}")
    lines.append(f"Покрыто TE-аннотацией (te/): {n_covered}")
    lines.append(f"Порог качества poales_complete_pct >= {args.min_completeness}%")
    lines.append("")
    lines.append("-- По genome_type --")
    lines.append(by_type.to_string(index=False))
    lines.append("")
    missing_types = by_type.loc[by_type["n_with_te"] == 0, "genome_type"].tolist()
    lines.append("-- genome_type БЕЗ te/-покрытия (нужен новый panEDTA+TEsorter) --")
    lines.append(", ".join(missing_types) if missing_types else "(нет)")
    lines.append("")
    lines.append("-- Representative set --")
    lines.append(representative_set.to_string(index=False))
    lines.append("")
    lines.append("РЕКОМЕНДАЦИЯ при дедлайне до ноября:")
    lines.append("  READY_te_dir -> сразу в TEsorter-парсинг / стандартизацию Wicker.")
    lines.append("  REQUIRES_NEW_panEDTA_TEsorter -> обсудить с Зукколо/Вингом: "
                  "включать в scope или явно ограничиться 40 rv2-геномами.")
    report = "\n".join(lines)
    (outdir / "step1_report.txt").write_text(report)
    print("\n" + report)
    print(f"\n[OK] Результаты в: {outdir.resolve()}")


if __name__ == "__main__":
    main()
