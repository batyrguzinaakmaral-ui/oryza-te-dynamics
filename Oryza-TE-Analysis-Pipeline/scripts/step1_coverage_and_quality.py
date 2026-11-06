#!/usr/bin/env python3
"""
Step 1.1 + 1.2 — Oryza TE Annotation Dataset & Atlas
======================================================
Задача: понять, какие из 198 геномов /xdisk/dcopetti/mdalam/rte уже имеют
TE-аннотацию (panEDTA v3/v4 + TEsorter domains, 40 меток из rv2), отфильтровать
по completeness (poales_odb12), и предложить кандидатов на "representative set"
для дипломного deliverable (TE Annotation Dataset & Atlas).

Запускать НА HPC (там, где смонтирован /xdisk/dcopetti/mdalam/rte), python3
с pandas. Если pandas не установлен:
    pip install pandas --break-system-packages
    (или module load python / conda activate <env> — в зависимости от кластера)

Использование:
    python3 step1_coverage_and_quality.py --base /xdisk/dcopetti/mdalam/rte

Выход (в --outdir, по умолчанию ./te_atlas_step1/):
    coverage_full.tsv          - полная таблица 198 геномов + флаги/метрики
    coverage_by_type.tsv       - агрегат по genome-type (сколько покрыто te/)
    candidates_per_species.tsv - лучший геном на (genome_type, species)
    representative_set.tsv     - финальный шорт-лист кандидатов
    step1_report.txt           - текстовая сводка для быстрого чтения
"""

import argparse
import glob
import os
import re
import sys
from pathlib import Path

import pandas as pd


# --------------------------------------------------------------------------
# Утилиты поиска колонок: схема oryza_genomes.tsv (198 x 81) не зафиксирована
# в README дословно, поэтому ищем по шаблонам, а не по жёстким именам.
# --------------------------------------------------------------------------

def find_col(df: pd.DataFrame, patterns, required=True, label=""):
    """Найти первую колонку, чьё имя (в нижнем регистре) содержит один из
    паттернов. Возвращает имя колонки или None."""
    cols_lower = {c: c.lower() for c in df.columns}
    for pat in patterns:
        for orig, low in cols_lower.items():
            if pat in low:
                return orig
    if required:
        print(f"[ОШИБКА] Не найдена колонка для '{label}'. "
              f"Искали паттерны: {patterns}", file=sys.stderr)
        print(f"Доступные колонки:\n{list(df.columns)}", file=sys.stderr)
        sys.exit(1)
    return None


def parse_compleasm_summary(path: Path):
    """Парсит <basename>.poales_odb12.summary.txt -> dict с S/D/F/I/M/N.
    Формат compleasm summary обычно вида:
        S:95.50%, 6001
        D:2.10%, 132
        F:1.00%, 63
        I:0.00%, 0
        M:1.40%, 88
        N:6282
    Но допускаем вариации (S:xx.xx%,N:...  или построчно 'Single ...').
    """
    result = {"S_pct": None, "D_pct": None, "F_pct": None,
              "I_pct": None, "M_pct": None, "N_total": None}
    try:
        text = path.read_text(errors="replace")
    except Exception as e:
        return result

    # Ищем проценты по каждой букве категории
    for key, letter in [("S_pct", "S"), ("D_pct", "D"), ("F_pct", "F"),
                         ("I_pct", "I"), ("M_pct", "M")]:
        m = re.search(rf"\b{letter}\s*[:=]\s*([\d.]+)\s*%", text)
        if m:
            result[key] = float(m.group(1))

    m = re.search(r"\bN\s*[:=]\s*(\d+)", text)
    if m:
        result["N_total"] = int(m.group(1))

    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True,
                     help="Путь к /xdisk/dcopetti/mdalam/rte (или локальная копия)")
    ap.add_argument("--outdir", default="./te_atlas_step1",
                     help="Куда писать результаты (default: ./te_atlas_step1)")
    ap.add_argument("--min-completeness", type=float, default=90.0,
                     help="Минимальный %% Single+Duplicated (poales_odb12) "
                          "для прохождения фильтра качества (default: 90.0)")
    args = ap.parse_args()

    base = Path(args.base)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genomes_tsv = base / "oryza_genomes.tsv"
    crossref_tsv = base / "rv2_crossref.tsv"
    seqmap_tsv = base / "rv2_seqname_map.tsv"
    cm_dir = base / "cm"

    for p in [genomes_tsv, crossref_tsv, cm_dir]:
        if not p.exists():
            print(f"[ОШИБКА] Не найден путь: {p}", file=sys.stderr)
            sys.exit(1)

    # ----------------------------------------------------------------------
    # 1) Загрузка основной таблицы геномов
    # ----------------------------------------------------------------------
    genomes = pd.read_csv(genomes_tsv, sep="\t", dtype=str, low_memory=False)
    print(f"[OK] oryza_genomes.tsv: {genomes.shape[0]} строк, "
          f"{genomes.shape[1]} колонок")

    col_accession = find_col(genomes, ["accession", "assembly_accession", "gca"],
                              label="accession")
    col_species = find_col(genomes, ["species", "organism"], label="species")
    col_genome_type = find_col(genomes, ["genome_type", "genometype", "geno_type"],
                                label="genome_type")
    col_basename = find_col(genomes, ["basename", "file_base", "filename"],
                             required=False, label="basename")

    # completeness columns (poales_*); может быть несколько (S/D/F/I/M/N pct)
    poales_cols = [c for c in genomes.columns if "poales" in c.lower()]
    if not poales_cols:
        print("[ПРЕДУПРЕЖДЕНИЕ] Не найдено колонок с 'poales' в "
              "oryza_genomes.tsv — completeness возьмём только из cm/*.summary.txt",
              file=sys.stderr)

    # Если basename не найден — попробуем собрать его из accession+species+type
    # (см. схему в README: <accession>_<species>_<type>.fa.gz)
    if col_basename is None:
        def guess_basename(row):
            sp = str(row[col_species]).replace(" ", "_")
            return f"{row[col_accession]}_{sp}_{row[col_genome_type]}"
        genomes["_guessed_basename"] = genomes.apply(guess_basename, axis=1)
        col_basename = "_guessed_basename"
        print("[INFO] basename не найден в таблице — сконструирован из "
              "accession_species_genometype (может не совпасть с реальными "
              "именами файлов для rv2_ записей — проверить вручную).")

    # ----------------------------------------------------------------------
    # 2) rv2_crossref.tsv — какие 40 меток уже покрыты TE-аннотацией в te/
    # ----------------------------------------------------------------------
    crossref = pd.read_csv(crossref_tsv, sep="\t", dtype=str, low_memory=False)
    print(f"[OK] rv2_crossref.tsv: {crossref.shape[0]} строк")
    print(f"     колонки: {list(crossref.columns)}")

    # Пытаемся найти колонку с accession в crossref для джойна с genomes
    cross_acc_col = find_col(crossref, ["accession", "gca", "ncbi"],
                              required=False, label="accession (в crossref)")
    cross_label_col = find_col(crossref, ["rv2", "label"], required=False,
                                label="rv2_label (в crossref)")

    if cross_acc_col:
        covered_accessions = set(crossref[cross_acc_col].dropna().unique())
    else:
        covered_accessions = set()
        print("[ПРЕДУПРЕЖДЕНИЕ] Не нашли колонку accession в rv2_crossref.tsv — "
              "сверься с реальными именами колонок (см. вывод выше) и подправь "
              "cross_acc_col вручную в скрипте.", file=sys.stderr)

    genomes["has_te_annotation"] = genomes[col_accession].isin(covered_accessions)

    # rv2_only геномы (5 шт, без NCBI accession) — у них TE/gene аннотации
    # применяются напрямую (см. README), считаем их тоже покрытыми, если
    # basename начинается с rv2_
    rv2_mask = genomes[col_basename].astype(str).str.startswith("rv2_")
    genomes.loc[rv2_mask, "has_te_annotation"] = True

    n_covered = genomes["has_te_annotation"].sum()
    print(f"[OK] Геномов с существующей TE-аннотацией (te/): {n_covered} / "
          f"{len(genomes)}")

    # ----------------------------------------------------------------------
    # 3) Completeness из cm/*.summary.txt (Step 1.2 — фильтр качества)
    # ----------------------------------------------------------------------
    summary_files = list(cm_dir.glob("*.poales_odb12.summary.txt"))
    print(f"[OK] Найдено {len(summary_files)} файлов "
          f"*.poales_odb12.summary.txt в cm/")

    # basename файла summary = <basename>.poales_odb12.summary.txt
    summary_map = {}
    for f in summary_files:
        bn = f.name.replace(".poales_odb12.summary.txt", "")
        summary_map[bn] = parse_compleasm_summary(f)

    def lookup_completeness(basename):
        d = summary_map.get(str(basename))
        if d is None:
            return pd.Series({"S_pct": None, "D_pct": None, "F_pct": None,
                               "I_pct": None, "M_pct": None, "N_total": None,
                               "SD_pct": None, "summary_found": False})
        sd = None
        if d["S_pct"] is not None and d["D_pct"] is not None:
            sd = round(d["S_pct"] + d["D_pct"], 2)
        d2 = dict(d)
        d2["SD_pct"] = sd
        d2["summary_found"] = True
        return pd.Series(d2)

    comp_df = genomes[col_basename].apply(lookup_completeness)
    genomes = pd.concat([genomes, comp_df], axis=1)

    n_summary_found = genomes["summary_found"].sum()
    print(f"[OK] Completeness-файл найден и распарсен для {n_summary_found} / "
          f"{len(genomes)} геномов "
          f"(остальные — basename не совпал, см. предупреждение выше "
          f"если использовалась эвристика _guessed_basename).")

    # Флаг качества
    genomes["passes_quality"] = genomes["SD_pct"].apply(
        lambda x: (x is not None) and (x >= args.min_completeness)
    )

    # ----------------------------------------------------------------------
    # 4) Агрегация по genome-type
    # ----------------------------------------------------------------------
    by_type = genomes.groupby(col_genome_type).agg(
        n_total=(col_accession, "count"),
        n_with_te=("has_te_annotation", "sum"),
        n_pass_quality=("passes_quality", "sum"),
        n_species=(col_species, "nunique"),
        mean_SD_pct=("SD_pct", "mean"),
    ).reset_index().sort_values("n_total", ascending=False)

    # ----------------------------------------------------------------------
    # 5) Кандидаты: лучший геном на (genome_type, species)
    #    приоритет: has_te_annotation=True, затем максимальный SD_pct
    # ----------------------------------------------------------------------
    genomes_sorted = genomes.sort_values(
        by=["has_te_annotation", "SD_pct"],
        ascending=[False, False]
    )
    candidates = genomes_sorted.groupby([col_genome_type, col_species],
                                         as_index=False).first()

    keep_cols = [col_accession, col_species, col_genome_type, col_basename,
                 "has_te_annotation", "SD_pct", "passes_quality",
                 "summary_found"]
    keep_cols = [c for c in keep_cols if c in candidates.columns]
    candidates_out = candidates[keep_cols].sort_values(
        [col_genome_type, "has_te_annotation", "SD_pct"],
        ascending=[True, False, False]
    )

    # ----------------------------------------------------------------------
    # 6) Финальный representative_set: приоритет te/-покрытым + passes_quality
    #    Если для genome-type НЕТ ни одного te/-покрытого генома — geno-type
    #    помечается как "REQUIRES_NEW_ANNOTATION" (кандидат на пропуск при
    #    сжатых сроках до ноября, если только это не критичный тип, напр. AA).
    # ----------------------------------------------------------------------
    rep_rows = []
    for gtype, grp in candidates_out.groupby(col_genome_type):
        covered = grp[grp["has_te_annotation"]]
        if len(covered) > 0:
            best = covered.sort_values("SD_pct", ascending=False).iloc[0]
            status = "READY_te_dir"
        else:
            best = grp.sort_values("SD_pct", ascending=False).iloc[0] \
                if len(grp) > 0 else None
            status = "REQUIRES_NEW_panEDTA_TEsorter"
        if best is not None:
            row = best.to_dict()
            row["status"] = status
            row["n_species_in_type"] = grp[col_species].nunique()
            rep_rows.append(row)

    representative_set = pd.DataFrame(rep_rows).sort_values(
        ["status", col_genome_type], ascending=[True, True]
    )

    # ----------------------------------------------------------------------
    # 7) Запись результатов
    # ----------------------------------------------------------------------
    genomes.to_csv(outdir / "coverage_full.tsv", sep="\t", index=False)
    by_type.to_csv(outdir / "coverage_by_type.tsv", sep="\t", index=False)
    candidates_out.to_csv(outdir / "candidates_per_species.tsv", sep="\t",
                           index=False)
    representative_set.to_csv(outdir / "representative_set.tsv", sep="\t",
                               index=False)

    # ----------------------------------------------------------------------
    # 8) Текстовый отчёт
    # ----------------------------------------------------------------------
    lines = []
    lines.append("=" * 70)
    lines.append("STEP 1.1 + 1.2 — Coverage & Quality Report")
    lines.append("=" * 70)
    lines.append(f"Всего геномов в oryza_genomes.tsv: {len(genomes)}")
    lines.append(f"Уже покрыто TE-аннотацией (te/, 40 rv2-меток): {n_covered}")
    lines.append(f"Completeness-файл найден: {n_summary_found}")
    lines.append(f"Порог качества (S+D poales_odb12): >= {args.min_completeness}%")
    lines.append("")
    lines.append("-- По genome-type --")
    lines.append(by_type.to_string(index=False))
    lines.append("")
    lines.append("-- Genome-types БЕЗ существующей TE-аннотации "
                  "(потребуют новый прогон panEDTA+TEsorter) --")
    missing_types = by_type[by_type["n_with_te"] == 0][col_genome_type].tolist()
    if missing_types:
        lines.append(", ".join(missing_types))
    else:
        lines.append("(нет — все типы имеют хотя бы один te/-покрытый геном)")
    lines.append("")
    lines.append("-- Итоговый representative set (см. representative_set.tsv) --")
    lines.append(representative_set.to_string(index=False))
    lines.append("")
    lines.append("РЕКОМЕНДАЦИЯ при дедлайне до ноября:")
    lines.append("  1) Из representative_set.tsv взять все строки со статусом "
                  "READY_te_dir — это можно сразу вести в TEsorter-парсинг "
                  "и стандартизацию номенклатуры (Wicker), без затрат времени "
                  "на новый прогон panEDTA.")
    lines.append("  2) Строки REQUIRES_NEW_panEDTA_TEsorter обсудить с "
                  "Зукколо/Вингом: включать ли их в финальный atlas, или "
                  "явно ограничить scope дипломной работы уже аннотированными "
                  "40 rv2-геномами (это самодостаточный и полностью "
                  "воспроизводимый deliverable).")
    report_text = "\n".join(lines)
    (outdir / "step1_report.txt").write_text(report_text)
    print("\n" + report_text)
    print(f"\n[OK] Все файлы сохранены в: {outdir.resolve()}")


if __name__ == "__main__":
    main()
