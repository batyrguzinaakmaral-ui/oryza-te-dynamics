#!/usr/bin/env python3
"""
Step 2 (v2) — CC-clade TE composition, построено ИЗ ГОТОВЫХ агрегатов
=========================================================================
Использует уже посчитанные таблицы edta/testats/ (инструмент tecompile),
а не парсит TEanno.gff3 заново. Фокус — 9(факт. 6 найденных) видов
CC-клады (CC/BBCC/CCDD), запрошенных Michael/Rod.

Вход:
    edta/testats/genome_annot.tsv       genome, species, genome_type,
                                          mating_system, n_intact, bp_intact,
                                          genome_size   (БЕЗ заголовка!)
    edta/testats/per_genome_class.tsv   genome, classification, n_elements,
                                          total_bp, mean_len, median_len,
                                          min_len, max_len   (С заголовком)
    edta/testats/class_to_type.tsv      classification, feature_type, n
    edta/testats/across_genomes_class.tsv  genус-wide базовая статистика

Запуск:
    python3 step2_cc_clade_te_composition.py --base /xdisk/dcopetti/mdalam/rte \
        --outdir ./data/te_atlas_step2

Выход:
    cc_clade_genome_annot.tsv       - метаданные 6 найденных CC-clade геномов
    cc_clade_per_genome_class.tsv   - их TE-классы (длинный формат)
    cc_clade_class_wide_pct.tsv     - широкий формат: % генома по классам,
                                       геномы в столбцах, + genus-wide среднее
    cc_clade_report.txt             - что нашли / что отсутствует
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


# Известные CC-клада геномы (найдено вручную по ls edta/, см. переписку).
# Ключ — точный genome-label, как он записан в testats/*.tsv (это basename
# папки в edta/, БЕЗ .fa.mod... суффиксов).
CC_CLADE_GENOMES = {
    "GCA_008326285.1_Oryza_officinalis_CC":   ("CC",   "Oryza_officinalis"),
    "GCA_048166525.1_Oryza_minuta_BBCC":      ("BBCC", "Oryza_minuta"),
    "GCA_048564985.1_Oryza_malampuzhaensis_BBCC": ("BBCC", "Oryza_malampuzhaensis"),
    "GCA_047899615.1_Oryza_alta_CCDD":        ("CCDD", "Oryza_alta"),
    "GCA_048174585.1_Oryza_latifolia_CCDD":   ("CCDD", "Oryza_latifolia"),
    "GCA_048188845.1_Oryza_grandiglumis_CCDD": ("CCDD", "Oryza_grandiglumis"),
}

# Виды, которых по брифу Michael ожидалось ещё найти, но их физически нет
# в делайвери. ИСПРАВЛЕНО после проверки литературы (см. переписку):
# BBCC действительно включает 3 вида, а не 2 — третий это ТЕТРАПЛОИДНАЯ
# форма O. punctata (2n=48), отдельная от диплоидной BB-формы (2n=24).
# У неё независимое происхождение от minuta/malampuzhaensis: тетраплоидная
# punctata произошла от диплоидных punctata (B) x eichingeri (C), тогда как
# minuta/malampuzhaensis — от punctata (B) x officinalis (C).
# См.: Multiple origins of BBCC allopolyploid species (Sci Rep 2015),
#      PLOS One BBCC microsatellite markers (2014).
CC_CLADE_MISSING_NOTE = (
    "CC (диплоид, 3 признанных вида: officinalis, rhizomatis, eichingeri): "
    "в делайвери найден только O. officinalis. Сборок O. rhizomatis и "
    "O. eichingeri в /xdisk/dcopetti/mdalam/rte нет вообще (физически "
    "отсутствуют, не ошибка поиска).\n"
    "BBCC (3 признанных вида: minuta, malampuzhaensis, И тетраплоидная "
    "форма punctata): найдено 2 из 3 (minuta, malampuzhaensis). "
    "ВАЖНО — уточнено после проверки литературы: третий вид, тетраплоидная "
    "O. punctata (2n=48), это НЕ то же самое, что диплоидная O. punctata "
    "(2n=24, тип BB), которая уже есть в коллекции (GCA_000573905.2). Это "
    "разные цитотипы одного вида с независимым происхождением. Нужно "
    "отдельно проверить, есть ли в делайвери именно тетраплоидная сборка "
    "punctata (скорее всего её пока нет — увидели только BB-форму).\n"
    "CCDD (3 признанных вида): полный триплет найден (O. alta, O. "
    "latifolia, O. grandiglumis)."
)


def read_no_header_tsv(path: Path, fieldnames):
    rows = []
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            vals = line.split("\t")
            rows.append(dict(zip(fieldnames, vals)))
    return rows


def read_header_tsv(path: Path):
    with open(path, newline="", encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader), reader.fieldnames


def write_tsv(path: Path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                            extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="/xdisk/dcopetti/mdalam/rte")
    ap.add_argument("--outdir", default="./te_atlas_step2")
    args = ap.parse_args()

    base = Path(args.base)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    testats = base / "edta" / "testats"

    genome_annot_path = testats / "genome_annot.tsv"
    per_genome_class_path = testats / "per_genome_class.tsv"
    class_to_type_path = testats / "class_to_type.tsv"
    across_genomes_path = testats / "across_genomes_class.tsv"

    for p in [genome_annot_path, per_genome_class_path,
              class_to_type_path, across_genomes_path]:
        if not p.exists():
            sys.exit(f"[ОШИБКА] Не найден файл: {p}")

    # ------------------------------------------------------------------
    # 1) Загрузка
    # ------------------------------------------------------------------
    genome_annot_fields = ["genome", "species", "genome_type",
                            "mating_system", "n_intact_elements",
                            "total_bp_intact", "genome_size"]
    genome_annot = read_no_header_tsv(genome_annot_path, genome_annot_fields)
    genome_annot_by_id = {r["genome"]: r for r in genome_annot}
    print(f"[OK] genome_annot.tsv: {len(genome_annot)} геномов")

    per_class, per_class_fields = read_header_tsv(per_genome_class_path)
    print(f"[OK] per_genome_class.tsv: {len(per_class)} строк")

    class_to_type, _ = read_header_tsv(class_to_type_path)
    class_to_type_map = {r["classification"]: r["feature_type"]
                          for r in class_to_type}
    print(f"[OK] class_to_type.tsv: {len(class_to_type)} классов")

    across_genomes, _ = read_header_tsv(across_genomes_path)
    across_by_class = {r["classification"]: r for r in across_genomes}
    print(f"[OK] across_genomes_class.tsv: {len(across_genomes)} классов "
          f"(genus-wide baseline)")

    # ------------------------------------------------------------------
    # 2) Проверка наличия CC-clade геномов в данных
    # ------------------------------------------------------------------
    found = []
    missing = []
    for gid, (gtype, species) in CC_CLADE_GENOMES.items():
        if gid in genome_annot_by_id:
            found.append(gid)
        else:
            missing.append(gid)

    print(f"\n[OK] Найдено CC-clade геномов в genome_annot.tsv: "
          f"{len(found)} / {len(CC_CLADE_GENOMES)}")
    if missing:
        print(f"[ПРЕДУПРЕЖДЕНИЕ] НЕ найдены в testats/ (проверь написание "
              f"или отсутствие EDTA-прогона): {missing}")

    # Отдельная проверка: есть ли в коллекции ТЕТРАПЛОИДНАЯ O. punctata
    # (BBCC), а не только уже известная диплоидная BB-форма. Ищем по
    # подстроке "punctata" во всех genome_id, независимо от жёсткого
    # списка CC_CLADE_GENOMES выше.
    punctata_hits = [gid for gid in genome_annot_by_id
                     if "punctata" in gid.lower()]
    print(f"\n[INFO] Все записи с 'punctata' в датасете "
          f"({len(punctata_hits)} шт.) — проверь тип генома у каждой:")
    for gid in punctata_hits:
        gt = genome_annot_by_id[gid].get("genome_type", "?")
        print(f"    {gid}  [genome_type={gt}]"
              + ("  <-- ТЕТРАПЛОИД, добавь в CC_CLADE_GENOMES!" if gt == "BBCC" else ""))

    # ------------------------------------------------------------------
    # 3) cc_clade_genome_annot.tsv
    # ------------------------------------------------------------------
    cc_annot_rows = []
    for gid in found:
        row = dict(genome_annot_by_id[gid])
        gtype, species = CC_CLADE_GENOMES[gid]
        row["cc_clade_group"] = gtype
        cc_annot_rows.append(row)
    write_tsv(outdir / "cc_clade_genome_annot.tsv", cc_annot_rows,
              genome_annot_fields + ["cc_clade_group"])

    # ------------------------------------------------------------------
    # 4) cc_clade_per_genome_class.tsv (длинный формат) +
    #    добавляем feature_type (человекочитаемое имя) и % от генома
    # ------------------------------------------------------------------
    genome_size_map = {gid: float(genome_annot_by_id[gid]["genome_size"])
                        for gid in found}

    cc_class_rows = []
    for row in per_class:
        if row["genome"] not in found:
            continue
        r = dict(row)
        r["feature_type"] = class_to_type_map.get(row["classification"], "")
        gsize = genome_size_map.get(row["genome"])
        if gsize:
            r["pct_of_genome"] = f"{(float(row['total_bp']) / gsize * 100):.4f}"
        else:
            r["pct_of_genome"] = ""
        cc_class_rows.append(r)

    write_tsv(outdir / "cc_clade_per_genome_class.tsv", cc_class_rows,
              per_class_fields + ["feature_type", "pct_of_genome"])
    print(f"\n[OK] cc_clade_per_genome_class.tsv: {len(cc_class_rows)} строк")

    # ------------------------------------------------------------------
    # 5) Широкий формат: classification x genome (% от генома) +
    #    genus-wide baseline (mean_bp_per_genome / mean genome size,
    #    приблизительно, для ориентира)
    # ------------------------------------------------------------------
    all_classes = sorted({r["classification"] for r in cc_class_rows})
    pct_lookup = defaultdict(dict)  # classification -> genome -> pct
    for r in cc_class_rows:
        pct_lookup[r["classification"]][r["genome"]] = r["pct_of_genome"]

    wide_rows = []
    for cls in all_classes:
        row = {
            "classification": cls,
            "feature_type": class_to_type_map.get(cls, ""),
        }
        for gid in found:
            row[gid] = pct_lookup[cls].get(gid, "0.0000")
        # genus-wide ориентир (не строго сопоставимо с % генома конкретного
        # генома, но полезно для контекста масштаба)
        ag = across_by_class.get(cls)
        row["genus_mean_bp_per_genome"] = ag["mean_bp_per_genome"] if ag else ""
        row["genus_pct_genomes_present"] = ag["pct_genomes"] if ag else ""
        wide_rows.append(row)

    wide_fields = (["classification", "feature_type"] + found
                   + ["genus_mean_bp_per_genome", "genus_pct_genomes_present"])
    write_tsv(outdir / "cc_clade_class_wide_pct.tsv", wide_rows, wide_fields)
    print(f"[OK] cc_clade_class_wide_pct.tsv: {len(wide_rows)} классов x "
          f"{len(found)} геномов")

    # ------------------------------------------------------------------
    # 6) Текстовый отчёт
    # ------------------------------------------------------------------
    lines = []
    lines.append("=" * 70)
    lines.append("STEP 2 (v2) — CC-clade TE composition report")
    lines.append("=" * 70)
    lines.append(f"Найдено геномов CC-клады: {len(found)} / {len(CC_CLADE_GENOMES)}")
    for gid in found:
        gtype, species = CC_CLADE_GENOMES[gid]
        lines.append(f"  [{gtype:<5}] {species:<25} {gid}")
    lines.append("")
    lines.append("-- Таксономическая заметка (см. также в тексте выше) --")
    lines.append(CC_CLADE_MISSING_NOTE)
    lines.append("")
    lines.append("-- Файлы результата --")
    lines.append("  cc_clade_genome_annot.tsv      - метаданные (species, "
                  "mating_system, genome_size, intact TE summary)")
    lines.append("  cc_clade_per_genome_class.tsv  - TE-классы, длинный "
                  "формат, с % от генома")
    lines.append("  cc_clade_class_wide_pct.tsv    - широкий формат для "
                  "графиков/таблиц сравнения между видами")
    report = "\n".join(lines)
    (outdir / "cc_clade_report.txt").write_text(report, encoding="utf-8")
    print("\n" + report)
    print(f"\n[OK] Результаты в: {outdir.resolve()}")


if __name__ == "__main__":
    main()
