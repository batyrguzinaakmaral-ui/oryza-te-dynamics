#!/bin/bash
# inspect_busco_for_synteny.sh
# Run from /xdisk/dcopetti/mdalam/rte
#
# Goal: check whether a full per-gene BUSCO/compleasm table (not just the
# summary % already in oryza_genomes.tsv) exists for our 6 CC-clade genomes.
# A full_table.tsv gives single-copy ortholog coordinates per genome, which
# is enough to build a macro-synteny dot plot and spot large inversions --
# much cheaper than a full genome-vs-genome alignment with minimap2/SyRI.

set -uo pipefail

CC_GENOMES=(
  "GCA_008326285.1_Oryza_officinalis_CC"
  "GCA_048166525.1_Oryza_minuta_BBCC"
  "GCA_048564985.1_Oryza_malampuzhaensis_BBCC"
  "GCA_047899615.1_Oryza_alta_CCDD"
  "GCA_048174585.1_Oryza_latifolia_CCDD"
  "GCA_048188845.1_Oryza_grandiglumis_CCDD"
)

echo "=== 1) What's in cm/ ? ==="
if [ -d cm ]; then
  ls -la cm | head -20
  echo "total entries: $(ls cm | wc -l)"
  echo "--- is it flat files or per-genome subfolders? ---"
  find cm -maxdepth 1 -type d | head -10
else
  echo "cm/ NOT FOUND"
fi
echo

echo "=== 2) Look for full_table (BUSCO/compleasm per-gene output) anywhere in cm/ for our 6 genomes ==="
for g in "${CC_GENOMES[@]}"; do
  echo "--- $g ---"
  find cm -iname "*${g}*" 2>/dev/null | head -10
  acc=$(echo "$g" | cut -d'_' -f1,2)
  find cm -iname "*${acc}*full_table*" 2>/dev/null | head -10
  find cm -iname "*${acc}*" -iname "*.tsv" 2>/dev/null | head -10
done
echo

echo "=== 3) Broad search for 'full_table' anywhere under rte/ (compleasm/BUSCO naming) ==="
find . -maxdepth 6 -iname "*full_table*" 2>/dev/null | head -30
echo

echo "=== 4) If found, show header + a few rows of one full_table for sanity ==="
sample=$(find cm -iname "*full_table*" 2>/dev/null | head -1)
if [ -n "$sample" ]; then
  echo "Sample file: $sample"
  head -10 "$sample"
else
  echo "No full_table file found to preview."
fi
echo

echo "=== 5) Conda availability check (for installing minimap2/SyRI as fallback plan) ==="
which conda 2>/dev/null && echo "conda FOUND" || echo "conda not in PATH"
conda env list 2>/dev/null || echo "conda env list failed"
ls _setup/.conda_pkgs/ 2>/dev/null | grep -i -E "minimap|syri|mummer" || echo "no relevant cached packages found"

echo
echo "Done. Paste the full output back."
