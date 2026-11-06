#!/bin/bash
# inspect_cm_and_nlr.sh  -- run from /xdisk/dcopetti/mdalam/rte
# 1) Do the 6 CC-clade genomes have poales_odb12 per-gene tables in cm/? What columns?
# 2) What does nlr/ already contain (README/METHODS) and does it cover our genomes?

CC_GENOMES=(
  "GCA_008326285.1_Oryza_officinalis_CC"
  "GCA_048166525.1_Oryza_minuta_BBCC"
  "GCA_048564985.1_Oryza_malampuzhaensis_BBCC"
  "GCA_047899615.1_Oryza_alta_CCDD"
  "GCA_048174585.1_Oryza_latifolia_CCDD"
  "GCA_048188845.1_Oryza_grandiglumis_CCDD"
)

echo "=== 1) cm/ files for our 6 genomes ==="
for g in "${CC_GENOMES[@]}"; do
  echo "--- $g ---"
  ls -la cm/ | grep "$g" | grep poales_odb12
done

echo
echo "=== 2) Header + first rows of one table (officinalis) ==="
f="cm/GCA_008326285.1_Oryza_officinalis_CC.poales_odb12.tsv.gz"
if [ -f "$f" ]; then
  zcat "$f" | head -6
  echo "total rows: $(zcat "$f" | wc -l)"
  echo "--- summary.txt ---"
  cat "${f%.tsv.gz}.summary.txt" 2>/dev/null
else
  echo "NOT FOUND: $f"
fi

echo
echo "=== 3) nlr/README.md (first 60 lines) ==="
head -60 nlr/README.md

echo
echo "=== 4) nlr/METHODS.md (first 50 lines) ==="
head -50 nlr/METHODS.md

echo
echo "=== 5) nlr/miniprot_all198_summary.tsv: header + rows for our 6 genomes ==="
head -2 nlr/miniprot_all198_summary.tsv
for g in "${CC_GENOMES[@]}"; do grep "$g" nlr/miniprot_all198_summary.tsv; done

echo
echo "=== 6) nlr/ranking/ and nlr/tracker/ : entries for our genomes ==="
for g in "${CC_GENOMES[@]}"; do
  echo "--- $g ---"
  ls nlr/ranking nlr/tracker 2>/dev/null | grep -i "$g" | head -5
done

echo
echo "Done. Paste the full output back."
