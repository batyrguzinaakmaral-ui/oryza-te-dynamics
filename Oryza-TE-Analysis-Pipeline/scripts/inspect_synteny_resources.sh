#!/bin/bash
# inspect_synteny_resources.sh
# Run from /xdisk/dcopetti/mdalam/rte
#
# Goal: find out whether synteny / structural-rearrangement / inversion data
# already exists for our 6 CC-clade genomes before building an alignment
# pipeline from scratch. Checks the gfa/, gfai/, combined_annotations/
# directories plus a broad filename search for common synteny-tool outputs
# (PAF, SyRI, MUMmer delta, anchor files, SV VCFs).

set -uo pipefail

CC_GENOMES=(
  "GCA_008326285.1_Oryza_officinalis_CC"
  "GCA_048166525.1_Oryza_minuta_BBCC"
  "GCA_048564985.1_Oryza_malampuzhaensis_BBCC"
  "GCA_047899615.1_Oryza_alta_CCDD"
  "GCA_048174585.1_Oryza_latifolia_CCDD"
  "GCA_048188845.1_Oryza_grandiglumis_CCDD"
)

echo "=== 1) Top-level look at gfa/ and gfai/ ==="
for d in gfa gfai; do
  if [ -d "$d" ]; then
    echo "--- $d/ (first 20 entries) ---"
    ls -la "$d" | head -20
    echo "total entries: $(ls "$d" | wc -l)"
    echo "file extensions present:"
    find "$d" -maxdepth 2 -type f 2>/dev/null | sed -E 's/.*\.([A-Za-z0-9]+)$/\1/' | sort | uniq -c | sort -rn | head -15
  else
    echo "$d/ NOT FOUND"
  fi
  echo
done

echo "=== 2) Do gfa/gfai contain anything for our 6 CC-clade genomes specifically? ==="
for g in "${CC_GENOMES[@]}"; do
  echo "--- $g ---"
  find gfa gfai -iname "*${g}*" 2>/dev/null | head -10
  # also try just the GCA accession in case naming truncates the rest
  acc=$(echo "$g" | cut -d'_' -f1,2)
  find gfa gfai -iname "*${acc}*" 2>/dev/null | head -10
done
echo

echo "=== 3) combined_annotations/ -- might hold gene-based synteny anchors ==="
if [ -d combined_annotations ]; then
  ls -la combined_annotations | head -20
else
  echo "combined_annotations/ NOT FOUND"
fi
echo

echo "=== 4) oryza_chromosomes.tsv -- chromosome-level info (needed for any synteny plot axis) ==="
if [ -f oryza_chromosomes.tsv ]; then
  head -3 oryza_chromosomes.tsv
  echo "..."
  echo "rows for our 6 genomes:"
  for g in "${CC_GENOMES[@]}"; do
    grep -c "$g" oryza_chromosomes.tsv 2>/dev/null | xargs echo "  $g:"
  done
else
  echo "oryza_chromosomes.tsv NOT FOUND"
fi
echo

echo "=== 5) Broad filename search across the whole rte/ tree for synteny-tool outputs ==="
echo "(.paf, .delta, SyRI output names, anchors, .sv.vcf, synteny keyword in filename)"
find . -maxdepth 4 \( \
  -iname "*.paf" -o \
  -iname "*.delta" -o \
  -iname "*syri*" -o \
  -iname "*anchors*" -o \
  -iname "*.sv.vcf*" -o \
  -iname "*synteny*" -o \
  -iname "*inversion*" -o \
  -iname "*minimap*" \
  \) 2>/dev/null | grep -v "^./edta" | head -40
echo

echo "=== 6) Do we have the raw genome FASTAs for the 6 CC-clade genomes? (needed if we must align ourselves) ==="
for g in "${CC_GENOMES[@]}"; do
  f=$(find . -maxdepth 3 -iname "${g}*.fa*" ! -iname "*.fa.mod*" 2>/dev/null | head -3)
  echo "--- $g ---"
  if [ -n "$f" ]; then
    echo "$f"
  else
    echo "  no top-level .fa found (may only exist inside edta/ as .fa.mod)"
    find edta/"$g" -maxdepth 1 -iname "*.fa*" 2>/dev/null
  fi
done
echo

echo "=== 7) Check for minimap2 / mummer / syri already installed in this environment ==="
for tool in minimap2 nucmer syri mummer; do
  which "$tool" 2>/dev/null && echo "  -> $tool FOUND" || echo "  -> $tool not in PATH"
done

echo
echo "Done. Paste the full output back."
