# plot5_intact_vs_full_genome_gap.R
# Compares step2 (intact-only, cc_clade_class_wide_pct.tsv) against
# step3 (full-genome, full_genome_class_wide_pct.tsv) per TE class/species.
# The gap (full % - intact %) is the "decayed / non-structurally-intact"
# fraction of each TE class -- i.e. older, degraded copies that EDTA's
# homology search still finds but that no longer pass the intact-structure
# filter. Direct evidence for deliverable #4 (TE activity/age dynamics).
#
# Inputs (both from your data folder):
#   cc_clade_class_wide_pct.tsv       (step2, intact-only; columns are full
#                                       accession IDs like GCA_..._CC)
#   full_genome_class_wide_pct.tsv    (step3, full genome; columns are plain
#                                       species names like Oryza_officinalis)
# Note the two files use DIFFERENT column-naming conventions for species --
# this script normalizes both to a bare species name before joining.

library(tidyverse)
try(Sys.setlocale("LC_ALL", "en_US.UTF-8"), silent = TRUE)

# ---------------------------------------------------------------------------
# 1) Load + reshape step2 (intact-only). Columns are accession IDs; strip to
#    species name the same way as earlier scripts.
# ---------------------------------------------------------------------------
intact_wide <- read_tsv("cc_clade_class_wide_pct.tsv", col_types = cols(.default = "c"))

cat("Columns in cc_clade_class_wide_pct.tsv:\n")
print(names(intact_wide))

intact_long <- intact_wide %>%
  select(-any_of(c("genus_mean_bp_per_genome", "genus_pct_genomes_present"))) %>%
  pivot_longer(cols = -c(classification, feature_type),
               names_to = "genome_col", values_to = "intact_pct") %>%
  mutate(
    species = str_extract(genome_col, "(?<=Oryza_)[A-Za-z_]+(?=_[A-Z]+$)"),
    intact_pct = as.numeric(intact_pct)
  ) %>%
  select(classification, feature_type, species, intact_pct)

# ---------------------------------------------------------------------------
# 2) Load + reshape step3 (full genome). Columns are already bare species
#    names (Oryza_officinalis, etc.) -- strip the "Oryza_" prefix to match.
# ---------------------------------------------------------------------------
full_wide <- read_tsv("full_genome_class_wide_pct.tsv", col_types = cols(.default = "c"))

cat("\nColumns in full_genome_class_wide_pct.tsv:\n")
print(names(full_wide))

full_long <- full_wide %>%
  pivot_longer(cols = -c(classification, feature_type),
               names_to = "species_col", values_to = "full_pct") %>%
  mutate(
    species = str_remove(species_col, "^Oryza_"),
    full_pct = as.numeric(full_pct)
  ) %>%
  select(classification, feature_type, species, full_pct)

# ---------------------------------------------------------------------------
# 3) Join on (classification, species). full_join + fill 0, since a class
#    might be entirely absent from one side (e.g. a class EDTA's intact
#    filter never calls for a given genome).
# ---------------------------------------------------------------------------
df <- full_join(intact_long, full_long, by = c("classification", "species")) %>%
  mutate(
    feature_type = coalesce(feature_type.x, feature_type.y),
    intact_pct = replace_na(intact_pct, 0),
    full_pct = replace_na(full_pct, 0),
    decayed_pct = pmax(full_pct - intact_pct, 0),   # guard against rounding noise
    gap_flag = if_else(full_pct < intact_pct, "intact > full (check!)", "ok")
  ) %>%
  select(classification, feature_type, species, intact_pct, full_pct, decayed_pct, gap_flag)

cat("\nRows where full % < intact % (should be none / rounding only):\n")
print(df %>% filter(gap_flag != "ok"))

stopifnot(nrow(df) > 0, sum(is.na(df$species)) == 0)

cat("\nFull merged table (sample):\n")
print(df %>% arrange(classification, species) %>% head(20))

# ---------------------------------------------------------------------------
# 4) Stacked bar: Intact (bottom) + Decayed (top) = Full %, per species,
#    faceted by classification. Free y-scale per facet since magnitudes
#    range from ~0.01% (MITE) to ~20%+ (Gypsy/Helitron).
# ---------------------------------------------------------------------------
plot_df <- df %>%
  select(classification, feature_type, species, Intact = intact_pct, Decayed = decayed_pct) %>%
  pivot_longer(cols = c(Intact, Decayed), names_to = "component", values_to = "pct") %>%
  mutate(component = factor(component, levels = c("Decayed", "Intact")))

plot_gap <- ggplot(plot_df, aes(x = species, y = pct, fill = component)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Intact" = "#2166AC", "Decayed" = "#D1D1D1")) +
  facet_wrap(~ feature_type, scales = "free_y", ncol = 4) +
  labs(
    title = "Intact vs. Decayed TE Fraction by Class and Species",
    subtitle = "Full-genome % TE (step3) split into structurally-intact (step2) and decayed/degraded remainder",
    x = NULL, y = "% of genome", fill = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 7),
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("plot5_intact_vs_decayed_gap.png", plot_gap, width = 13, height = 10, dpi = 300)

# ---------------------------------------------------------------------------
# 5) Summary table: which classes are MOST decayed (lowest intact/full
#    ratio) vs. MOST recently active (highest ratio) -- genus-wide average
#    across the 6 CC-clade genomes.
# ---------------------------------------------------------------------------
summary_tbl <- df %>%
  group_by(classification, feature_type) %>%
  summarise(
    mean_full_pct = mean(full_pct),
    mean_intact_pct = mean(intact_pct),
    pct_still_intact = if_else(mean_full_pct > 0,
                                100 * mean_intact_pct / mean_full_pct, NA_real_),
    .groups = "drop"
  ) %>%
  arrange(pct_still_intact)

cat("\n'% of TE still structurally intact' per class (6-genome average, sorted",
    "lowest = most decayed first):\n")
print(summary_tbl)

write_tsv(summary_tbl, "plot5_pct_still_intact_summary_by_classification.tsv")

# ---------------------------------------------------------------------------
# 6) SAME summary, but aggregated by feature_type (not raw classification).
#    Needed because the full-genome GFF3 mixes two naming conventions for
#    the same biological superfamily, e.g. "DNA/DTC" (EDTA-native, can be
#    structurally intact) and "DNA/CMC-EnSpm" (Repbase-style synonym for
#    CACTA, homology-only, always 0% intact by construction). Grouping by
#    raw classification fragments one real category into misleading rows
#    (one at 0% intact). feature_type is shared across both naming
#    conventions, so THIS table is the biologically meaningful one -- use
#    it in the thesis, not the raw-classification version above.
# ---------------------------------------------------------------------------
summary_by_type <- df %>%
  group_by(feature_type) %>%
  summarise(
    mean_full_pct = sum(full_pct) / n_distinct(species),
    mean_intact_pct = sum(intact_pct) / n_distinct(species),
    pct_still_intact = if_else(mean_full_pct > 0,
                                100 * mean_intact_pct / mean_full_pct, NA_real_),
    classification_variants = n_distinct(classification),
    .groups = "drop"
  ) %>%
  arrange(pct_still_intact)

cat("\n'% of TE still structurally intact' per FEATURE TYPE (merges naming",
    "variants of the same superfamily -- use this table, not the per-",
    "classification one, for the thesis):\n")
print(summary_by_type)

write_tsv(summary_by_type, "plot5_pct_still_intact_summary_by_feature_type.tsv")
cat("\nSaved: plot5_intact_vs_decayed_gap.png,",
    "plot5_pct_still_intact_summary_by_classification.tsv,",
    "plot5_pct_still_intact_summary_by_feature_type.tsv\n")
