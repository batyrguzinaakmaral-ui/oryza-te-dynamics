# plot4_te_length_distributions.R
# Element-length summary (min/median/mean/max) per TE classification, across
# the 6 CC-clade genomes. A rough proxy for "how decayed/fragmented" each
# TE class's intact copies are per species -- useful for the recent-activity
# / evolutionary-dynamics angle (deliverable #4).
#
# Input: cc_clade_per_genome_class.tsv (long format, per genome x class)
# Columns expected: genome, classification, feature_type, n_elements,
#                    total_bp, mean_len, median_len, min_len, max_len
#
# IMPORTANT CAVEAT (read before interpreting the plot):
# This file only has SUMMARY stats (min/median/mean/max) per genome x class,
# not the raw per-element length list. So this is NOT a true boxplot/violin
# (no quartiles, no distribution shape) -- it's a min-max range with the
# median marked. Don't call it a "distribution" in the thesis without this
# caveat; call it a length-range summary.

library(tidyverse)
try(Sys.setlocale("LC_ALL", "en_US.UTF-8"), silent = TRUE)

# ---------------------------------------------------------------------------
# 1) Load, keep only rows where the class actually has elements in that
#    genome (n_elements > 0) -- a 0-element class has no length to show.
# ---------------------------------------------------------------------------
raw <- read_tsv("cc_clade_per_genome_class.tsv", col_types = cols(.default = "c"))

cat("Columns in cc_clade_per_genome_class.tsv:\n")
print(names(raw))

df <- raw %>%
  mutate(
    n_elements = as.numeric(n_elements),
    mean_len   = as.numeric(mean_len),
    median_len = as.numeric(median_len),
    min_len    = as.numeric(min_len),
    max_len    = as.numeric(max_len),
    species = str_extract(genome, "(?<=Oryza_)[A-Za-z_]+(?=_[A-Z]+$)"),
    genome_type_tag = str_extract(genome, "[A-Z]+$"),
    ploidy = if_else(genome_type_tag == "CC", "Diploid (CC)", "Allotetraploid")
  ) %>%
  filter(n_elements > 0, !is.na(min_len), !is.na(max_len))

cat("\nRows with length data (n_elements > 0):", nrow(df), "of", nrow(raw), "total\n")
cat("Species found:", paste(unique(df$species), collapse = ", "), "\n")
cat("Classes found:", paste(unique(df$feature_type), collapse = ", "), "\n\n")

stopifnot(nrow(df) > 0, sum(is.na(df$median_len)) == 0)

# ---------------------------------------------------------------------------
# 2) Range plot: one row per (species x classification), segment = min-max,
#    point = median. Faceted by TE classification so all 14 classes are
#    readable at once; x-axis log10 (element lengths span orders of
#    magnitude, e.g. MITEs ~100 bp vs. LTR retrotransposons ~5-10 kb).
# ---------------------------------------------------------------------------
plot_len <- ggplot(df, aes(y = species, color = ploidy)) +
  geom_segment(aes(x = min_len, xend = max_len, yend = species), linewidth = 1.1, alpha = 0.7) +
  geom_point(aes(x = median_len), size = 2.2) +
  scale_x_log10(labels = scales::label_number(big.mark = ",")) +
  scale_color_manual(values = c("Diploid (CC)" = "#2166AC", "Allotetraploid" = "#B2182B")) +
  facet_wrap(~ feature_type, scales = "free_x", ncol = 4) +
  labs(
    title = "Intact TE Element Length Range by Class and Species",
    subtitle = "Min-max range per genome x class, point = median length (log10 bp scale)",
    x = "Element length, bp (log10 scale)", y = NULL, color = "Ploidy"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y = element_text(face = "italic", size = 8),
    strip.text = element_text(face = "bold", size = 8),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("plot4_te_length_ranges.png", plot_len, width = 13, height = 9, dpi = 300)

cat("Saved: plot4_te_length_ranges.png\n")
cat("\nReminder: this is a min/median/max RANGE summary, not a true\n")
cat("distribution (no raw per-element lengths available in this file).\n")
