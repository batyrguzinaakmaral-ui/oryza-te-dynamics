# plot6_pct_still_intact_ranked.R
# Ranked bar chart: "% of TE still structurally intact" per feature_type,
# from plot5's feature_type-level summary (naming-convention variants
# already merged). Class I (LTR retrotransposons) vs Class II (DNA/TIR +
# Helitron) are colored separately -- the story this plot tells: DNA
# transposons retain far more recently-active copies than LTR
# retrotransposons in the CC-clade.
#
# Input: plot5_pct_still_intact_summary_by_feature_type.tsv (from plot5 script)

library(tidyverse)
try(Sys.setlocale("LC_ALL", "en_US.UTF-8"), silent = TRUE)

df <- read_tsv("plot5_pct_still_intact_summary_by_feature_type.tsv",
                col_types = cols(.default = "c")) %>%
  mutate(
    mean_full_pct = as.numeric(mean_full_pct),
    mean_intact_pct = as.numeric(mean_intact_pct),
    pct_still_intact = as.numeric(pct_still_intact),
    classification_variants = as.integer(classification_variants)
  )

cat("Rows loaded:", nrow(df), "\n")

# ---------------------------------------------------------------------------
# Classify each feature_type into a TE superclass for coloring. Anything
# not recognized falls into "Other / non-TE repeat" (LINE, SINE, rRNA,
# snRNA, pararetrovirus, repeat_fragment -- these have 0% intact by
# construction, see plot5 notes).
# ---------------------------------------------------------------------------
ltr_types <- c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon",
               "LTR_retrotransposon")
dna_types <- c("CACTA_TIR_transposon", "Mutator_TIR_transposon",
               "PIF_Harbinger_TIR_transposon", "Tc1_Mariner_TIR_transposon",
               "hAT_TIR_transposon", "helitron")

df <- df %>%
  mutate(
    superclass = case_when(
      feature_type %in% ltr_types ~ "Class I (LTR retrotransposon)",
      feature_type %in% dna_types ~ "Class II (DNA transposon)",
      TRUE ~ "Other / non-TE repeat"
    ),
    # drop rows with NA pct_still_intact (mean_full_pct == 0, if any)
  ) %>%
  filter(!is.na(pct_still_intact)) %>%
  arrange(pct_still_intact)

df$feature_type <- factor(df$feature_type, levels = df$feature_type)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
plot_rank <- ggplot(df, aes(x = pct_still_intact, y = feature_type, fill = superclass)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.2f%%  (n=%d naming variant%s)",
                                 pct_still_intact, classification_variants,
                                 if_else(classification_variants == 1, "", "s"))),
            hjust = -0.05, size = 3.2) +
  scale_fill_manual(values = c(
    "Class I (LTR retrotransposon)" = "#B2182B",
    "Class II (DNA transposon)" = "#2166AC",
    "Other / non-TE repeat" = "#999999"
  )) +
  expand_limits(x = max(df$pct_still_intact) * 1.6) +
  labs(
    title = "Share of TE Content Still Structurally Intact, by Class",
    subtitle = "6-genome average across the CC-clade; naming-convention duplicates merged (see plot5 notes)",
    x = "% of this TE class's genome content that is still structurally intact",
    y = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("plot6_pct_still_intact_ranked.png", plot_rank, width = 10, height = 7, dpi = 300)

cat("\nSaved: plot6_pct_still_intact_ranked.png\n")
cat("\nKey contrast -- Class I (LTR) vs Class II (DNA) mean % still intact:\n")
print(df %>% filter(superclass != "Other / non-TE repeat") %>%
        group_by(superclass) %>%
        summarise(mean_pct_still_intact = mean(pct_still_intact), .groups = "drop"))
