# plot2_te_burden_vs_ploidy.R
# TE burden (intact TE % of genome) vs. ploidy across the 6 CC-clade genomes.
# Input: cc_clade_genome_annot.tsv  (output of step2_cc_clade_te_composition.py)
# Addresses deliverable #4: TE expansion/contraction linked to polyploidization.

library(tidyverse)
try(Sys.setlocale("LC_ALL", "en_US.UTF-8"), silent = TRUE)

# ---------------------------------------------------------------------------
# 1) Load. Read everything as character first so a header mismatch doesn't
#    crash the script silently — we print the columns so you can eyeball them
#    before anything downstream runs.
# ---------------------------------------------------------------------------
raw <- read_tsv("cc_clade_genome_annot.tsv", col_types = cols(.default = "c"))

cat("Columns found in cc_clade_genome_annot.tsv:\n")
print(names(raw))
cat("\nRows:", nrow(raw), "\n\n")

# ---------------------------------------------------------------------------
# 2) Normalize. Derive ploidy from genome_type (CC = diploid, BBCC/CCDD =
#    allotetraploid), and a clean species label the same way as Plot 1.
# ---------------------------------------------------------------------------
df <- raw %>%
  mutate(
    n_intact_elements = as.numeric(n_intact_elements),
    total_bp_intact   = as.numeric(total_bp_intact),
    genome_size       = as.numeric(genome_size),
    pct_intact_te     = 100 * total_bp_intact / genome_size,
    species = str_extract(genome, "(?<=Oryza_)[A-Za-z_]+(?=_[A-Z]+$)"),
    ploidy = case_when(
      genome_type == "CC"              ~ "Diploid (CC)",
      genome_type %in% c("BBCC", "CCDD") ~ "Allotetraploid",
      TRUE ~ paste0("Unknown (", genome_type, ")")
    ),
    ploidy = factor(ploidy, levels = c("Diploid (CC)", "Allotetraploid",
                                        unique(ploidy[str_detect(ploidy, "^Unknown")])))
  ) %>%
  arrange(desc(pct_intact_te))

cat("Sanity check — one row per genome, should be 6:\n")
print(df %>% select(genome, species, genome_type, ploidy,
                     n_intact_elements, total_bp_intact, genome_size, pct_intact_te))

stopifnot(nrow(df) == 6, sum(is.na(df$pct_intact_te)) == 0)

# ---------------------------------------------------------------------------
# 3) PLOT A — bar chart: intact TE burden (% genome) per species, colored by
#    ploidy. Direct answer to "did TE burden change after polyploidization?"
# ---------------------------------------------------------------------------
plot_burden <- ggplot(df, aes(x = reorder(species, pct_intact_te),
                               y = pct_intact_te, fill = ploidy)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", pct_intact_te)),
            hjust = -0.15, size = 3.6) +
  coord_flip() +
  scale_fill_manual(values = c("Diploid (CC)" = "#2166AC",
                                "Allotetraploid" = "#B2182B")) +
  expand_limits(y = max(df$pct_intact_te) * 1.15) +
  labs(
    title = "Intact TE Burden across CC-clade Oryza Species",
    subtitle = "Share of genome occupied by structurally intact TE copies, by ploidy",
    x = NULL, y = "% of genome (intact TE)", fill = "Ploidy"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold"))

ggsave("plot2_te_burden_bar.png", plot_burden, width = 8, height = 5, dpi = 300)

# ---------------------------------------------------------------------------
# 4) PLOT B — scatter: genome size vs. intact TE bp, labeled by species,
#    colored by ploidy. Shows whether bigger (tetraploid) genomes carry
#    proportionally more or less intact TE, not just more in absolute bp.
# ---------------------------------------------------------------------------
plot_scatter <- ggplot(df, aes(x = genome_size / 1e6, y = total_bp_intact / 1e6,
                                color = ploidy, label = species)) +
  geom_point(size = 4) +
  ggrepel::geom_text_repel(size = 3.8, fontface = "italic", show.legend = FALSE) +
  scale_color_manual(values = c("Diploid (CC)" = "#2166AC",
                                 "Allotetraploid" = "#B2182B")) +
  labs(
    title = "Genome Size vs. Intact TE Content",
    subtitle = "Mb of structurally intact TE sequence as a function of assembly size",
    x = "Genome size (Mb)", y = "Intact TE content (Mb)", color = "Ploidy"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold"))

ggsave("plot2_genomesize_vs_te.png", plot_scatter, width = 8, height = 6, dpi = 300)

cat("\nSaved: plot2_te_burden_bar.png, plot2_genomesize_vs_te.png\n")
cat("Mean intact-TE % — diploid:",
    round(mean(df$pct_intact_te[df$ploidy == "Diploid (CC)"]), 3),
    "| allotetraploid:",
    round(mean(df$pct_intact_te[df$ploidy == "Allotetraploid"]), 3), "\n")
