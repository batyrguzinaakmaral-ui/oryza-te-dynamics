# plot8_karyotype_and_totals.R
# Figures from step11 (+ step9 / step10 for subgenomes and exchange points).
# NOT tested in the Claude sandbox (no R there): if a line fails, send the error text.
#
# 1) fig_total_TE_per_genome.png    one bar per genome: total % TE (not split by class)
# 2) fig_karyotype_<species>_density.png   TE density in 1 Mb windows along every chromosome
#                                          + subgenome strip (step9) + exchange points (step10)
# 3) fig_karyotype_<species>_recent.png    share of TE that is structurally intact (recent) per window
#
# Set data_root below to the folder that contains the downloaded step folders (any layout works):
#   te_atlas_step11/te_density_windows.tsv, te_atlas_step11/genome_te_totals.tsv
#   te_atlas_step9/subgenome_segments.tsv
#   te_atlas_step10/breakpoints_subgenome.tsv

library(tidyverse)

data_root <- "D:/Windows/Загрузки"          # <- folder that contains te_atlas_step8, _step9, _step10, _step11, _step12
out_dir   <- file.path(data_root, "figures")
dir.create(out_dir, showWarnings = FALSE)

# finds a file by name anywhere under data_root (works whether the files are in step folders or in one folder)
find_file <- function(name) {
  hits <- list.files(data_root, pattern = paste0("^", gsub("\\.", "\\\\.", name), "$"),
                     recursive = TRUE, full.names = TRUE)
  hits <- hits[!grepl("figures", hits)]
  if (length(hits) == 0) stop("File not found under ", data_root, ": ", name)
  if (length(hits) > 1) message("Several matches for ", name, ", using: ", hits[1])
  hits[1]
}

win  <- read_tsv(find_file("te_density_windows.tsv"), show_col_types = FALSE)
tot  <- read_tsv(find_file("genome_te_totals.tsv"), show_col_types = FALSE)
seg  <- read_tsv(find_file("subgenome_segments.tsv"), show_col_types = FALSE)
junc <- read_tsv(find_file("breakpoints_subgenome.tsv"), show_col_types = FALSE)

# ---- 1) total %TE per genome ---------------------------------------------------------
p1 <- tot %>%
  mutate(species = fct_reorder(species, pct_all)) %>%
  ggplot(aes(species, pct_all)) +
  geom_col(fill = "#4C72B0", width = 0.7) +
  geom_text(aes(label = sprintf("%.1f%%", pct_all)), hjust = -0.1, size = 3.5) +
  coord_flip(ylim = c(0, max(tot$pct_all) * 1.12)) +
  labs(x = NULL, y = "Transposable elements, % of assembly",
       title = "Total TE content per genome (EDTA, all classes)") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "fig_total_TE_per_genome.png"), p1, width = 7, height = 3.5, dpi = 300)

# ---- 2/3) karyotypes -------------------------------------------------------------------
sub_cols <- c(B = "#C44E52", C = "#4C72B0", D = "#55A868")

karyo <- function(sp, metric = c("density", "recent")) {
  metric <- match.arg(metric)
  w <- win %>% filter(species == sp)
  if (nrow(w) == 0) return(invisible(NULL))
  # order chromosomes by length, label them 1..n
  ord <- w %>% group_by(chrom) %>% summarise(len = max(win_end)) %>% arrange(desc(len))
  w <- w %>% mutate(chrom = factor(chrom, levels = ord$chrom),
                    mb0 = win_start / 1e6, mb1 = win_end / 1e6,
                    value = if (metric == "density") frac_all
                            else ifelse(frac_all > 0, frac_intact / frac_all, NA_real_))
  s <- seg %>% filter(species == sp, seq %in% ord$chrom) %>%
    mutate(chrom = factor(seq, levels = ord$chrom))
  j <- junc %>% filter(query == sp, q_chr %in% ord$chrom) %>%
    mutate(chrom = factor(q_chr, levels = ord$chrom), mb = (bp_start + bp_end) / 2e6)
  lab <- if (metric == "density") "TE density" else "Intact share of TE\n(recent copies)"
  p <- ggplot(w) +
    geom_rect(aes(xmin = mb0, xmax = mb1, ymin = 0, ymax = 1, fill = value)) +
    scale_fill_viridis_c(option = if (metric == "density") "magma" else "viridis",
                         name = lab, na.value = "grey90", limits = c(0, 1)) +
    facet_grid(chrom ~ ., switch = "y") +
    labs(x = "Position, Mb", y = NULL,
         title = paste0(sp, ": ", if (metric == "density") "TE density (1 Mb windows)"
                        else "share of recent (intact) TE")) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), panel.grid = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 6), panel.spacing = unit(0.1, "lines"))
  if (nrow(s) > 0) {
    p <- p + geom_segment(data = s, aes(x = start / 1e6, xend = end / 1e6, y = -0.3, yend = -0.3,
                                        colour = label), linewidth = 2.5, inherit.aes = FALSE) +
      scale_colour_manual(values = sub_cols, name = "Subgenome")
  }
  if (nrow(j) > 0) {
    p <- p + geom_vline(data = j, aes(xintercept = mb), colour = "black", linetype = "dashed", linewidth = 0.4)
  }
  ggsave(file.path(out_dir, paste0("fig_karyotype_", sp, "_", metric, ".png")),
         p, width = 9, height = 1.1 + 0.55 * nrow(ord), dpi = 300, limitsize = FALSE)
}

for (sp in unique(win$species)) {
  karyo(sp, "density")
  karyo(sp, "recent")
}
message("Done. Figures in: ", out_dir)
