# plot9_tree_te.R  (v2: bars are aligned to the real vertical position of each tip)
# BUSCO subgenome tree (step8) with TE content on the leaves (step12).
# NOT tested in the Claude sandbox (no R). Needs the package 'ape'.
#
# Files (found automatically under data_root): tree.nwk (step8), leaf_te_table.tsv (step12)
# Output: figures/fig_tree_TE.png
# The tree is unrooted (no outgroup): where the root is drawn is arbitrary.

library(ape)
library(readr)
library(dplyr)

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

tree <- read.tree(find_file("tree.nwk"))
tree <- ladderize(tree)
leaf <- read_tsv(find_file("leaf_te_table.tsv"), show_col_types = FALSE)

tips <- tree$tip.label
n <- length(tips)
d <- tibble(leaf = tips) %>%
  left_join(leaf %>% select(leaf, subgenome, pct_all, intact_share), by = "leaf") %>%
  mutate(subgenome = ifelse(is.na(subgenome), sub(".*_", "", leaf), subgenome),   # punctata_B: no TE data
         intact_pct = 100 * intact_share)

cols <- c(B = "#C44E52", C = "#4C72B0", D = "#55A868")
tip_col <- unname(cols[d$subgenome])

png(file.path(out_dir, "fig_tree_TE.png"), width = 2400, height = 1300, res = 300)
layout(matrix(1:3, nrow = 1), widths = c(3.2, 1.5, 1.5))

par(mar = c(4, 1, 3, 0))
plot(tree, type = "phylogram", tip.color = tip_col, cex = 0.8, y.lim = c(0.5, n + 0.5),
     label.offset = 0.0004, main = "BUSCO subgenome tree (unrooted)")
add.scale.bar(cex = 0.7)
lp <- get("last_plot.phylo", envir = .PlotPhyloEnv)
yy <- lp$yy[seq_len(n)]          # vertical position of tip i (= tree$tip.label[i]) in the drawing

bar_panel <- function(values, xlab, xmax) {
  par(mar = c(4, 0.5, 3, 1))
  plot(NA, xlim = c(0, xmax), ylim = c(0.5, n + 0.5), yaxt = "n", xlab = xlab, ylab = "", bty = "n")
  for (i in seq_len(n)) {
    if (!is.na(values[i])) {
      rect(0, yy[i] - 0.4, values[i], yy[i] + 0.4, col = tip_col[i], border = NA)
      text(values[i], yy[i], sprintf(" %.1f", values[i]), adj = 0, cex = 0.6)
    } else {
      text(0.02 * xmax, yy[i], "n/a (no TE annotation)", cex = 0.55, col = "grey40", adj = 0)
    }
  }
}
bar_panel(d$pct_all, "Total TE, % of subgenome", 100)
bar_panel(d$intact_pct, "Intact (recent) TE, % of TE bp", max(5, max(d$intact_pct, na.rm = TRUE) * 1.25))
dev.off()
message("Saved: ", file.path(out_dir, "fig_tree_TE.png"))
