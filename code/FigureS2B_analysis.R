# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(data.table)

library(Hmisc)

load(repo_input("Inputs/1616_microbe.rdata"))

{
    genus_sum2 = data
    group_vec <- mapid_vec[rownames(genus_sum2)]
    prev <- rowSums(genus_sum2 > 0)
    bac_df <- genus_sum2[group_vec == "Bacteria", , drop = FALSE]
    fun_df <- genus_sum2[group_vec == "Fungi", , drop = FALSE]
    virus_df <- genus_sum2[group_vec == "Viruses", , drop = FALSE]
    virus_ordered <- virus_df[order(prev[rownames(virus_df)], decreasing = TRUE), ]
    bac_ordered <- bac_df[order(prev[rownames(bac_df)], decreasing = TRUE), ]
    virus_top10 <- head(virus_ordered, 10)
    bac_top10 <- head(bac_ordered, 10)
    genus_sum3 <- rbind(bac_top10, fun_df, virus_top10)
}

rownames(genus_sum3)[rownames(genus_sum3) == "HHV-4"] = "EBV"

overallres = rcorr(t(genus_sum3), type = "spearman")

diag(overallres$P) <- 0

pdf(repo_output("Outputs/03_microbe_cor.pdf", "FigureS2B_analysis.R"), width = 6, height = 4.5)

corrplot::corrplot(overallres$r, p.mat = overallres$P, sig.level = 0.05, insig = "blank", tl.cex = 0.8, 
    type = "lower", method = "color", outline = TRUE, col = rev(RColorBrewer::brewer.pal(n = 10, name = "PuOr")), 
    tl.srt = 45, tl.offset = 0.3, tl.col = "black", mar = c(0.1, 0.1, 0.1, 0.1))

{
    species <- rownames(overallres$r)
    n <- length(species)
    row_num <- which(species == "Aspergillus")
    x_left <- 0.5
    x_right <- row_num + 0.5
    y_pos_upper <- n - row_num + 1 + 0.5
    y_pos_lower <- n - row_num + 1 - 0.5
    rect(x_left, y_pos_lower, x_right, y_pos_upper, border = "black", lwd = 0.3)
}

dev.off()

