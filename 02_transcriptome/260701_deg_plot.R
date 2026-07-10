rm(list = ls())

library(tidyverse)
library(ggpubr)

# This script starts from the three already generated DEG tables.
# Positive logFC = higher in 28-day non-survivors; negative logFC = lower.
script_dir <- "/Users/huangsisi/workspace/MASS/sepsis_microbiome/MASS_mortality-main"
in_dir  <- file.path(script_dir, "Outputs/260701_deg_D1")
out_dir <- file.path(script_dir, "Outputs/260701_deg_D1/immune_remodeling_plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

timepoints <- c("D1", "D4", "D7")

# ---------- 1. Read and combine DEG results ----------
read_deg <- function(tp) {
  f <- file.path(in_dir, paste0("deg_", tp, ".csv"))
  if (!file.exists(f)) stop("Cannot find: ", f)

  x <- read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (any(names(x) == "")) names(x)[names(x) == ""] <- "ensembl_id"
  required <- c("gene_symbol", "logFC", "P.Value", "adj.P.Val")
  missing_cols <- setdiff(required, names(x))
  if (length(missing_cols) > 0) {
    stop(basename(f), " is missing: ", paste(missing_cols, collapse = ", "))
  }

  x %>%
    filter(!is.na(gene_symbol), gene_symbol != "") %>%
    # If several Ensembl IDs map to one symbol, retain the most significant row.
    arrange(adj.P.Val, P.Value) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>%
    transmute(
      Timepoint = tp,
      gene = gene_symbol,
      logFC = as.numeric(logFC),
      p_value = as.numeric(P.Value),
      FDR = as.numeric(adj.P.Val)
    )
}

deg_all <- map_dfr(timepoints, read_deg)

# Focused display: adaptive immunity, NF-kB feedback, and IFI27 only.
axis_genes <- tribble(
  ~axis,                         ~gene,      ~gene_order,
  "Adaptive immunity",          "CD4",       1,
  "Adaptive immunity",          "CD3D",      2,
  "Adaptive immunity",          "CD3E",      3,
  "Adaptive immunity",          "TRAC",      4,
  "Adaptive immunity",          "IL7R",      5,
  "NF-kB feedback",             "TNFAIP3",   1,
  "NF-kB feedback",             "NFKBIA",    2,
  "IFI27",                      "IFI27",     1
) %>%
  mutate(axis = factor(axis, levels = c(
    "Adaptive immunity", "NF-kB feedback", "IFI27"
  )))

plot_df <- axis_genes %>%
  left_join(deg_all, by = "gene") %>%
  mutate(
    Timepoint = factor(Timepoint, levels = timepoints),
    evidence = case_when(
      FDR < 0.05     ~ "FDR < 0.05",
      p_value < 0.05 ~ "P < 0.05",
      TRUE           ~ "P >= 0.05"
    ),
    evidence = factor(evidence, levels = c(
      "FDR < 0.05", "P < 0.05", "P >= 0.05"
    )),
    neglog10_FDR = pmin(-log10(pmax(FDR, 1e-300)), 6),
    gene_facet = factor(gene, levels = rev(unique(axis_genes$gene)))
  ) %>%
  mutate(
    tp_num = as.numeric(Timepoint),
    gene_y = unname(c(CD4 = 8, CD3D = 7, CD3E = 6, TRAC = 5, IL7R = 4,
                      TNFAIP3 = 2.6, NFKBIA = 1.6, IFI27 = 0)[gene])
  )

write.csv(
  plot_df %>% dplyr::select(axis, gene, Timepoint, logFC, p_value, FDR, evidence),
  file.path(out_dir, "immune_axis_gene_plotting_data.csv"),
  row.names = FALSE
)

# Shared journal-style design for two same-height, side-by-side panels.
blue <- "#3F6FA6"
red  <- "#C75B4D"
gene_cols <- c(
  "CD4" = "#2F5597", "CD3D" = "#5577A8", "CD3E" = "#7691B4",
  "TRAC" = "#93A8C5", "IL7R" = "#AAB9D0",
  "TNFAIP3" = "#B94A48", "NFKBIA" = "#DE796D", "IFI27" = "#D49A2A"
)

label_df <- plot_df %>%
  filter(Timepoint == "D7") %>%
  mutate(label_y = case_when(
    gene == "CD4"    ~ -1.02,
    gene == "CD3D"   ~ -0.28,
    gene == "CD3E"   ~ -0.40,
    gene == "TRAC"   ~ -0.52,
    gene == "IL7R"   ~ -0.64,
    gene == "TNFAIP3" ~ 0.55,
    gene == "NFKBIA"  ~ 0.44,
    TRUE ~ logFC
  ))

theme_paper <- theme_classic(base_size = 9.5, base_family = "sans") +
  theme(
    axis.line = element_line(color = "#5F6870", linewidth = 0.35),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 11, color = "#1E2830"),
    plot.subtitle = element_text(color = "#5E6870", size = 8.5),
    plot.caption = element_text(color = "#6D747A", size = 7.5, hjust = 0),
    legend.title = element_text(face = "bold", size = 8.5),
    legend.text = element_text(size = 8),
    plot.margin = margin(8, 8, 6, 8)
  )

# ---------- 2. Main figure: minimalist effect-size matrix ----------
p_heatmap <- ggplot(plot_df, aes(tp_num, gene_y)) +
  geom_vline(xintercept = 1:3, color = "#E5E8EA", linewidth = 0.45) +
  annotate("segment", x = 0.72, xend = 0.72, y = 3.7, yend = 8.3,
           color = blue, linewidth = 2.2, lineend = "round") +
  annotate("segment", x = 0.72, xend = 0.72, y = 1.35, yend = 2.85,
           color = red, linewidth = 2.2, lineend = "round") +
  annotate("segment", x = 0.72, xend = 0.72, y = -0.25, yend = 0.25,
           color = "#D49A2A", linewidth = 2.2, lineend = "round") +
  geom_point(
    data = filter(plot_df, evidence == "FDR < 0.05"),
    aes(fill = logFC),
    shape = 24, size = 4.6, color = "#30363B", stroke = 0.45
  ) +
  geom_point(
    data = filter(plot_df, evidence == "P < 0.05"),
    aes(fill = logFC),
    shape = 21, size = 4.3, color = "#30363B", stroke = 0.45
  ) +
  geom_point(
    data = filter(plot_df, evidence == "P >= 0.05"),
    aes(color = logFC),
    shape = 21, size = 4.3, fill = "white", stroke = 0.8
  ) +
  scale_fill_gradient2(
    low = "#3B75AF", mid = "white", high = "#C7473A",
    midpoint = 0, limits = c(-1.25, 1.25), oob = scales::squish,
    name = "log2FC\n(Mortality vs Survival)"
  ) +
  scale_color_gradient2(
    low = "#3B75AF", mid = "grey55", high = "#C7473A",
    midpoint = 0, limits = c(-1.25, 1.25), oob = scales::squish,
    guide = "none"
  ) +
  labs(
    x = NULL, y = NULL,
    caption = "Triangle: FDR<0.05; filled circle: P<0.05 but FDR>=0.05; open circle: P>=0.05."
  ) +
  scale_x_continuous(breaks = 1:3, labels = timepoints, limits = c(0.62, 3.12)) +
  scale_y_continuous(
    breaks = c(8, 7, 6, 5, 4, 2.6, 1.6, 0),
    labels = c("CD4", "CD3D", "CD3E", "TRAC", "IL7R", "TNFAIP3", "NFKBIA", "IFI27"),
    limits = c(-0.45, 8.45), expand = expansion(mult = 0)
  ) +
  theme_paper +
  theme(
    axis.text.y = element_text(face = "italic"),
    axis.text.x = element_text(face = "bold", color = "#3F464C"),
    legend.position = "bottom",
    legend.key.width = unit(25, "pt"),
    axis.line.y = element_blank()
  )

# ---------- 3. Companion figure: trajectories of representative genes ----------
p_trajectory <- ggplot(
  plot_df,
  aes(tp_num, logFC, group = gene, color = gene)
) +
  geom_hline(yintercept = 0, color = "#7C858C", 
             linetype = "dashed", linewidth = 0.35) +
  geom_line(linewidth = 0.85, alpha = 0.9) +
  geom_point(aes(shape = evidence), size = 2.7, stroke = 0.75) +
  geom_segment(
    data = label_df,
    aes(x = 3.02, xend = 3.11, y = logFC, yend = label_y, color = gene),
    inherit.aes = FALSE, linewidth = 0.45
  ) +
  geom_text(
    data = label_df,
    aes(x = 3.14, y = label_y, label = gene, color = gene),
    inherit.aes = FALSE, hjust = 0, size = 2.65, fontface = "italic"
  ) +
  facet_grid(axis ~ ., scales = "free_y", switch = "y") +
  scale_shape_manual(
    values = c("FDR < 0.05" = 17, "P < 0.05" = 16, "P >= 0.05" = 1),
    drop = FALSE
  ) +
  scale_color_manual(values = gene_cols) +
  labs(
    x = NULL, y = "log2FC (Mortality vs Survival)",
    shape = "Statistical evidence"
  ) +
  theme_paper +
  scale_x_continuous(breaks = 1:3, labels = timepoints,
                     limits = c(0.9, 3.6), expand = expansion(mult = c(0, 0))) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    axis.text.x = element_text(face = "bold", color = "#3F464C"),
    panel.spacing.y = unit(0.12, "lines"),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 8.5,
                                     color = "#4E5860", margin = margin(r = 4)),
    axis.line.x = element_line(color = "#69737B"),
    axis.line.y = element_blank()
  ) +
  guides(
    color = "none",
    shape = guide_legend(nrow = 1)
  )

# ---------- 4. Save publication-ready outputs ----------
ggsave(file.path(out_dir, "Fig_immune_axis_dot_heatmap.pdf"),
       p_heatmap, width = 3, height = 4)

ggsave(file.path(out_dir, "Fig_immune_gene_trajectories.pdf"),
       p_trajectory, width = 3.2, height = 4)

ggsave(file.path(out_dir, "Fig_immune_axis_dot_heatmap_lgd.pdf"),
       p_heatmap, width = 5, height = 4)

ggsave(file.path(out_dir, "Fig_immune_gene_trajectories_lgd.pdf"),
       p_trajectory, width = 5, height = 4)

message("Done. Figures and plotting data are in: ", out_dir)
