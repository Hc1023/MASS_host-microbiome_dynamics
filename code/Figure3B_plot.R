# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggnewscale)
    library(dplyr)
})

if (Sys.getlocale("LC_CTYPE") %in% c("C", "POSIX")) try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)

self <- grep("^--file=", commandArgs(FALSE), value = TRUE)

root <- "."

out <- repo_path(root, "Outputs")

input <- repo_path(out, "Inputs/Figure3B_refined_source_data.csv")

md5_before <- tools::md5sum(repo_input(input))

d <- read.csv(repo_input(input), check.names = FALSE)

d <- d[order(d$Display_order), ]

for (nm in c("D1_FDR_star", "D4_FDR_star", "D7_FDR_star", "Trajectory_FDR_star")) d[[nm]] <- ifelse(is.na(d[[nm]]), 
    "", as.character(d[[nm]]))

stars <- function(q) ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "")))

stopifnot(nrow(d) == 9L, identical(as.integer(table(d$Block_ID)), c(6L, 1L, 2L)))

for (tp in c("D1", "D4", "D7")) stopifnot(identical(d[[paste0(tp, "_FDR_star")]], stars(d[[paste0(tp, 
    "_p_adj_BH")]])))

stopifnot(identical(d$Trajectory_FDR_star, stars(d$overall_trajectory_FDR)))

d$y <- 9:1 - c(rep(0, 6), 0.18, 0.36, 0.36)

effects <- bind_rows(lapply(seq_along(c("D1", "D4", "D7")), function(j) {
    tp <- c("D1", "D4", "D7")[j]
    data.frame(Pathway = d$Pathway, Timepoint = tp, x = j, y = d$y, effect = d[[paste0(tp, "_estimate")]], 
        FDR = d[[paste0(tp, "_p_adj_BH")]], star = d[[paste0(tp, "_FDR_star")]])
}))

trajectory <- d %>% transmute(Pathway, x = 4.2, y, direction = factor(Trajectory_direction, levels = c("Trajectory-up", 
    "Trajectory-down")), FDR = overall_trajectory_FDR, star = Trajectory_FDR_star)

groups <- d %>% group_by(Block_ID) %>% summarise(top = max(y) + 0.43, bottom = min(y) - 0.43, center = mean(y), 
    .groups = "drop") %>% mutate(label = c("Innate sensing /\ninflammatory", "Phagolysosomal", "Ribosomal /\ntranslational"))

brackets <- bind_rows(groups %>% transmute(x = 4.81, xend = 4.81, y = bottom, yend = top), groups %>% 
    transmute(x = 4.69, xend = 4.81, y = top, yend = top), groups %>% transmute(x = 4.69, xend = 4.81, 
    y = bottom, yend = bottom))

separators <- data.frame(y = c(mean(d$y[c(6, 7)]), mean(d$y[c(7, 8)])))

header_rules <- data.frame(x = c(0.52, 3.77), xend = c(3.48, 4.63), y = 10.05)

headers <- data.frame(x = c(1, 2, 3, 4.2), y = 9.79, label = c("D1", "D4", "D7", "Direction"))

title_headers <- data.frame(x = c(2, 4.2), y = 10.64, label = c("Cross-sectional\nmortality effect", 
    "Overall\ntrajectory"))

pt_to_mm <- function(pt) pt/(72.27/25.4)

p <- ggplot() + geom_tile(data = effects, aes(x, y, fill = effect), width = 0.98, height = 0.98, color = "white", 
    linewidth = 0.18) + scale_fill_gradient2(name = "Standardized mortality effect", low = "#6487AD", 
    mid = "#FFFFFF", high = "#C77772", midpoint = 0, limits = c(-0.7, 0.7), breaks = c(-0.7, 0, 0.7), 
    labels = c("-0.7", "0", "+0.7"), space = "Lab", guide = guide_colorbar(order = 1, title.position = "top", 
        barwidth = grid::unit(1.3, "in"), barheight = grid::unit(0.075, "in"), frame.colour = NA, ticks.colour = "grey55")) + 
    ggnewscale::new_scale_fill() + geom_tile(data = trajectory, aes(x, y, fill = direction), width = 0.86, 
    height = 0.98, color = "white", linewidth = 0.18) + scale_fill_manual(name = "Trajectory direction", 
    values = c(`Trajectory-up` = "#C57E79", `Trajectory-down` = "#7896B7"), labels = c("Up", "Down"), 
    drop = FALSE, guide = guide_legend(order = 2, title.position = "top", nrow = 1, override.aes = list(color = NA))) + 
    geom_text(data = bind_rows(effects %>% select(x, y, star), trajectory %>% select(x, y, star)), aes(x, 
        y, label = star), family = "Helvetica", size = pt_to_mm(8.8), fontface = "bold", color = "#222222", 
        vjust = 0.65) + geom_segment(data = separators, aes(x = 0.5, xend = 4.63, y = y, yend = y), color = "#DDE1E5", 
    linewidth = 0.18) + geom_segment(data = brackets, aes(x, y, xend = xend, yend = yend), color = "#9AA2AC", 
    linewidth = 0.25) + geom_text(data = groups, aes(x = 5.04, y = center, label = label), hjust = 0, 
    family = "Helvetica", size = pt_to_mm(7.1), lineheight = 1.05, color = "#525B65") + geom_segment(data = header_rules, 
    aes(x, y, xend = xend, yend = y), color = "#BDC3CB", linewidth = 0.18) + geom_text(data = headers, 
    aes(x, y, label = label), size = pt_to_mm(7.1), family = "Helvetica", color = "#262B31") + geom_text(data = title_headers, 
    aes(x, y, label = label), size = pt_to_mm(7.2), family = "Helvetica", fontface = "bold", lineheight = 1.05, 
    color = "#262B31") + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = d$y, labels = d$Short_label) + 
    coord_cartesian(xlim = c(0.5, 8.15), ylim = c(0.1, 11.18), expand = FALSE, clip = "off") + labs(x = NULL, 
    y = NULL) + theme_void(base_family = "Helvetica", base_size = 7.5) + theme(axis.text.y = element_text(size = 7.6, 
    color = "#262B31", hjust = 1, margin = margin(r = 4)), legend.position = "bottom", legend.box = "horizontal", 
    legend.location = "plot", legend.justification = "left", legend.box.just = "left", legend.title = element_text(size = 6.8), 
    legend.text = element_text(size = 6.8), legend.key.width = grid::unit(0.12, "in"), legend.key.height = grid::unit(0.11, 
        "in"), legend.spacing.x = grid::unit(0.3, "in"), legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(t = 2, 
        r = 0, b = 0, l = 9), plot.margin = margin(5, 3, 3, 4), plot.background = element_rect(fill = "white", 
        color = NA))

W <- 4.2

H <- 3.6

ggsave(repo_output(repo_path(out, "Figure3B_compact.pdf"), "Figure3B_plot.R"), plot = p, width = W, 
    height = H, units = "in", device = grDevices::pdf, family = "Helvetica", encoding = "WinAnsi.enc", 
    useDingbats = FALSE, bg = "white")

ggsave(repo_output(repo_path(out, "Figure3B_compact.png"), "Figure3B_plot.R"), plot = p, width = W, 
    height = H, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")

saveRDS(p, repo_output(repo_path(out, "Figure3B_compact_ggplot.rds"), "Figure3B_plot.R"))

for (ext in c("pdf", "png")) file.copy(repo_path(out, paste0("Figure3B_compact.", ext)), repo_path(out, 
    paste0("Figure3B_compact_nofootnote.", ext)), overwrite = TRUE)

stopifnot(inherits(p, "ggplot"), nrow(effects) == 27, nrow(trajectory) == 9, identical(md5_before, tools::md5sum(repo_input(input))))

built <- ggplot_build(p)

stopifnot(nrow(built$data[[1]]) == 27, nrow(built$data[[2]]) == 9, nrow(built$data[[3]]) == 36, sum(nzchar(built$data[[3]]$label)) == 
    20)

writeLines(c("ggplot2 compact heatmap, 4 x 3 inches for PDF and PNG.", "geom_tile: three standardized effect columns and one categorical trajectory column.", 
    "ggnewscale: independent continuous/categorical fill scales with native ggplot legends.", "geom_text: FDR stars and headers; geom_segment: group brackets and separators.", 
    "All nine pathways, values, FDR stars and color direction are reused from verified source data.", 
    "Abbreviated labels do not alter the positive/negative regulation semantics of the full GO terms.", 
    "PDF device: grDevices::pdf(); no cairo_pdf or manual grid drawing.", "Editable plot object: Figure3B_compact_ggplot.rds. Previous grid version: archive_grid_compact/.", 
    "Source data and statistical legend: Figure3B_refined_source_data.csv and Figure3B_refined_legend.txt."), 
    repo_output(repo_path(out, "Figure3B_compact_notes.txt"), "Figure3B_plot.R"))

message("Created ggplot2 PDF/PNG and saved plot object: ", out)

