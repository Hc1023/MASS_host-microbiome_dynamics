# Run from the repository root. Panel sources are in manifest.csv.
source("code/shared_paths.R")
invisible(NULL)

library(tidyverse)

library(data.table)

library(ggtree)

library(data.tree)

library(ape)

library(treeio)

library(ggtreeExtra)

library(ggstar)

library(ggnewscale)

load(repo_input("Inputs/1616_microbe.rdata"))

load(repo_input("Inputs/1211_metadata.rdata"))

data_sp = read.csv(repo_input("Inputs/1527_BLsp_fg.csv"), row.names = 1)

row_max = apply(data_sp[(ncol(data_sp) - 10):ncol(data_sp)], 1, max)

row_max[order(row_max, decreasing = T)][1:10]

data_sp = data_sp[, -c((ncol(data_sp) - 10):ncol(data_sp))]

data_sp1 <- data_sp - row_max

data_sp1[data_sp1 < 0] <- 0

data_sp1 = data_sp1 %>% select(all_of(df_long$SampleID))

data_sp1 = data_sp1[rowSums(data_sp1) > 0, ]

data_sp1 = data_sp1[!grepl("Lactobacillus", rownames(data_sp1)), ]

nrow(data_sp1)

mysheets <- sisiUtils::read_excel_allsheets("Inputs/BLmapid_1101.xlsx")

mapid = mysheets[[1]]

mapid <- mapid %>% select(Pathogen_type, Genus, Species)

mapid[mapid$Genus == "Influenza A", "Genus"] = "Influenza_A"

table(rownames(data_sp1) %in% mapid$Species)

mapids = setNames(mapid$Genus, mapid$Species)

unique(mapids[rownames(data_sp1)]) %>% length()

mapid1 = mapid %>% filter(Species %in% rownames(data_sp1))

tmp = data.frame(Microbial_mass = log2(rowSums(data_sp1) + 1)) %>% rownames_to_column(var = "Species")

mapid1 = mapid1 %>% left_join(tmp, by = "Species")

head(mapid1)

mapid1$pathString <- with(mapid1, paste("Root", Pathogen_type, Genus, Species, sep = "/"))

tree <- as.Node(mapid1, pathName = "pathString")

phylo_tree <- as.phylo(tree)

p <- ggtree(phylo_tree, layout = "radial")

unique(mapid1$Pathogen_type)

cls <- list(Clade_1 = mapid1$Species[mapid1$Pathogen_type == unique(mapid1$Pathogen_type)[1]], Clade_2 = mapid1$Species[mapid1$Pathogen_type == 
    unique(mapid1$Pathogen_type)[2]], Clade_3 = mapid1$Species[mapid1$Pathogen_type == unique(mapid1$Pathogen_type)[3]])

phylo_tree <- groupOTU(phylo_tree, cls)

p = ggtree(phylo_tree, layout = "circular", aes(color = group))

p

g <- attr(phylo_tree, "group")

g[g == "0"] <- "Clade_2"

attr(phylo_tree, "group") <- droplevels(g)

p = ggtree(phylo_tree, layout = "circular", aes(color = group))

colors = RColorBrewer::brewer.pal(3, "Set2")

names(colors) = c("Clade_2", "Clade_1", "Clade_3")

p = ggtree(phylo_tree, layout = "radial", aes(color = group)) + scale_color_manual(values = colors, labels = c("Bacteria", 
    "Fungi", "Viruses"), name = "")

tmp = setNames(mapid1$Pathogen_type, mapid1$Species)

sample_dat = data.frame(label = phylo_tree$tip.label)

sample_dat$group2 = tmp[sample_dat$label]

sample_dat = sample_dat %>% left_join(mapid1 %>% select(Species, Microbial_mass), by = c(label = "Species"))

p1 = p + scale_fill_manual(values = colors, labels = c("Bacteria", "Fungi", "Viruses"), name = "") + 
    geom_fruit(geom = geom_col, data = sample_dat, mapping = aes(x = Microbial_mass, y = label, fill = group), 
        alpha = 0.6, grid.params = list(linewidth = 0.2), axis.params = list(axis = TRUE, text.size = 3, 
            nbreak = 5), pwidth = 1) + geom_tiplab(color = "black", size = 1.8, offset = 4)

pdf(file = repo_output("Outputs/03_taxtree.pdf", "FigureS2A_plot.R"), width = 8, height = 8)

print(p1)

dev.off()

