rm(list = ls())
library(tidyverse)
library(magrittr)
library(data.table)
library(ConsensusTranscriptomicSubtype)
run_subtype_classifier
library(edgeR)
library(rstatix)
library(ggpubr)
library(gridExtra)

load('Inputs/1211_microbe.rdata')
load('Inputs/1211_transcriptome.rdata')
load('Inputs/1211_metadata.rdata')
# D1/4/7 have samples or death
meta_analysis = meta %>% filter(!is.na(D1) & !is.na(D4.death) & !is.na(D7.death))
table(meta_analysis$Mortality28d)
# 0   1 
# 110 107

meta_analysis %<>%
  dplyr::select(
    HumanID,
    Mortality28d,
    D1,
    D4 = D4.death,
    D7 = D7.death
  ) %>%
  pivot_longer(
    cols = c(D1, D4, D7),
    names_to = "Timepoint",
    values_to = "SampleID"
  )

df1 = meta_analysis
data("exp_core_g", package = "ConsensusTranscriptomicSubtype", 
     envir = environment())
g <- data.frame(gene = rownames(exp_core_g))

# new_expr_data must be a data.frame or matrix
# First column = gene identifiers (optional), rownames = Ensembl IDs

dge <- DGEList(counts)
dge <- calcNormFactors(dge)
# cpm_norm <- cpm(dge, normalized.lib.sizes = TRUE)
cpm_norm <- cpm(dge, normalized.lib.sizes = TRUE, log =T)
new_expr <- cpm_norm
result <- run_subtype_classifier(new_expr_data = new_expr)
head(result$predictions)
tmp = result[["predictions"]]

tmp = tmp[df1$SampleID,]
df1$CTS = tmp$CTS

new_expr_cts = new_expr[g$gene,]
df_expr = data.frame(t(new_expr_cts))

df_expr = df_expr[df1$SampleID,]

df2 = bind_cols(df1, df_expr)


myfun = function(cts = 3, mortality = 0){
  tmp = df2 %>% filter(Timepoint == 'D1', Mortality28d == mortality, CTS == cts)
  df3 = df2[df2$HumanID %in% tmp$HumanID,]
  
  friedman_test_result_list = list()
  for (g1 in g$gene) {
    outlier_ids <- df3 %>%
      group_by(Timepoint) %>%
      identify_outliers(all_of(g1)) %>%   # rstatix 提供
      filter(is.outlier) %>%
      pull(HumanID) %>%
      unique()
    df_nonoutlier =  df3 %>% filter(!HumanID %in% outlier_ids) %>% na.omit()
    
    friedman_test_result <- df_nonoutlier %>% 
      rstatix::friedman_test(as.formula(paste(g1, "~ Timepoint | HumanID")) )
    friedman_test_result_list[[g1]] = friedman_test_result
  }
  friedman_all <- dplyr::bind_rows(friedman_test_result_list, .id = "Gene")
  rownames(gene_attr) = gene_attr$gene_id
  friedman_all$SYMBOL = gene_attr[friedman_all$Gene,][["SYMBOL"]]
  # friedman_all$p
  friedman_all$FDR <- p.adjust(friedman_all$p, method = "BH")
  friedman_all %<>% arrange(FDR)
  
  mycol = RColorBrewer::brewer.pal(6,'Set3')[-2][1:3]
  names(mycol) = c('D1','D4','D7')
  
  plist = list()
  for (i in 1:nrow(friedman_all)) {
    g1 = friedman_all$Gene[i]
    outlier_ids <- df3 %>%
      group_by(Timepoint) %>%
      identify_outliers(all_of(g1)) %>%   # rstatix 提供
      filter(is.outlier) %>%
      pull(HumanID) %>%
      unique()
    df_nonoutlier =  df3 %>% filter(!HumanID %in% outlier_ids)
    
    plist[[i]] = ggplot(df_nonoutlier, 
                        aes(x = Timepoint, y = .data[[g1]], group = HumanID)) +
      geom_boxplot(aes(x = Timepoint, y = .data[[g1]], color = Timepoint), 
                   outliers = F, inherit.aes = FALSE,
                   linewidth = 0.8) +
      scale_color_manual(values = mycol) +
      geom_point(color = "grey50", alpha = 0.5, size = 1) +
      geom_line(color = "grey50", alpha = 0.2) +
      theme_bw() +
      annotate(
        "text",
        x = Inf,
        y = Inf,
        label = paste0("FDR = ", signif(friedman_all$FDR[i], 2)),
        hjust = 1.1, vjust = 1.5
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
      theme(panel.grid.minor = element_blank(),
            legend.position = 'none') +
      ylab(gene_attr$SYMBOL[gene_attr$gene_id == g1])
  }
  return(plist)
}

plist1 = myfun(cts = 1, mortality = 0)
plist2 = myfun(cts = 2, mortality = 0)
pdf(paste0("Outputs/04_cts_genes.pdf"), width = 11, height = 5)
grid.arrange(grobs = plist1, ncol = 6)
grid.arrange(grobs = plist2, ncol = 6)
dev.off()


