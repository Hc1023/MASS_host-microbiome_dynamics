rm(list = ls())
library(tidyverse)
library(magrittr)
library(rstatix)
library(stringr)
library(ggplot2)

load('Inputs/1616_meta_model.rdata')

tmp = names(meta_model)
tmp = tmp[grepl('D1_mi_', tmp)]
tmp = gsub('D1_mi_', '', tmp)
pathogen_virus = tmp[11:20]
pathogen_bacf = tmp[1:10]

df = meta_model %>% dplyr::select(Mortality28d, matches("D[147]_mi_"))


res <- df %>%
  pivot_longer(
    cols = -Mortality28d,
    names_to = "pathogen",
    values_to = "value"
  ) %>%
  group_by(pathogen) %>%
  summarise(
    p_value = wilcox.test(value ~ Mortality28d, exact = FALSE)$p.value,
    d = -rstatix::cohens_d(data = dplyr::pick(Mortality28d, value), 
                          formula = value ~ Mortality28d)$effsize,
    .groups = "drop"
  )


df_heat <- res %>%
  mutate(
    Day = str_extract(pathogen, "^D[147]"),
    pathogen = str_remove(pathogen, "^D[147]_mi_")
  )
df_heat <- df_heat %>%
  select(pathogen, Day, d)

df_heat$Day <- factor(df_heat$Day, levels = c("D1", "D4", "D7"))

df_heat <- df_heat %>%
  left_join(
    res %>%
      mutate(
        Day = str_extract(pathogen, "^D[147]"),
        pathogen = str_remove(pathogen, "^D[147]_mi_"),
        sig = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          p_value < 0.1   ~ "\u00B7",
          TRUE ~ ""
        )
      ) %>%
      select(pathogen, Day, sig),
    by = c("pathogen", "Day")
  )
df_heat$label = paste0(sprintf("%.3f", df_heat$d), df_heat$sig)


order_pathogen <- df_heat %>%
  filter(Day == "D1") %>%
  arrange(d) %>%
  pull(pathogen)
df_heat$pathogen <- factor(df_heat$pathogen, levels = order_pathogen)
df_heat$type = 'Bacf'
df_heat$type[df_heat$pathogen %in% pathogen_virus] = 'Virus'

levels(df_heat$pathogen) = gsub('\\.','-',levels(df_heat$pathogen))
levels(df_heat$pathogen)[levels(df_heat$pathogen) == 'Influenza-A'] = 'Influenza A'
levels(df_heat$pathogen)[levels(df_heat$pathogen) == 'Rhinovirus-A'] = 'Rhinovirus A'
levels(df_heat$pathogen)[levels(df_heat$pathogen) == 'HHV-4'] = 'EBV'

df_heat_virus = df_heat %>% filter(type == 'Virus') %>% droplevels()

p1 = ggplot(df_heat_virus,
       aes(x = Day, y = pathogen, fill = d)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "Cohen's d"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL) +
  geom_text(aes(label = label), size = 3)

p1
df_heat_bacf = df_heat %>% filter(type != 'Virus') %>% droplevels()

p2 = ggplot(df_heat_bacf,
            aes(x = Day, y = pathogen, fill = d)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "Cohen's d"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL) +
  geom_text(aes(label = label), size = 3)
p2

pdf(paste0("Outputs/03_micobe_cohenD.pdf"), width = 4, height = 3.2)
print(p1)
print(p2)
dev.off()


#### Total ####
df = meta_model %>% dplyr::select(Mortality28d, matches("D[147]_mi2_"))

res <- df %>%
  pivot_longer(
    cols = -Mortality28d,
    names_to = "pathogen",
    values_to = "value"
  ) %>%
  group_by(pathogen) %>%
  summarise(
    p_value = wilcox.test(value ~ Mortality28d, exact = FALSE)$p.value,
    d = -rstatix::cohens_d(data = dplyr::pick(Mortality28d, value), 
                           formula = value ~ Mortality28d)$effsize,
    .groups = "drop"
  )
df_heat <- res %>%
  mutate(
    Day = str_extract(pathogen, "^D[147]"),
    pathogen = str_remove(pathogen, "^D[147]_mi2_")
  )
df_heat <- df_heat %>%
  dplyr::select(pathogen, Day, d)

df_heat$Day <- factor(df_heat$Day, levels = c("D1", "D4", "D7"))

df_heat <- df_heat %>%
  left_join(
    res %>%
      mutate(
        Day = str_extract(pathogen, "^D[147]"),
        pathogen = str_remove(pathogen, "^D[147]_mi2_"),
        sig = case_when(
          p_value < 0.001 ~ "***",
          p_value < 0.01  ~ "**",
          p_value < 0.05  ~ "*",
          p_value < 0.1   ~ "\u00B7",
          TRUE ~ ""
        )
      ) %>%
      dplyr::select(pathogen, Day, sig),
    by = c("pathogen", "Day")
  )
df_heat$label = paste0(sprintf("%.3f", df_heat$d), df_heat$sig)


order_pathogen <- rev(c('Total', 'Bacf', 'Viruses'))
df_heat$pathogen <- factor(df_heat$pathogen, levels = order_pathogen)
levels(df_heat$pathogen)[2] = "Bacteria/Fungi"
p1 = ggplot(df_heat,
            aes(x = Day, y = pathogen, fill = d)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = ""
  ) +
  theme_classic() +
  theme(
    legend.position = 'right',
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(x = NULL, y = NULL) +
  geom_text(aes(label = label), size = 3)

p1

pdf(paste0("Outputs/03_micobe_cohenD_total.pdf"), width = 3.6, height = 1.5)
print(p1)
dev.off()



