# Alpha diversity metrics
# Data processing, analysis and plotting
# Authors: Moritz Müller, Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(forcats)
library(ggtext)
library(ggsignif)
library(vegan)
library(rstatix)
library(ggpubr)
library(ggtree)
library(picante)
library(cld)

# Data input --------------------------------------------------------------

syncom20.df <- 
  readRDS(here("input", "processed", "otu_df.rds")) |>
  remove_rownames() |> 
  column_to_rownames(var = "species") |> 
  select(-phylum:-genus)

#   Check sample coverage
df_long <-
  syncom20.df |> 
  rownames_to_column(var = 'OTU') |> 
  pivot_longer(cols = -OTU, names_to = "sample")

sample_coverage <-
  df_long  |>  
  group_by(sample)  |>  
  summarise(n_seqs = sum(value))

#   Rarefying
##  seq limit used as alternative too goods coverage limit as there are too many low seq, 
##  high goods coverage samples (only controls fall out)
keep_samples <- 
  sample_coverage |> 
  filter(n_seqs >= 300) |> 
  pull(sample)

filter_df <- t(syncom20.df[, keep_samples])

min_seq <- 
  syncom20.df[, keep_samples] |> 
  pivot_longer(everything(), names_to = "sample") |> 
  group_by(sample) |> 
  summarise(n_seqs = sum(value)) |> 
  summarise(min = min(n_seqs)) |> 
  pull(min)

##  Diversity metrics
rarefied_list <- list()

for(i in 1:100){
  rarefied_list[[i]] <- 
    rrarefy(filter_df, sample = min_seq) |>
    as_tibble(rownames = "sample") |> 
    pivot_longer(-sample)
}

rarefied_data <-
  rarefied_list |> 
  bind_rows(.id = "loop") |> 
  group_by(sample, name) |> 
  summarise(mean = mean(value),
            sd = sd(value),
            .groups = "drop") |> 
  mutate(value = round(mean, 0))

count_table_rarefied <-
  rarefied_data |> 
  pivot_wider(id_cols = sample, 
              names_from = name,
              values_from = value,
              values_fill = 0) |> 
  column_to_rownames(var = "sample")

#  Alpha diversity
filter_df1 <- as.data.frame(filter_df)

## Determining statistical significance for Alpha Diversity measurements

alphadivstats <-
  filter_df1 |> 
  rownames_to_column(var = "sample") |> 
  rowwise() |> 
  mutate(
    Richness = specnumber(c_across(-sample)),
    Shannon = diversity(c_across(-sample), index = "shannon"),
    Evenness = Shannon/log(Richness)) |> 
  ungroup() |>
  select(sample, Richness, Shannon, Evenness) |> 
  separate(sample, 
           into = c("Condition", "rep"), 
           sep = "_",
           remove = FALSE) |> 
  mutate(Treatment = factor(Condition,
                            levels = c("EC", "IN", "SF", "SNF", "LS"),
                            labels = c("Kitome", "Inoculum", "Leaf + Feeding", "Leaf - Feeding", "Larvae"))) |>
  filter(Treatment != "Kitome")


# MDP ---------------------------------------------------------------------

#Calculation of mean pairwise phylogenetic distance per sample
##Clean up Dataset
forMPD <- syncom20.df
forMPD <- forMPD[, -c(1:26)]
forMPD <- forMPD[, -c(9:11)]
forMPD <- t(forMPD) 
forMPD <- as.data.frame(forMPD)

##Import tree and rename nodelabels
My.tree <- 
  read.tree(here("input", "SpeciesTree_rooted_node_labelsmanual.txt"))
My.tree$tip.label<- c("AC266","Leaf225","Leaf354","Leaf145","AC273","AC026","AC267","radiotolerans","AC433","AC832","Leaf396","Leaf357","melonis","AC435","AC133","eucalypti","AC369","AC021","syringae","AC640")

##Calculate MPDs
forMPD$MPDscore <- mpd(forMPD, cophenetic(My.tree), abundance.weighted = TRUE)

###Normalisation of MPD of samples by MPD of homogeneous S20 community
forMPD$MPDscore <- (forMPD$MPDscore)/1.047995

#Mean pairwise phylogenetic difference of bacterial communities grouped by treatment
##rename treatments
forMPD$Treatment <- rownames(forMPD)
forMPD$Treatment[startsWith(forMPD$Treatment, 'EC') == TRUE] <- 'Kitome'
forMPD$Treatment[startsWith(forMPD$Treatment, 'IN') == TRUE] <- 'Inoculum'
forMPD$Treatment[startsWith(forMPD$Treatment, 'SF') == TRUE] <- 'Leaf + Feeding'
forMPD$Treatment[startsWith(forMPD$Treatment, 'SNF') == TRUE] <- 'Leaf - Feeding'
forMPD$Treatment[startsWith(forMPD$Treatment, 'LS') == TRUE] <- 'Larvae'


# Combine Diversity Metrics -----------------------------------------------

diversity_metrics <-
  forMPD |> 
  filter(Treatment != "Kitome") |> 
  rownames_to_column(var = "sample") |> 
  left_join(alphadivstats, forMPD, by = c("sample", "Treatment")) |> 
  pivot_longer(c(MPDscore, Richness, Shannon, Evenness), names_to = "metric") |> 
  select(Treatment:value) |> 
  mutate(metric = factor(metric, levels = c("Richness", "Evenness", "Shannon", "MPDscore")))


# Stats -------------------------------------------------------------------
##  Descriptive
diversity_metrics |> 
  group_by(metric, Treatment) |> 
  summarise(mean = mean(value),
            sd = sd(value))

## stat test alpha diversity
##  Normality
diversity_metrics |> 
  filter(Treatment != "Inoculum") |> 
  group_by(metric, Treatment) |> 
  shapiro_test(value) |> 
  mutate(pass = ifelse(p > 0.05, "Y", "N"))

##  HETEROSCEDASTICITY
diversity_metrics |> 
  group_by(metric) |> 
  levene_test(value ~ factor(Treatment)) |> 
  mutate(pass = ifelse(p > 0.05, "Y", "N"))

##  ANOVA
diversity_metrics |> 
  group_by(metric) |> 
  anova_test(value ~ factor(Treatment), detailed = TRUE) 

##  TUKEY POST HOC
diversity_stat <-
  diversity_metrics |> 
  group_by(metric) |> 
  tukey_hsd(value ~ Condition) |> 
  add_xy_position(stack = TRUE) |> 
  mutate(y.position = case_when(metric == "Richness" ~ 19, 
                                metric == "Shannon" ~ 1.8, 
                                metric == "MPDscore" ~ 0.5, 
                                metric == "Evenness" ~ 0.8))

list_metrics <- split(diversity_metrics, diversity_metrics$metric)

list_tukey <- lapply(list_metrics, tukey_hsd, value ~ Condition)

list_cld <- lapply(list_tukey, cld::make_cld, reversed = TRUE)

diversity_cld <- 
  bind_rows(list_cld, .id = "metric") |> 
  as.data.frame() |> 
  select(metric, Condition = group, cld) |> 
  mutate(y.position = case_when(metric == "Richness" ~ 20, 
                                metric == "Shannon" ~ 2.6, 
                                metric == "MPDscore" ~ 0.9, 
                                metric == "Evenness" ~ 0.9),
         metric = factor(metric, levels = levels(diversity_metrics$metric))) |> 
  left_join(diversity_metrics |> select(Treatment, Condition) |> unique(), 
            by = "Condition") |> 
  mutate(Treatment = factor(Treatment,
                            levels = c("Inoculum", "Leaf + Feeding", "Leaf - Feeding", "Larvae")))


# Plot --------------------------------------------------------------------

metric_lab <- c("Richness",
                "Evenness",
                "Shannon Diversity",
                "Mean Pairwise Distance")

names(metric_lab) <- levels(diversity_metrics$metric)

plt_figure7 <-
  diversity_metrics |> 
  ggplot(aes(x = Treatment, y = value)) +
  facet_wrap(~metric, ncol = 4, scales = "free_y", axes = "all_x", labeller = labeller(metric = metric_lab)) +
  geom_boxplot(outlier.alpha = 0, fill = "gray90") +
  geom_point(position = position_jitter(width = 0.05, height = 0)) +
  geom_text(data = diversity_cld, aes(y = y.position, label = cld)) +
  labs(y = "Diversity Index") +
  theme_light() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 10, colour = "black", hjust = 0),
    text = element_text(size = 10, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title.x = element_blank()
  )


# Output ------------------------------------------------------------------

mapply(function(x) 
  ggsave(x, 
         plot = plt_figure7, 
         dpi = 300, width = 10, height = 4),
  x = c(here("output", "Figure7.png"), here("output", "Figure7.eps")))
