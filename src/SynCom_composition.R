# SynCom Composition
# Plotting
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(forcats)
library(ggtext)

# Data --------------------------------------------------------------------

##  Species labels (for element_markdown)
sp_label <-
  c("Methylobacterium radiotolerans 0-1" = "*Methylobacterium radiotolerans* 0-1",
    "Pantoea eucalypti 299R" = "*Pantoea eucalypti* 299R",
    "Methylobacterium sp. AC832" = "*Methylobacterium* sp. AC832",
    "Sphingomonas melonis FR1" = "*Sphingomonas melonis* FR1",
    "Arthrobacter sp. Leaf145" = "*Arthrobacter* sp. Leaf145",
    "Microbacterium sp. AC026" = "*Microbacterium* sp. AC026",
    "Williamsia sp. Leaf354" = "*Williamsia* sp. Leaf354",
    "Bradyrhizobium sp. Leaf396" = "*Bradyrhizobium* sp. Leaf396",
    "Luteimonas sp. AC640" = "*Luteimonas* sp. AC640",
    "Bacillus sp. AC266" = "*Bacillus* sp. AC266",
    "Curtobacterium sp. AC273" = "*Curtobacterium* sp. AC273",
    "Flavobacterium sp. AC267" = "*Flavobacterium* sp. AC267",
    "Massilia sp. AC133" = "*Massilia* sp. AC133",
    "Pseudomonas sp. AC369" = "*Pseudomonas* sp. AC369",
    "Pseudomonas syringae B728a" = "*Pseudomonas syringae* B728a",
    "Methylobacterium sp. AC433" = "*Methylobacterium* sp. AC433",
    "Sphingomonas sp. AC435" = "*Sphingomonas* sp. AC435",
    "Sphingomonas sp. Leaf357" = "*Sphingomonas* sp Leaf357",
    "Pseudomonas sp. AC021" = "*Pseudomonas* sp. AC021",
    "Rhodococcus sp. Leaf225" = "*Rhodococcus* sp. Leaf225")
  
sp.df <-
  as.data.frame(sp_label) |> 
  rownames_to_column(var = "Strain")

syncom.df <- read.csv(here("input", "syncom_composition.csv"))


# Plot --------------------------------------------------------------------

plt_figureS1 <-
  syncom.df |> 
  mutate(present = 1) |> 
  pivot_wider(id_cols = c(experiment, SynCom), 
              names_from = Strain, 
              values_from = present,
              values_fill = 0) |> 
  pivot_longer(cols = c(-experiment, -SynCom), names_to = "Strain") |> 
  mutate(SynCom = factor(SynCom, 
                         levels = c("Syncom5", "Syncom10", "SynCom20"),
                         labels = c("SynCom5", "SynCom10", "SynCom20"))) |> 
  inner_join(sp.df, by = "Strain") |> 
  ggplot(aes(x = SynCom, y = sp_label, fill = as.factor(value))) +
  facet_wrap(~ experiment, ncol = 4) +
  geom_tile(color = "black") +
  scale_fill_manual(name = "SynCom\nmember", labels = c("No", "Yes"), values = c("white", "grey50")) +
  theme(
    panel.background = element_blank(),
    strip.background = element_blank(),
    axis.text = element_text(size = 8),
    axis.text.y = element_markdown(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_blank(),
    axis.ticks = element_blank())

# Output ------------------------------------------------------------------

mapply(function(x) 
  ggsave(x, 
         plot = plt_figureS1, 
         dpi = 300, width = 8, height = 3),
  x = c(here("output", "FigureS1.png"), here("output", "FigureS1.eps")))
