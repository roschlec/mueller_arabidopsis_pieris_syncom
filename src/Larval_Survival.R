# Larval growth systems
# Data analysis and plotting
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(patchwork)

# Data --------------------------------------------------------------------

pot_box_data <-
  read_csv(here("input", "larvae_weight_mortality.csv"), show_col_types = FALSE) |> 
  separate(col = replicate, 
           into = c("type", "rep"), 
           sep = "(?<=[a-z])(?=[0-9])", 
           remove = FALSE) |> 
  mutate(Npercent = (numberoflarvae - numberofdeadlarvae)/numberoflarvae * 100)

# Stats -------------------------------------------------------------------
# Survival
survival <- 
  pot_box_data |> 
  group_by(type) |> 
  summarise(
    total = sum(numberoflarvae),
    dead = sum(numberofdeadlarvae),
    alive = total - dead)

xtab <-
  as.table(rbind(
    survival$dead,
    survival$alive
  ))

dimnames(xtab) <- list(
  Survival = c("dead", "alive"),
  GrowthType = c("box", "pot")
)

fisher_survival <- 
  fisher_test(xtab) |> 
  mutate(group1 = "box", group2 = "pot", y.position = 105)

# Larval weight
# Shapiro-Wilk Normality Test
pot_box_data |> 
  shapiro_test(meanbiomass_mg)

# Levene's test for homogeneity of variances
pot_box_data |> 
  levene_test(meanbiomass_mg ~ as.factor(type))

# Wilcoxon rank sum test
test_mg <-
  pot_box_data |> 
  wilcox_test(meanbiomass_mg ~ type) |> 
  add_xy_position() |> 
  mutate(psignif = ifelse(p < 0.05, "*", "ns"))


# Plot --------------------------------------------------------------------

plt_figS3a <-
  pot_box_data |> 
  ggplot(aes(x = type, y = Npercent)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "firebrick3") +
  geom_point(position = position_jitter(width = 0.05, height = 0), size = 1) +
  theme_bw() +
  scale_x_discrete(name = "Growth system", labels = c("Magenta Box", "Pot")) +
  scale_y_continuous(name = "Larval survival [%]", limits = c(0, 110)) +
  stat_pvalue_manual(data = fisher_survival, label = "p = {p}") +
  theme(
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 12))
  
plt_figS3b <-
  pot_box_data |> 
  ggplot(aes(x = type, y = meanbiomass_mg)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "firebrick3") +
  geom_point(position = position_jitter(width = 0.05, height = 0), size = 1) +
  theme_bw() +
  scale_x_discrete(name = "Growth system", labels = c("Magenta Box", "Pot")) +
  scale_y_continuous(name = "Mean larval weight [mg]", limits = c(0, 10)) +
  stat_pvalue_manual(data = test_mg, label = "psignif", size = 5) +
  theme(
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 12))

plt_figureS3 <-
  plt_figS3a + plt_figS3b +
  plot_annotation(tag_levels = "A")


# Output ------------------------------------------------------------------

mapply(function(x) 
  ggsave(x, 
         plot = plt_figureS3, 
         dpi = 300, width = 6, height = 3),
  x = c(here("output", "FigureS3.png"), here("output", "FigureS3.eps")))
