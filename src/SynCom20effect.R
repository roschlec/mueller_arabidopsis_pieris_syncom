# SynCom20 abundance on fed leaves ----------------------------------------
# Statistics and plotting
# Authors: Moritz Müller, Rudolf Schlechter

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(patchwork)

# Data input --------------------------------------------------------------

larval_weight <-
  read_delim(here("input", "larval_weight_syncom20.txt"), delim = "\t", show_col_types = FALSE)
  
cfu <- read_delim(here("input", "syncom20_leaf.txt"), delim = "\t", show_col_types = FALSE)

# Stats -------------------------------------------------------------------
# Larval weight

# N replicates
larval_weight |> 
  group_by(treatment) |> 
  tally()

# Shapiro-Wilk Normality test
larval_weight |> 
  group_by(treatment) |> 
  shapiro_test(weightperlarva)

# Levene test for homogeneity of variances
larval_weight |>  
  levene_test(weightperlarva ~ as.factor(treatment))

# Student's t test
larval_weight_stat <-
  larval_weight |> 
  t_test(weightperlarva ~ treatment, paired = FALSE, var.equal = TRUE, ref.group = "control") |> 
  add_xy_position() |> 
  mutate(p.signif = ifelse(p < 0.05, "*", "ns"))


## CFU data
# Shapiro-Wilk Normality test
cfu |> 
  group_by(treatment) |> 
  shapiro_test(logcfu)

# Levene test for homogeneity of variances
cfu |> 
  levene_test(logcfu ~ as.factor(treatment))

# Student's t test
cfu_stat <-
  cfu |> 
  t_test(logcfu ~ treatment, paired = FALSE, var.equal = TRUE, ref.group = "SynCom20") |> 
  add_xy_position() |> 
  mutate(p = format(p, scientific = FALSE),
         p.signif = ifelse(p < 0.05, "*", "ns"))

# Plot --------------------------------------------------------------------

plt_fig5a <-
  larval_weight |> 
  ggplot(aes(x = treatment, y = weightperlarva)) +
  geom_boxplot(fill = "gray90", width = 0.4) +
  geom_point(position = position_jitter(width = 0.05)) +
  theme_light() +
  labs( x = "Treatment", y = "Larval weight (mg)")+
  scale_x_discrete(limit = c('control', 'SynCom20'),
                   labels = c('Control', 'SynCom20'))+
  stat_pvalue_manual(data = larval_weight_stat, label = "p.signif", size = 5) +
  ylim(0, 32) +
  theme(
    text = element_text(size = 10, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank()
  )

plt_fig5b <-
  cfu |> 
  ggplot(aes(x = treatment, y = logcfu)) +
  geom_boxplot(fill = "gray90", width = 0.4) +
  geom_point(position = position_jitter(width = 0.05)) +
  theme_light() +
  labs( x = "Treatment", y = expression('Log'[10]~'(CFU/g leaf)'))+
  scale_x_discrete(limit = c('SynCom20', 'SynCom20 + Larvae'),
                   labels = c('SynCom20\n- Feeding', 'SynCom20\n+ Feeding'))+
  stat_pvalue_manual(data = cfu_stat, label = "p.signif", size = 5) +
  ylim(7, 9.5) +
  theme(
    text = element_text(size = 10, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank()
  )

plt_fig5 <-
  plt_fig5a + plt_fig5b +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 10, colour = "black"),
    axis.text = element_text(size = 10, colour = "black"),
    axis.title.x = element_blank()
  )


# Output ------------------------------------------------------------------

mapply(function(x) 
  ggsave(x, 
         plot = plt_fig5, 
         dpi = 300, width = 6, height = 3),
  x = c(here("output", "Figure5.png"), here("output", "Figure5.eps")))
