# Effect of SynCom richness to larval performance
# Plotting
# Authors: Moritz Müller and Rudolf Schlechter


# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(smplot2)

# Function ----------------------------------------------------------------

plot_corr <- function(df, x_name, fill_name){
  min_x <- df |> summarise(min = min({{x_name}})) |> pull(min)
  
  df |> 
    ggplot(aes(x = {{x_name}}, y = normweight)) +
    geom_point(shape = 21, size = 3, color = "white", aes(fill = {{fill_name}}))+
    sm_statCorr(corr_method = "spearman", label_x = min_x, color = "black")
}

# Data Input --------------------------------------------------------------

source(here("src", "SynCom_richness.R"))

# Plots -------------------------------------------------------------------

# Larval survival
plt_figure3a <- 
  larvae_survival |> 
  ggplot(aes(x = treatment, y = Npercent)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "firebrick3") +
  geom_point(aes(fill = experiment_rep),
             pch = 21, size = 1.5, alpha = 0.8,
             position = position_jitter(width = 0.1, height = 0))+
  theme_light()+
  ylim(0, 105)+
  labs(x = "Treatment", y = "Larval survival [%]") +
  scale_x_discrete(limit = c('control', 'SynCom5', 'SynCom10', 'SynCom20'),
                   label = c('Control', 'SynCom5', 'SynCom10', 'SynCom20')) +
  theme(
    text = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  ) +
  guides(fill="none")

# CFU on Leaves
plt_figure3b <-
  weight |> 
  ggplot(aes(x = treatment, y = logleafcfu)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2, fill = "gray90")+
  geom_point(aes(fill = experiment_rep),
             pch = 21, size = 1.5, alpha = 1,
             position = position_jitter(width = 0.1, height = 0))+
  theme_light()+
  labs( x = "Treatment", y = expression('Log'[10]~'(CFU/g leaf)'))+
  scale_x_discrete(limit = c('control', 'SynCom5', 'SynCom10', 'SynCom20'),
                   labels = c('Control', 'SynCom5', 'SynCom10', 'SynCom20')) +
  theme(legend.position="none")+
  geom_text(data = cld_leafcfu, aes(x = treatment, y = 9.5, label = .group))

# CFU in Larvae
plt_figure3c <- 
  weight |> 
  ggplot(aes(x = treatment, y = loglarvaecfu)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2, fill = "gray90")+
  geom_point(aes(fill = experiment_rep),
             pch = 21, size = 1.5, alpha = 1,
             position = position_jitter(width = 0.1, height = 0))+
  theme_light()+
  labs( x = "Treatment", y = expression('Log'[10]~'(CFU/g larva)'))+
  scale_x_discrete(limit = c('control', 'SynCom5', 'SynCom10', 'SynCom20'),
                   labels = c('Control', 'SynCom5', 'SynCom10', 'SynCom20')) +
  theme(legend.position="none")+
  geom_text(data = cld_lvcfu, aes(x = treatment, y = 9.5, label = .group))

# Larval weight
plt_figure3d <-
  weight |> 
  ggplot(aes(x = treatment, y = normweight)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.2, fill = "gray90")+
  geom_point(aes(fill = experiment_rep),
             pch = 21, size = 1.5, alpha = 1,
             position = position_jitter(width = 0.1, height = 0))+
  theme_light()+
  labs( x = "Treatment", y = "Normalised larval\nweight to control")+
  scale_x_discrete(limit = c('control', 'SynCom5', 'SynCom10', 'SynCom20'),
                   labels = c('Control', 'SynCom5', 'SynCom10', 'SynCom20')) +
  scale_y_continuous(breaks = seq(0.4, 1.8, 0.2)) +
  geom_text(data = cld_lv_weight, aes(x = treatment, y = 1.7, label = .group)) +
  guides(fill = guide_legend(title = "Experiment", keyheight = unit(0.15, 'inch')))


# Combine plots -----------------------------------------------------------

plt_figure3 <-
  plt_figure3a + plt_figure3b + plt_figure3c + plt_figure3d +
  plot_layout(guides = "collect", ncol = 4) +
  plot_annotation(tag_levels = "A") &
  theme(
    text = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank()
  )

## Supplemental

# Correlation richness and normalised larval weight
plt_figureS4a <-
  plot_corr(weight, BacRichness, treatment) +
  labs(x = "Richness") +
  scale_fill_manual(labels = c("Control", "SynCom5", "SynCom10", "SynCom20"),
                    values = c("indianred1","forestgreen","purple","dodgerblue"))
  

# Correlation Faith's PD and normalised larval weight
plt_figureS4b <-
  plot_corr(weight, FaithsPD, treatment) +
  labs(x = "Faith's PD")+
  scale_fill_manual(labels = c("Control", "SynCom5", "SynCom10", "SynCom20"),
                    values = c("indianred1","forestgreen","purple","dodgerblue"))

# Correlation Larval weight and CFU in larvae
plt_figureS4c <-
  plot_corr(weight_onlySynCom, loglarvaecfu, treatment) + 
  labs(x = expression(paste(Log[10], " CFU/g leaf"))) +
  scale_fill_manual(labels = c("SynCom5", "SynCom10", "SynCom20"), 
                    values = c("forestgreen","purple","dodgerblue"))

# Correlation Larval weight and CFU on leaves
plt_figureS4d <-
  plot_corr(weight_onlySynCom, logleafcfu, treatment) + 
  labs(x = expression(paste(Log[10], " CFU/g larva"))) +
  scale_fill_manual(labels = c("SynCom5", "SynCom10", "SynCom20"), 
                    values = c("forestgreen","purple","dodgerblue"))

plt_figureS4 <-
  (plt_figureS4a + plt_figureS4b + plt_figureS4c + plt_figureS4d) +
  plot_layout(ncol = 4) +
  plot_annotation(tag_levels = "A") &
  labs(y = "Normalised larval weight")

# Output ------------------------------------------------------------------

mapply(function(x) 
  ggsave(x, 
         plot = plt_figure3, 
         dpi = 300, width = 12, height = 3.5),
  x = c(here("output", "Figure3.png"), here("output", "Figure3.eps")))

mapply(function(x) 
  ggsave(x, 
         plot = plt_figureS4, 
         dpi = 300, width = 12, height = 3.5),
  x = c(here("output", "FigureS4.png"), here("output", "FigureS4.eps")))
