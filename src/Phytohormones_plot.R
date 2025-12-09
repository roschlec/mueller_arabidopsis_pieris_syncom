# Phytohormone analysis
# Plotting
# Authors: Moritz Müller and Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(patchwork)
library(ggh4x)

# Data --------------------------------------------------------------------

source(here("src", "Phytohormones.R"))

# Plot --------------------------------------------------------------------

plt_figure4 <-
  hormone_plt_df |> 
  ggplot(aes(x = inoculated, y = norm_fw, fill = fed)) +
  facet_wrap(~hormone, scales = "free_y", axes = "all_x") +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(size = 1, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  labs(x = "SynCom inoculation", y = "Relative phytohormone change to T0") +
  geom_text(aes(label = .group, y = ypos), 
            position = position_dodge(width = 0.75), check_overlap = TRUE) +
  guides(fill = guide_legend(override.aes = list(linetype = 0, shape = 22, size = 6, colour = "black", alpha = 1))) +
  scale_fill_manual(name = "Larval feeding", values = c("#469990", "#fffac8"), labels = c("No", "Yes")) +
  theme_light() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, colour = "black", hjust = 0),
    text = element_text(size = 12, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    legend.background = element_rect(fill = "transparent")
  )


# Supplemental ------------------------------------------------------------

plt_figureS <-
  phytohormone.df |> 
  pivot_longer(cols = c(JA_ng.fw, JA_Ile_ng.fw, ABA_ng.fw, SA_ng.fw), 
               names_to = "hormone") |> 
  mutate(hormone = factor(hormone, 
                          levels = c("JA_ng.fw", "JA_Ile_ng.fw", "ABA_ng.fw", "SA_ng.fw"),
                          labels = c("JA", "JA-Ile", "ABA", "SA")),
         inoculated = factor(inoculated,
                             levels = c("NS", "S"),
                             labels = c("Control", "SynCom20")),
         fed = factor(fed,
                      levels = c("NF", "F"),
                      labels = c("Non-fed", "Fed"))) |> 
  select(run, timepoint, fed, inoculated, hormone, value) |> 
  ggplot(aes(x = inoculated, y = value, fill = fed)) +
  facet_nested(run ~ hormone + timepoint, scales = "free_y", independent = "y") +
  geom_boxplot(outlier.alpha = 0) +
  geom_point(size = 0.8, 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75)) +
  labs(x = "SynCom inoculation", y = "Phytohormone level [ng/g leaf]") +
  guides(fill = guide_legend(override.aes = list(linetype = 0, shape = 22, size = 6, colour = "black", alpha = 1))) +
  scale_fill_manual(name = "Larval feeding", values = c("#469990", "#fffac8"), labels = c("No", "Yes")) +
  theme_light() +
  theme(
    strip.background = element_rect(fill = NA, linewidth = 0.4, colour = "grey"),
    strip.text = element_text(size = 12, colour = "black", hjust = 0),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    legend.background = element_rect(fill = "transparent")
  )


# Output ------------------------------------------------------------------

mapply(function(x) 
  ggsave(x, 
         plot = plt_figure4, 
         dpi = 300, width = 12, height = 5),
  x = c(here("output", "Figure4.png"), here("output", "Figure4.eps")))

#mapply(function(x) 
#  ggsave(x, 
#         plot = plt_figureS, 
#         dpi = 300, width = 12, height = 5),
#  x=c(here("output", "FigureS_hormone.png"), here("output", "FigureS_hormone.eps")))
