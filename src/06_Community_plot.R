# Relative abundance analysis ---------------------------------------------

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(forcats)
library(ggtext)
library(cowplot)
library(colorspace)
library(patchwork)
library(ggrepel)
library(ggpp)

# Global options ----------------------------------------------------------
# Color palette
c25 <- c("dodgerblue2", "#FF7F00", "orchid1", "#E31A1C", "#6A3D9A", "black", "gold1",
         "skyblue2", "#FB9A99", "palegreen2", "#CAB2D6", "#FDBF6F","gray70", "khaki2",
         "maroon", "green4", "deeppink1", "blue1", "steelblue4", "darkturquoise", 
         "green1", "yellow4", "yellow3","darkorange4", "brown", "red", "yellow", "lightblue")

# Species labels (for element_markdown)
sp_label <-
  c("Methylobacterium radiotolerans" = "*Methylobacterium radiotolerans* 0-1",
    "Pantoea eucalypti" = "*Pantoea eucalypti* 299R",
    "Methylobacterium AC832" = "*Methylobacterium* sp. AC832",
    "Sphingomonas melonis" = "*Sphingomonas melonis* FR1",
    "Arthrobacter Leaf145" = "*Arthrobacter* sp. Leaf145",
    "Microbacterium AC026" = "*Microbacterium* sp. AC026",
    "Williamsia Leaf354" = "*Williamsia* sp. Leaf354",
    "Bradyrhizobium Leaf396" = "*Bradyrhizobium* sp. Leaf396",
    "Agrococcus AC640" = "*Luteimonas* sp. AC640",
    "Bacillus AC266" = "*Bacillus* sp. AC266",
    "Curtobacterium AC273" = "*Curtobacterium* sp. AC273",
    "Flavobacterium AC267" = "*Flavobacterium* sp. AC267",
    "Massilia AC133" = "*Massilia* sp. AC133",
    "Pseudomonas AC369" = "*Pseudomonas* sp. AC369",
    "Pseudomonas syringae" = "*Pseudomonas syringae* B728a",
    "Methylobacterium AC433" = "*Methylobacterium* sp. AC433",
    "Sphingomonas AC435" = "*Sphingomonas* sp. AC435",
    "Sphingomonas Leaf357" = "*Sphingomonas* sp Leaf357",
    "Pseudomonas AC021" = "*Pseudomonas* sp. AC021")

# Data Input --------------------------------------------------------------
# Relative abundance
source(here("src", "05_RelativeAbundance.R"))

# NMDS
nmds_df <- readRDS(here("input", "processed", "nmds.rds"))
nmds_fit <- readRDS(here("input", "processed", "fit.nmds.rds")) |> 
  mutate(genus_species = paste(genus, species, sep = " ")) |> 
  left_join(as.data.frame(sp_label) |> rownames_to_column(var = "genus_species"), 
            by = "genus_species")

# Differential abundance data
diff <- readRDS(here("input", "processed", "diffential.abundance.rds"))

label_contrast <- 
  c(feeding = "Leaf + Feeding\nvs\nLeaf - Feeding",
    filtering = "Larvae\nvs\nLeaf + Feeding",
    contrast_larva_nofeeding = "Larvae\nvs\nLeaf - Feeding")

# Plot Relative Abundance -------------------------------------------------

plt_fig6a1 <-
  df_syncom20_trt |> 
  # Plot
  ggplot(aes(x = sample, y = freq, fill = fct_reorder(genus_species, -freq))) +
  facet_wrap(~ sampleType + Feeding, nrow = 1, scale = "free_x") +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c25, labels = sp_label)+
  guides(fill = guide_legend(title = "Genus species", ncol = 2,
                             override.aes = list(stroke = 0))) +
  scale_y_continuous(lim = c(0, 1.01), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_minimal()+
  labs(x = "Replicate", y = "Relative abundance [%]")+
  theme(
    panel.border = element_rect(fill = NA, linewidth = 0.5),
    text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length.y = unit(1, "mm"),
    strip.background = element_rect(linewidth = 0.2, colour = "black"),
    legend.key.spacing.y = unit(0.1, "mm"),
    legend.key.size = unit(3, "mm"),
    legend.title = element_text(size = 10),
    legend.text = element_markdown(size = 8),
    legend.background = element_rect(linewidth = 0.2, colour = "black"),
    legend.box.spacing = unit(1, "mm"))

plt_fig6a <-
  ggdraw(plt_fig6a1) + 
  draw_label("Sample Type", x = 0.73, y = 0.94, size = 10, hjust = 0) +
  draw_label("Feeding", x = 0.73, y = 0.87, size = 10, hjust = 0)

# Plot NMDS ---------------------------------------------------------------

nmds_fit <- 
  nmds_fit |> 
  mutate(strain_label = case_when(
    species == "eucalypti" ~ "299R",
    species == "melonis" ~ "FR1",
    species == "radiotolerans" ~ "0-1",
    TRUE ~ species
  ))

plt_fig6b <-
  nmds_df |> 
  ggplot() +
  geom_point(aes(NMDS1, NMDS2, fill = Condition), 
             size = 3, stroke = 0.5, pch = 21) +
  geom_segment(data = nmds_fit, 
               aes(x = 0, xend = 0.5*NMDS1, y = 0, yend = 0.5*NMDS2),
               arrow = grid::arrow(angle = 20, length = unit(0.2, "cm")),
               alpha = 0.5) +
  geom_text_repel(data = nmds_fit, aes(x = 0.5*NMDS1, y = 0.5*NMDS2, label = strain_label),
                  xlim = c(NA, Inf),
                  ylim = c(-Inf, Inf),
                  seed = 41, 
                  box.padding = 0.3,
                  min.segment.length = Inf, size = 3) +
  theme_bw() +
  lims(x = c(-1, 1.3), y = c(-0.9, 1)) +
  scale_fill_discrete_qualitative(palette = "Dark3",
                                  labels = c("Larvae",
                                             "Leaf + Feeding",
                                             "Leaf - Feeding")) +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 0.5),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    legend.position.inside = c(0.8, 0.2),
    legend.key.spacing.y = unit(0.1, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.background = element_rect(linewidth = 0.2, colour = "black"),
    legend.box.spacing = unit(1, "mm"),
    legend.position = "bottom"
  )

# Plot Differential Abundance ---------------------------------------------

plt_fig6c <-
  diff |>
  mutate(Contrast = factor(Contrast, 
                           levels = c("contrast_feeding", 
                                      "contrast_larva_fedLeaf", 
                                      "contrast_larva_nonfedLeaf"))) |> 
  mutate(
    label_contrast1 = case_when(
      Contrast == "contrast_feeding" ~ "Leaf - Feeding",
      Contrast == "contrast_larva_fedLeaf" ~ "Leaf + Feeding",
      Contrast == "contrast_larva_nonfedLeaf" ~ "Leaf - Feeding"),
    label_contrast2 = case_when(
      Contrast == "contrast_feeding" ~ "Leaf + Feeding",
      Contrast == "contrast_larva_fedLeaf" ~ "Larvae",
      Contrast == "contrast_larva_nonfedLeaf" ~ "Larvae")
  ) |>
  ggplot(aes(x = log2FoldChange, 
             y = fct_reorder(genus_species, log2FoldChange), 
             fill = padj)) +
  facet_wrap(~Contrast,
             labeller = labeller(Contrast = label_contrast)) +
  geom_segment(aes(xend = 0), linewidth = 0.5)+
  geom_vline(xintercept = 0) + 
  geom_point(size = 2.5, pch = 21) +
  geom_text(aes(x = -6, y = 20, label = label_contrast1), 
            colour = "black", alpha = 1, vjust = "center", hjust = "center", size = 3, check_overlap = TRUE) +
  geom_text(aes(x = 6, y = 20, label = label_contrast2), 
            colour = "black", alpha = 1, vjust = "center", hjust = "center", size = 3, check_overlap = TRUE) +
  coord_cartesian(ylim = c(1,20), clip = "off") +
  guides(fill = guide_colourbar(title = expression(paste(italic("P"), "-adjusted")))) +
  scale_fill_continuous_sequential(palette = "Grays", 
                                     begin = 1,
                                     end = 0,
                                     breaks = c(0.05, 0.5, 1.0),
                                     labels = c("0.05", "0.5", "1.0")) +
  scale_y_discrete(name = NULL, labels = sp_label) + 
  scale_x_continuous(name = expression(paste("Relative abundance (",log[2], "FC)")),
                     limits = c(-12, 12), breaks = seq(-8, 9, 4)) +
  theme_bw() +
  theme(
    panel.border = element_rect(fill = NA, linewidth = 0.3),
    text = element_text(size = 12, color = "black"),
    axis.text.y = element_markdown(size = 9, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_blank(),
    legend.key.spacing.y = unit(0.1, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.background = element_rect(linewidth = 0.2, colour = "black"),
    legend.box.spacing = unit(5, "mm"),
    legend.position = "bottom",
    legend.margin = margin(t = 5, r = 15, b = 5, l = 10)
  )

# Combine plots -----------------------------------------------------------

plt_fig6 <-
  free(plt_fig6a) / wrap_plots(free(plt_fig6b), plt_fig6c, widths = c(1, 1)) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")

# Output ------------------------------------------------------------------

mapply(function(x) 
  ggsave(x, 
         plot = plt_fig6, 
         dpi = 300, width = 15, height = 8),
  x = c(here("output", "Figure6.png"), here("output", "Figure6.eps")))
