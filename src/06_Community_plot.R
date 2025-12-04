# Relative abundance analysis ---------------------------------------------

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(forcats)
library(ggtext)
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

plt_fig6a <-
  df_syncom20_trt |> 
  # Plot
  ggplot(aes(x = sample, y = freq, fill = fct_reorder(genus_species, -freq))) +
  facet_wrap(~ sampleType + SynCom + Feeding, nrow = 1, scale = "free_x") +
  geom_col(position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c25, labels = sp_label)+
  guides(fill = guide_legend(title = "Genus species", ncol = 2,
                             override.aes = list(stroke = 0))) +
  scale_y_continuous(limits = c(0, 1.01), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_minimal()+
  labs(x = "Replicate", y = "Relative abundance [%]")+
  theme(
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

# Plot NMDS ---------------------------------------------------------------

plt_fig6b <-
  nmds_df |> 
  ggplot() +
  geom_point(aes(NMDS1, NMDS2, fill = Condition), 
             size = 3, stroke = 0.5, pch = 21) +
  geom_segment(data = nmds_fit, 
               aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
               arrow = grid::arrow(angle = 20, length = unit(0.2, "cm")),
               alpha = 0.5) +
  geom_richtext(data = nmds_fit, aes(x = NMDS1, y = NMDS2, label = sp_label),
               fill = NA, label.color = NA, size = 3,
               position = position_nudge_center(0.3, 0.01, 0, 0)) +
  theme_bw() +
  scale_fill_discrete_qualitative(palette = "Dark3",
                                  labels = c("Larvae",
                                             "Leaf + Feeding",
                                             "Leaf - Feeding")) +
  theme(
    text = element_text(size = 12, color = "black"),
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
                           levels = c("contrast_feeding", "contrast_larva_fedLeaf", "contrast_larva_nonfedLeaf"))) |> 
  ggplot(aes(x = log2FoldChange, 
             y = fct_reorder(genus_species, log2FoldChange), 
             alpha = padj)) +
  facet_wrap(~Contrast,
             labeller = labeller(Contrast = label_contrast)) +
  geom_vline(xintercept = 0) + 
  geom_point(size = 2.5) +
  geom_segment(aes(xend = 0), linewidth = 1)+
  guides(colour = "none",
         alpha = guide_legend(title = expression(paste(italic("P"), "-adjusted")),
                              override.aes = list(linewidth = 0))) +
  scale_alpha_continuous(range = c(1, 0.2), 
                         breaks = c(0.001, 0.05, 0.5, 1),
                         labels = c("<0.001", "0.05", "0.5", "1.0")) +
  scale_y_discrete(name = NULL, labels = sp_label) + 
  scale_x_continuous(name = expression(paste("Relative abundance (",log[2], "FC)")),
                     limits = c(-11, 9), breaks = seq(-8, 9, 4)) +
  theme_bw() +
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text.y = element_markdown(size = 9),
    strip.background = element_blank(),
    legend.key.spacing.y = unit(0.1, "mm"),
    legend.key.size = unit(4, "mm"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.background = element_rect(linewidth = 0.2, colour = "black"),
    legend.box.spacing = unit(1, "mm"),
    legend.position = "bottom"
  )

# Combine plots -----------------------------------------------------------

plt_fig6 <-
  free(plt_fig6a) / wrap_plots(free(plt_fig6b), plt_fig6c, widths = c(1, 1)) +
  plot_layout(heights = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")

# Output ------------------------------------------------------------------

ggsave(plot = plt_fig6,
  here('output', 'figure6.pdf'), dpi = 600, width = 15, height = 9)
