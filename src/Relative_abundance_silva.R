# Relative abundance everything -------------------------------------------
# Plotting
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(phyloseq)
library(forcats)
library(cowplot)
library(patchwork)

# Data Input --------------------------------------------------------------

trt <- c("leaf_noSynCom_noFeeding", 
         "leaf_noSynCom_Feeding", 
         "leaf_SynCom_noFeeding",
         "leaf_SynCom_Feeding",
         "larvae_noSynCom_Feeding",
         "larvae_SynCom_Feeding",
         "inoculum",
         "extraction")

trt_label <- c("leaf_0_0", 
               "leaf_0_1",
               "leaf_1_0", 
               "leaf_1_1",
               "larvae_0_NA",
               "larvae_1_NA",
               "inoculum_NA_NA", 
               "kitome_NA_NA")

# Metadata
meta <- 
  read.csv(here("input", "metadata.csv")) |> 
  dplyr::rename(Sample = SampleName)

# SynCom0 reads
syncom20 <- read_rds(here("input", "processed", "otu_df.rds"))

# Not classified OTU
others <- read_rds(here("input", "processed", "ps_others.rds"))

# Rename and create OTU data frame
other_df <-
  others |> 
  psmelt() |> 
  select(phylum = Phylum, 
         class = Class,
         order = Order,
         family = Family,
         genus = Genus,
         species = OTU,
         Sample, Abundance) 

# OTU from mock samples (kitome)
mock_otu <-
  other_df |> 
  filter(!order %in% c("Chloroplast", "Rickettsiales")) |> 
  filter(grepl("Mock|^C[0-9]", Sample) & Abundance > 0) |>
  pull(species)

# Collate non-classified reads into Chloroplast, Mitochondria and rest
relabel_df <-
  other_df |> 
  filter(!species %in% mock_otu) |> 
  mutate(label = case_when(
    order == "Chloroplast" ~ "Cloroplast",
    family == "Mitochondria" ~ "Mitochondria",
    TRUE ~ "Other")) |> 
  group_by(label, Sample) |> 
  summarise(Abundance = sum(Abundance), .groups = "drop")


# Final Data Frames -------------------------------------------------------

# Total Reads
full_df <-
  syncom20 |> 
  pivot_longer(-phylum:-species, names_to = "Sample", values_to = "Abundance") |> 
  select(-species) |> 
  mutate(label = "SynCom20") |> 
  group_by(label, Sample) |> 
  summarise(Abundance = sum(Abundance), .groups = "drop") |> 
  rbind(relabel_df) |> 
  left_join(meta, by = "Sample") |> 
  mutate(Condition = factor(Condition, levels = trt, labels = trt_label)) |> 
  separate(Condition, into = c("sampleType", "SynCom", "Feeding"), sep = "_", remove = FALSE) |> 
  mutate(across(SynCom:Feeding, ~case_when(. == 1 ~ "+", . == 0 ~ "-", TRUE ~ ""))) |> 
  mutate(
    sampleType = factor(sampleType, levels = c("kitome", "inoculum", "leaf", "larvae"),
                        labels = c("Kitome", "Inoculum", "Leaf", "Larvae")),
    SynCom = factor(SynCom, levels = c("", "-", "+")),
    Feeding = factor(Feeding, levels = c("", "-", "+"))
  )

# Relative abundance
otu_rel <-
  full_df |> 
  pivot_wider(id_cols = label, names_from = "Sample", 
              values_from = "Abundance", values_fill = 0) |> 
  mutate(across(-label, ~ 100* .x/sum(.x))) |> 
  mutate(across(-label, ~ ifelse(is.nan(.x), 0, .x))) |> 
  pivot_longer(cols = -label, names_to = "Sample", values_to = "Abundance") |> 
  left_join(meta, by = "Sample") |> 
  mutate(Condition = factor(Condition, levels = trt, labels = trt_label)) |> 
  separate(Condition, into = c("sampleType", "SynCom", "Feeding"), sep = "_", remove = FALSE) |> 
  mutate(across(SynCom:Feeding, ~case_when(. == 1 ~ "+", . == 0 ~ "-", TRUE ~ ""))) |> 
  mutate(
    sampleType = factor(sampleType, levels = c("kitome", "inoculum", "leaf", "larvae"),
                        labels = c("Kitome", "Inoculum", "Leaf", "Larvae")),
    SynCom = factor(SynCom, levels = c("", "-", "+")),
    Feeding = factor(Feeding, levels = c("", "-", "+"))
  )

# Relative abundance without Chloroplasts and Mitochondria
otu_rel_noMit_noChl <-
  full_df |> 
  mutate(Abundance = ifelse(grepl("Cloroplast|Mitochondria", label), 0, Abundance)) |> 
  pivot_wider(id_cols = label, names_from = "Sample", 
              values_from = "Abundance", values_fill = 0) |> 
  mutate(across(-label, ~ 100* .x/sum(.x))) |> 
  mutate(across(-label, ~ ifelse(is.nan(.x), 0, .x))) |> 
  pivot_longer(cols = -label, names_to = "Sample", values_to = "Abundance") |> 
  left_join(meta, by = "Sample") |> 
  mutate(Condition = factor(Condition, levels = trt, labels = trt_label)) |> 
  separate(Condition, into = c("sampleType", "SynCom", "Feeding"), sep = "_", remove = FALSE) |> 
  mutate(across(SynCom:Feeding, ~case_when(. == 1 ~ "+", . == 0 ~ "-", TRUE ~ ""))) |> 
  mutate(
    sampleType = factor(sampleType, levels = c("kitome", "inoculum", "leaf", "larvae"),
                        labels = c("Kitome", "Inoculum", "Leaf", "Larvae")),
    SynCom = factor(SynCom, levels = c("", "-", "+")),
    Feeding = factor(Feeding, levels = c("", "-", "+"))
  )


# Plot --------------------------------------------------------------------

# Total Reads
plt_abundance <-
  full_df |> 
  ggplot(aes(x = Sample, y = Abundance, fill = label)) +
  facet_wrap(~ sampleType + SynCom + Feeding, scales = "free_x", ncol = 8) + 
  geom_bar(stat = "identity", color = "black") +
  labs(y = "Total reads")

# Relative abundance
plt_relabundance <-
  otu_rel |> 
  ggplot(aes(x = Sample, y = Abundance, fill = label)) +
  facet_wrap(~ sampleType + SynCom + Feeding, scales = "free_x", ncol = 8) + 
  geom_bar(stat = "identity", color = "black") +
  labs(y = "Relative abundance [%]")

# Relative abundance without Chloroplasts and Mitochondria
plt_relabundance_flt <-
  otu_rel_noMit_noChl |>
  ggplot(aes(x = Sample, y = Abundance, fill = label)) +
  facet_wrap(~ sampleType + SynCom + Feeding, scales = "free_x", ncol = 8) + 
  geom_bar(stat = "identity", color = "black") +
  labs(y = "Relative abundance [%]")

# Combine plots -----------------------------------------------------------

plt <- 
  plt_abundance/plt_relabundance/plt_relabundance_flt +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") & 
  theme_bw() &
  scale_fill_manual(values = c("darkolivegreen", "coral1", "grey", "cyan4")) &
  guides(fill = guide_legend(title = "Taxa", ncol = 1,
                             override.aes = list(stroke = 0))) &
  labs(x = "Replicate") &
  theme(
    text = element_text(size = 12, color = "black"),
    axis.text = element_text(size = 12, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_line(colour = "black", linewidth = 0.5),
    axis.ticks.length.y = unit(1, "mm"),
    legend.key.spacing.y = unit(0.1, "mm"),
    legend.key.size = unit(3, "mm"),
    legend.title = element_text(size = 10),
    legend.background = element_rect(linewidth = 0.2, colour = "black"),
    legend.box.spacing = unit(1, "mm"),
    panel.border = element_rect(linewidth = 0.8, color = "black"),
    strip.text = element_text(size = 12),
    strip.background = element_rect(linewidth = 0.8, color = "black"))

x_lab <- 0.915
y_A <- 0.96
y_B <- 0.63
y_C <- 0.3
y_dec <- 0.025

plt_figureS2 <- 
  ggdraw(plt) +
    # PLOT A
  draw_label("Sample Type", x = x_lab, y = y_A, size = 10, hjust = 0) +
  draw_label("SynCom", x = x_lab, y = y_A - y_dec, size = 10, hjust = 0) +
  draw_label("Feeding", x = x_lab, y = y_A - y_dec*2, size = 10, hjust = 0) +
    # PLOT B
  draw_label("Sample Type", x = x_lab, y = y_B, size = 10, hjust = 0) +
  draw_label("SynCom", x = x_lab, y = y_B - y_dec, size = 10, hjust = 0) +
  draw_label("Feeding", x = x_lab, y = y_B - y_dec*2, size = 10, hjust = 0) +
    # PLOT C
  draw_label("Sample Type", x = x_lab, y = y_C, size = 10, hjust = 0) +
  draw_label("SynCom", x = x_lab, y = y_C - y_dec, size = 10, hjust = 0) +
  draw_label("Feeding", x = x_lab, y = y_C - y_dec *2, size = 10, hjust = 0)

# Save --------------------------------------------------------------------

mapply(function(x) 
  ggsave(x, 
         plot = plt_figureS2, 
         dpi = 300, height = 11, width = 15),
  x = c(here("output", "FigureS2.png"), here("output", "FigureS2.eps")))
