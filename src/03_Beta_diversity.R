# Beta Diversity SynCom20 analysis on Arabidopsis thaliana under herbivory
# Data analysis
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(vegan)
library(pairwiseAdonis)

# Global Options ----------------------------------------------------------
set.seed(123)

# Data Input --------------------------------------------------------------

otu <- 
  readRDS(here("input", "processed", "otu_df.rds"))

metadata <- 
  read.csv(here("input", "metadata.csv"))
  

# Filter conditions -------------------------------------------------------
# Subset larva fed on SynCom20-leaves, SynCom20-treated control leaves
# and control fed leaves

subset_conditions <- c("larvae_SynCom_Feeding", 
                       "leaf_SynCom_noFeeding", 
                       "leaf_SynCom_Feeding")

# Create a data frame with the defined subset of conditions
meta.df <-
  metadata |> 
  filter(Condition %in% subset_conditions) |> 
  column_to_rownames(var = "SampleName")

samples <- rownames(meta.df)


# SynCom20 OTU table ------------------------------------------------------
# Subset
comm.otu <- 
  otu |> 
  select(species, all_of(samples)) |> 
  group_by(species) |> 
  summarise(across(everything(), sum), .groups = "drop") |> 
  column_to_rownames(var = "species")

# Relative abundances
comm.otu.rel <- decostand(t(comm.otu), method = "total")

# Bray-Curtis dissimilarity matrix
bray_dist <- vegdist(comm.otu.rel, method = "bray")


# NMDS --------------------------------------------------------------------

nmds <- metaMDS(bray_dist, k = 2, trymax = 100)

taxa_fit <- envfit(nmds, comm.otu.rel, permutations = 999)

# Merge NMDS results with metadata
nmds_df <- 
  as.data.frame(scores(nmds, "sites")) |> 
  rownames_to_column(var = "SampleName") |> 
  left_join(metadata, by = "SampleName")

fit_df <-
  as.data.frame(scores(taxa_fit, "vectors") * ordiArrowMul(taxa_fit)) |> 
  mutate(
    species = rownames(taxa_fit$vectors$arrows),
    R2 = taxa_fit$vectors$r,
    p = taxa_fit$vectors$pvals) |> 
  left_join(otu |> select(phylum:species), by = "species") |> 
  filter(p < 0.05)

# PERMANOVA using adonis2
adonis_result <- adonis2(bray_dist ~ Condition, data = meta.df)
print(adonis_result)

# Save Files --------------------------------------------------------------

saveRDS(nmds_df, here("input", "processed", "nmds.rds"))
saveRDS(fit_df, here("input", "processed", "fit.nmds.rds"))
