# Differential abundance SynCom20 members
# Data analysis
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
library(forcats)
library(DESeq2)
library(phyloseq)
library(apeglm)

# Data Input --------------------------------------------------------------
otu <- readRDS(here("input", "processed", "otu_df.rds"))
meta <- read.csv(here("input", "metadata.csv"))


# SynCom20 taxa -----------------------------------------------------------

taxa <-
  otu |> 
  select(phylum:species) |> 
  unique() |>
  as_tibble() |> 
  arrange(species) |> 
  column_to_rownames(var = "species") |> 
  as.matrix()
colnames(taxa) <- c("Phylum", "Class", "Order", "Family", "Genus")
  

# OTU table ---------------------------------------------------------------

comm.otu <- 
  otu |> 
  select(-phylum:-genus) |> 
  group_by(species) |> 
  summarise(across(everything(), sum), .groups = "drop") |> 
  column_to_rownames(var = "species")

meta <-
  meta |> 
  column_to_rownames(var = "SampleName")

otu_element <- otu_table(comm.otu, taxa_are_rows = TRUE)
taxa_element <- tax_table(taxa)

physeq <- 
  phyloseq(
    otu_element, 
    taxa_element,
    sample_data(meta))

physeq <- prune_samples(sample_sums(physeq) > 0, physeq)
sample_data(physeq)


# Differential abundance with DESeq2 --------------------------------------

# Create DESeq2 object
dds <- phyloseq_to_deseq2(physeq, ~ Condition)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds)

# Define contrasts
contrast_feeding <- c("Condition", "leaf_SynCom_Feeding", "leaf_SynCom_noFeeding")
contrast_larva_fedLeaf <- c("Condition", "larvae_SynCom_Feeding", "leaf_SynCom_Feeding")
contrast_larva_nonfedLeaf <- c("Condition", "larvae_SynCom_Feeding", "leaf_SynCom_noFeeding")

list_contrasts <- list(contrast_feeding, 
                       contrast_larva_fedLeaf, 
                       contrast_larva_nonfedLeaf)

# Obtain results of contrasts and log FC shrinkage
res_deseq <- 
  lapply(
    list_contrasts, function(con) {
      results(dds, contrast = con, alpha = 0.05)})
names(res_deseq) <- c("contrast_feeding", "contrast_larva_fedLeaf", "contrast_larva_nonfedLeaf")

lres_deseq <-
  lapply(
    list_contrasts, function(con){
      lfcShrink(dds, contrast = con, res = results(dds, contrast = con), type = "normal")
    }
  )
names(lres_deseq) <- c("contrast_feeding", "contrast_larva_fedLeaf", "contrast_larva_nonfedLeaf")

# Combine results
res.df <- 
  do.call(
    rbind, lapply(names(lres_deseq), function(name) {
      df <- lres_deseq[[name]]
      df$Contrast <- name
      df$species <- rownames(df)
      return(df)
      })) |> 
  as.data.frame() |> 
  mutate(signif = padj < 0.05,
         log_padj = log(padj)) |>
  left_join(otu |> select(phylum:species), by = "species", relationship = "many-to-many") |> 
  mutate(genus_species = paste(genus, species))

# Save File ---------------------------------------------------------------
saveRDS(res.df, here("input", "processed", "diffential.abundance.rds"))

