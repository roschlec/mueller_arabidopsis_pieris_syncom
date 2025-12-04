# Relative abundance analysis ---------------------------------------------
# Data processing
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)

# Data Input --------------------------------------------------------------
# Read sequencing data set analysed with pairwise alignment based on reference sequences
syncom20.df <- 
  readRDS(here("input", "processed", "otu_df.rds"))

# Relative abundance data set
df_syncom20 <- 
  syncom20.df |> 
  pivot_longer(cols = -phylum:-species, names_to = "sample") |> 
  group_by(genus, species, sample) |> 
  summarise(count = sum(value), .groups = "drop") |>
  unite(genus_species, genus:species, sep = " ") |> 
  separate(sample, into = c("condition", "rep"), sep = "_", remove = FALSE) |> 
  group_by(sample) |> 
  mutate(freq = count / sum(count),
        # define treatment as SampleType_SynCom_Larvae
         treatment = case_when(
           # Leaf samples
           condition == "CL"  ~ "leaf_0_0",
           condition == "CF"  ~ "leaf_0_1",
           condition == "SNF" ~ "leaf_1_0",
           condition == "SF"  ~ "leaf_1_1",
           # Larval samples
           condition == "LC"  ~ "larvae_0_1",
           condition == "LS"  ~ "larvae_1_NA",
           # Controls
           condition == "EC"  ~ "mock_0_0",
           condition == "IN"  ~ "inoculum_NA_NA"))

## Remove controls
control <- c("CL", "EC", "LC", "CF")

##  Final data set
df_syncom20_trt <-
  df_syncom20 |> 
  filter(!condition %in% control) |> 
  separate(treatment, 
           into = c("sampleType", "SynCom", "Feeding"), sep = "_", remove = FALSE) |> 
  mutate(across(SynCom:Feeding, ~case_when(. == 1 ~ "+", . == 0 ~ "-", TRUE ~ ""))) |> 
  mutate(
    sampleType = factor(sampleType, levels = c("inoculum", "leaf", "larvae"),
                        labels = c("Inoculum", "Leaf", "Larvae")),
    SynCom = factor(SynCom, levels = c("", "-", "+")),
    Feeding = factor(Feeding, levels = c("", "-", "+"))
  )
