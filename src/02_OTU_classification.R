# OTU Classification
# Author: Rudolf Schlechter

# Load libraries ----------------------------------------------------------

library(here)
library(tidyverse)
library(DECIPHER)
library(dada2)
library(phyloseq)

# Input Data --------------------------------------------------------------

# Sequence table
comm_matrix <- 
  readRDS(here("input", "processed", "otu_others.rds")) |> 
  column_to_rownames(var = "OTU")

seq <- readDNAStringSet(here("input", "otu_notclassified.fa"))

# Classification ----------------------------------------------------------
taxa <- 
  assignTaxonomy(
    seq,
    here("input", "silva_nr99_v138.2_toGenus_trainset.fa"),
    multithread = TRUE)

rownames(taxa) <- seq@ranges@NAMES

# Analysis ----------------------------------------------------------------

ps <- 
  phyloseq(
    otu_table(comm_matrix, taxa_are_rows=TRUE), 
    tax_table(taxa))

# Output ------------------------------------------------------------------

saveRDS(taxa, here("input", "processed", "taxa_silva.rds"))

saveRDS(ps, here("input", "processed", "ps_others.rds"))
