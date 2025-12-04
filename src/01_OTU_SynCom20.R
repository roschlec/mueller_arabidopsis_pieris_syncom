# OTU Classification
# Data processing
# Author: Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)

# Data input --------------------------------------------------------------

otu_classification <- 
  read_delim(here('input', 'otu_syncom20.txt'), 
             col_names = FALSE, show_col_types = FALSE) |> 
  select(OTU = X1, X2) |> 
  separate(X2, 
           into = c("root", "domain", "phylum", "class", "order", "family", "genus", "species"),
           sep = ";")

otu_table <-
  read_delim(here('input', 'otu_table.txt'), 
             delim = "\t", show_col_types = FALSE)


# OTU associated with SynCom20 members ------------------------------------

OTU_SynCom20 <- otu_classification |> pull(OTU)


# OTU not classified as SynCom20 members ----------------------------------

otu_others <-
  otu_table |> 
  filter(!OTU %in% OTU_SynCom20)


# OTU Table with SynCom20 classification ----------------------------------

otu_df <-
  left_join(otu_table, otu_classif, by = "OTU") |> 
  select(-OTU, -root, -domain) |> 
  group_by(phylum, class, order, family, genus, species) |> 
  summarise(across(everything(), sum), .groups = "drop") |> 
  na.omit()


# Output ------------------------------------------------------------------

##  Save data frame

write_rds(otu_df, here("input", "processed", "otu_df.rds"))

write_rds(otu_others, here("input", "processed", "otu_others.rds"))

