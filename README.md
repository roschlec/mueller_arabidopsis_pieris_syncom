<p align="left">
<!-- Zenodo DOI badge -->
<a href="https://10.5281/zenodo.17822324"> <img src="https://img.shields.io/badge/DOI-pending-blue.svg" alt="DOI"/> </a>
<!-- License badge -->
<a href="https://github.com/USERNAME/REPO/blob/main/LICENSE"> <img src="https://img.shields.io/badge/license-MIT-green.svg" alt="License: MIT"/> </a>
<!-- Built with R -->
<img src="https://img.shields.io/badge/R-%3E%3D%204.3.0-blue.svg" alt="R &gt;=4.3"/>
</p>

# From the leaf to the gut and back again: the fate and influence of phyllosphere bacteria in a gnotobiotic *Arabidopsis* - *Pieris brassicae* system

Analysis code and documentation for the manuscript on leaf-associated synthetic microbial communities shaping *Pieris brassicae* herbivory outcomes.

## Overview

This repository includes all the R scripts to reproduce the analysis of the manuscript XX published in. The scripts include:

-   Effect of SynCom richness on larval performance
-   Phytohormone analysis
-   Bacterial community analysis

The study investigates how defined synthetic leaf microbial communities influence herbivore feeding outcomes and how herbivory reshapes leaf-associated microbial communities.

All datasets to reproduce the findings of our publication are stored in the associated <a href="https://10.5281/zenodo.17822324">Zenodo</a> repository. Raw 16S rRNA gene amplicon sequencing data are available in the EMBL-ENI ENA BioProject accession XXX.

## Repository Structure

```         
📂 data_input/                # Must be created locally and include the datasets from the Zenodo repository
 📂 processed/               # Must be created locally
📂 output/                    # Must be created locally to store output objects (plots)
📂 src/                       # Core analysis scripts and helpers
 📄 SynCom_richness.R        # Effect of SynCom richness on larval performance
 📄 SynCom_richness_plot.R   # Figure 3
 📄 SynCom_composition.R     # Figure S2 Community structure summaries
 📄 Phytohormones.R          # Analysis of phytohormones quantification
 📄 Phytohormones_plot.R     # Figure 4
 📄 SynCom20effect.R         # Figure 5 SynCom20 effect on Pieris
 📄 01_OTU_SynCom20.R        # OTU table processing for SynCom20
 📄 02_OTU_classification.R  # Taxonomy assignment workflows
 📄 03_Beta_diversity.R      # Beta diversity metrics and ordination
 📄 04_Differential_abundance.R  # Differential abundance of SynCom20 members
 📄 05_RelativeAbundance.R   # Abundance aggregation and genus-level summaries
 📄 06_Community_plot.R.     # Figure 6 Final community-level visualisation
 📄 Alpha_diversity.R        # Figure 7 Richness/diversity assessment
 📄 Relative_abundance_silva.R   # SILVA taxonomy-based relative abundance summarisation
 📄 Larval_Survival.R        # Larval survival in sterile and open growth systems
```

## Workflow

1.  Clone this repository

```         
git clone https://github.com/relab-fuberlin/mueller_arabidopsis_pieris_syncom.git
cd mueller_arabidopsis_pieris_syncom
```

2.  Create an input and output directory

```         
mkdir -p data_input
mkdir -p data_input/processed
mkdir -p output
```

3.  Download the datasets and metadata available in <a href="https://10.5281/zenodo.17822324">Zenodo</a>. Download the contents into the `data_input/` directory.

4.  Open RStudio and run the workflow

Execute scripts sequentially:

```         
# Figure 3
source("src/SynCom_richness.R")
source("src/SynCom_richness_plot.R")

# Figure 4
source("src/Phytohormones.R")
source("src/Phytohormones_plot.R")

# Figure 5
source("src/SynCom20effect.R")

# Figure 6
source("src/01_OTU_SynCom20.R")
source("src/02_OTU_classification.R")
source("src/03_Beta_diversity.R")
source("src/04_Differential_abundance.R")
source("src/05_RelativeAbundance.R")
source("src/06_Community_plot.R")

# Figure 7
source("src/Alpha_diversity.R")
```
To recreate the supplementary material:
```
# Figure S1
source("src/Larval_Survival.R")

# Figure S2
source("src/SynCom_composition.R")

# Figure S3
# Generated from main Figure 3 scripts

# Figure S4
source("src/Relative_abundance_silva.R")
```

Figures and processed outputs will be saved automatically into the `output/` directory.

## Reference

If you use this repository, please cite:

> Add Publication link

> Müller, M., Huve, M. A. P., Kunzler, M., Remus-Emsermann, M., Schlechter, R., & Paniagua Voirol, L. R. (2025). Gnotobiotic Arabidopsis–Pieris brassicae–Microbiota System Dataset [Data set]. Zenodo. https://doi.org/10.5281/zenodo.17822324
