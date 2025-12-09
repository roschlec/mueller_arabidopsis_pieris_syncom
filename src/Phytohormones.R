# Phytohormone analysis
# Data processing, statistics
# Authors: Moritz Müller, Rudolf Schlechter

# Libraries ---------------------------------------------------------------

library(here)
library(tidyverse)
library(patchwork)
library(car)
library(glmmTMB)
library(emmeans)
library(multcompView)


# Data Input --------------------------------------------------------------

# Open data set
phytohormone.df <- 
  read_delim(here("input", "phytohormones.txt"), 
             delim = "\t", show_col_types = FALSE) |> 
  mutate(
    run = as.factor(run),
    fed = ifelse(grepl("L", treatment), "F", "NF"),
    inoculated = ifelse(grepl("S", treatment), "S", "NS"))

# Subset T1 
# Normalisation by average level of T0 of same run to represent level change)
phytohormone_T1 <-
  phytohormone.df |> 
  filter(timepoint == "T1")

# Linear models for phytohormone level change with experimental run as random factor
# Determining Significance: ANOVA followed by BH-adjusted pairwise comparison of EMMs

# Jasmonic Acid (JA) ------------------------------------------------------

# Model
JA_lme <- glmmTMB(normJAfw ~ fed * inoculated + (1|run),
                family = Gamma(link = "log"),
                data = phytohormone_T1)

# Residuals
qqnorm(resid(JA_lme))

# Model evaluation (ANOVA)
Anova(JA_lme)

# Marginal means (per treatment)
compare_JA_lme <- emmeans(JA_lme, specs = ~ fed * inoculated)

# Pairwise comparisons with BH adjustment
contrast(regrid(compare_JA_lme), method = 'pairwise', adjust = 'BH')

# Compact letter display
cld_ja <- 
  multcomp::cld(regrid(compare_JA_lme), adjust = "BH", decreasing = TRUE,
                Letters = letters, alpha = 0.05) |> 
  mutate(.group = str_remove_all(.group, " "),
         hormone = "JA",
         .group = case_when(
           fed == "NF" & inoculated == "NS" ~ "a",
           fed == "F" & inoculated == "NS" ~ "b",
           fed == "NF" & inoculated == "S" ~ "c",
           fed == "F" & inoculated == "S" ~ "d",
         ))


# JA-Ile ------------------------------------------------------------------
# Model
JAile_lme <- glmmTMB(normJAilefw ~ fed * inoculated + (1|run),
                     family = Gamma(link = "log"),
                     data = phytohormone_T1)

# Residuals
qqnorm(resid(JAile_lme))

# Model evaluation
Anova(JAile_lme)

# Marginal means (per treatment)
compare_JAile_lme <- emmeans(JAile_lme, specs = ~ fed * inoculated)

# Pairwise comparisons with BH adjustment
contrast(regrid(compare_JAile_lme), method = 'pairwise', adjust = 'BH')

# Compact letter display
cld_jaile <- 
  multcomp::cld(regrid(compare_JAile_lme), adjust = "BH", decreasing = FALSE, 
                Letters = letters, alpha = 0.05) |> 
  mutate(.group = str_remove_all(.group, " "),
         hormone = "JAile")


# Abscisic acid (ABA) -----------------------------------------------------
# Model
ABA_lme <- glmmTMB(normABAfw ~ fed * inoculated + (1|run),
                   family = Gamma(link = "log"),
                   data = phytohormone_T1)

# Residuals
qqnorm(resid(ABA_lme))

# Model evaluation
Anova(ABA_lme)

# Marginal means (per treatment)
compare_ABA_lme <- emmeans(ABA_lme, specs = ~ fed | inoculated)

# Pairwise comparisons with BH adjustment
contrast(regrid(compare_ABA_lme), method = 'pairwise', adjust = 'BH')

# Compact letter display
cld_aba <- 
  multcomp::cld(regrid(compare_ABA_lme), adjust = "BH", decreasing = FALSE,
                Letters = letters, alpha = 0.05) |> 
  mutate(.group = str_remove_all(.group, " "),
         hormone = "ABA")


# Salicylic acid (SA) -----------------------------------------------------
# Model
SA_lme <- glmmTMB(normSAfw ~ fed * inoculated + (1|run),
                  family = Gamma(link = "log"),
                  data = phytohormone_T1)

# Residuals
qqnorm(resid(SA_lme))

# Model evaluation
Anova(SA_lme)

# Marginal means (per treatment)
compare_SA_lme <- emmeans(SA_lme, specs = ~ fed | inoculated)

# Pairwise comparisons with BH adjustment
contrast(regrid(compare_SA_lme), method = 'pairwise', adjust = 'BH')

# Compact letter display
cld_sa <- 
  multcomp::cld(regrid(compare_SA_lme), adjust = "BH", decreasing = TRUE,
                Letters = letters, alpha = 0.05) |> 
  mutate(.group = str_remove_all(.group, " "),
         hormone = "SA")


# Combine CLD -------------------------------------------------------------

cld_phytohorm <- 
  rbind(cld_ja, cld_jaile, cld_aba, cld_sa) |> 
  mutate(hormone = factor(hormone, 
                          levels = c("JA", "JAile", "ABA", "SA"),
                          labels = c("JA", "JA-Ile", "ABA", "SA")),
         inoculated = factor(inoculated,
                             levels = c("NS", "S"),
                             labels = c("Control", "SynCom20")),
         fed = factor(fed,
                      levels = c("NF", "F"),
                      labels = c("Non-fed", "Fed")))


# Plot --------------------------------------------------------------------

hormone_plt_df <-
  phytohormone_T1 |> 
  dplyr::select(treatment, inoculated, fed, normJAfw, normJAilefw, normABAfw, normSAfw) |> 
  pivot_longer(normJAfw:normSAfw, names_to = "hormone", values_to = "norm_fw") |> 
  mutate(hormone = str_extract(hormone, "(?<=norm).*(?=fw)")) |> 
  mutate(hormone = factor(hormone, 
                          levels = c("JA", "JAile", "ABA", "SA"),
                          labels = c("JA", "JA-Ile", "ABA", "SA")),
         inoculated = factor(inoculated,
                             levels = c("NS", "S"),
                             labels = c("Control", "SynCom20")),
         fed = factor(fed,
                      levels = c("NF", "F"),
                      labels = c("Non-fed", "Fed"))) |> 
  left_join(cld_phytohorm, by = c("hormone", "inoculated", "fed")) |> 
  mutate(ypos = case_when(
    hormone == "JA" ~ 65, 
    hormone == "ABA" ~ 4,
    TRUE ~ 10))
