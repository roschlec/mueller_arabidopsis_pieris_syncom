# Effect of SynCom richness to larval performance
# Data processing and statistical analyses
# Authors: Moritz Müller and Rudolf Schlechter

# Libraries ---------------------------------------------------------------
library(here)
library(tidyverse)
# Stats
library(nlme)
library(glmmTMB)
library(car)
library(rstatix)
library(emmeans)

# Data Input --------------------------------------------------------------

weight <-
  read_delim(here("input", "larval_weight_cfu.txt"), 
             delim = "\t", show_col_types = FALSE) |> 
  mutate(treatment = as.factor(treatment),
         experiment_rep = as.factor(experiment_rep))

# Larval survival ---------------------------------------------------------

larvae_survival <-
  weight |> 
  filter(Nlarvae <= 3) |> 
  group_by(date, treatment, experiment_rep) |> 
  mutate(
    initial = 3,
    Nfraction = Nlarvae/max(Nlarvae),
    Npercent = Nfraction * 100) |> 
  ungroup()

larvae_survival |>
  group_by(treatment) |>
  summarise(mean = mean(Npercent))

# Regression Model
glm_survival <- glmmTMB(cbind(Nlarvae, initial - Nlarvae) ~ treatment + (1|experiment_rep),
                        family = binomial,
                        data = larvae_survival)
# Model evaluation (ANOVA)
Anova(glm_survival)


# CFU on Leaves -----------------------------------------------------------
# Shapiro-Wilk Normality Test 
weight |> 
  group_by(experiment_rep) |> 
  shapiro_test(logleafcfu)

# Levene Test for Homogeneity of Variances
weight |>
  group_by(experiment_rep) |> 
  levene_test(logleafcfu ~ treatment)

# Regression Model
# Linear mixed-effect model
LMM_leafcfu_trt <- 
  lme(logleafcfu ~ treatment,
      random = ~ 1 | experiment_rep,
      data = weight |> filter(logleafcfu > 0))

# Model Residuals
resid(LMM_leafcfu_trt) |> qqnorm()

# Model Evaluation (ANOVA)
anova(LMM_leafcfu_trt)

# Marginal means (per treatment)
emm_leafcfu_trt <- emmeans(LMM_leafcfu_trt, ~ treatment)

# Pairwise comparisons with BH adjustment
pairs(emm_leafcfu_trt, adjust = "BH")

# Compact letter display
cld_leafcfu <- 
  multcomp::cld(emm_leafcfu_trt, adjust = "BH", Letters = letters, alpha = 0.05) |> 
  mutate(.group = str_remove_all(.group, " "))


# CFU in Larvae -----------------------------------------------------------
# Regression Model
# Linear mixed-effect model
model_lvcfu_trt <- 
  lme(loglarvaecfu ~ treatment,
      random = ~ 1 | experiment_rep,
      data = weight)

# Model Residuals
resid(model_lvcfu_trt) |> qqnorm()

# Model Evaluation (ANOVA)
anova(model_lvcfu_trt)

# Marginal means (per treatment)
emm_lvcfu_trt <- emmeans(model_lvcfu_trt, ~ treatment)

# Pairwise comparisons with BH adjustment
pairs(emm_lvcfu_trt, adjust = "BH")

# Compact letter display
cld_lvcfu <- 
  multcomp::cld(emm_lvcfu_trt, adjust = "BH", Letters = letters, alpha = 0.05) |> 
  mutate(.group = str_remove_all(.group, " "))

# Larval Weight -----------------------------------------------------------
# Model
model_lvweight_trt <- 
  lme(weightperlarva ~ treatment,
      random = ~ 1 | experiment_rep,
      data = weight)

# Model Residuals
resid(model_lvweight_trt) |> qqnorm()

# Model Evaluation (ANOVA)
anova(model_lvweight_trt)

# Marginal means (per treatment)
emm_lvweight_trt <- emmeans(model_lvweight_trt, ~ treatment)

# Pairwise comparisons with BH adjustment
pairs(emm_lvweight_trt, adjust = "BH")

# Compact letter display
cld_lv_weight <- 
  multcomp::cld(emm_lvweight_trt, adjust = "BH", Letters = letters, alpha = 0.05,
                decreasing = TRUE) |> 
  mutate(.group = str_remove_all(.group, " "))


# Detailed stats ----------------------------------------------------------

# Larval weight -----------------------------------------------------------
### model definition
model_lm_trt <- lm(weightperlarva ~ treatment, data = weight)
model_lm_day <- lm(weightperlarva ~ experiment_rep, data = weight)
model_lm_trt_day <- lm(weightperlarva ~ treatment*experiment_rep, data = weight)
model_lme_trt <- lme(weightperlarva ~ treatment, random = ~ 1 | experiment_rep, data = weight)

###model definition normalized weight (accuracy of normalized data in representing model with y axis correction)
normmodel_lm_trt <- lm(normweight ~ treatment, data = weight)
normmodel_lm_day <- lm(normweight ~ experiment_rep, data = weight)
normmodel_lm_trt_day <- lm(normweight ~ treatment*experiment_rep, data = weight)
normmodel_lme_trt <- lme(normweight ~ treatment, random = ~ 1 | experiment_rep, data = weight)

### check residuals
plot(resid(model_lm_trt))
plot(resid(model_lm_day))
plot(resid(model_lm_trt_day))
plot(resid(model_lme_trt))

### normalised weight
plot(resid(normmodel_lm_trt))
plot(resid(normmodel_lm_day))
plot(resid(normmodel_lm_trt_day))
plot(resid(normmodel_lme_trt))

### anova
anova(model_lm_trt)
anova(model_lm_day)
anova(model_lm_trt_day) #-> As both treatment and sampling day have an effect on larval weight the lme with sampling day as a random effect is used
anova(model_lme_trt) #important

### anova normweight
anova(normmodel_lm_trt)
anova(normmodel_lm_day)
anova(normmodel_lm_trt_day) #-> sampling day (experimental run) does not effect the larval weight in the data normalized by each experimental runs average control weight -> norm data used as visual representation of model
anova(normmodel_lme_trt) 

### summary
summary(model_lm_trt)
summary(model_lm_day)
summary(model_lm_trt_day)
summary(model_lme_trt)

### group comparison
compare_lme <- emmeans(model_lme_trt, specs = ~ treatment)
compare_lme |>  as.data.frame()

contrast(compare_lme, method = 'pairwise', adjust = 'BH')


# Supplemental ------------------------------------------------------------

# Correlation between CFU on leaves or in larvae and larval weight
# Subset data

weight_onlySynCom <- weight |> filter(treatment != "control")