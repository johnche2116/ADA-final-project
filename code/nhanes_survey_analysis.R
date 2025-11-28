# nhanes_survey_analysis.R
# NHANES 2019–2020 Survey Analysis: Physical Activity and BP Control
# Outputs CSV files with descriptive stats and regression results


#1 Load packages
library(tidyverse)
library(survey)
library(foreign)
library(broom)
library(readr)

#2 Load NHANES datasets (place XPT files in data/ folder)
demo <- read.xport("data/P_DEMO.xpt")
bpx  <- read.xport("data/P_BPXO.xpt")
paq  <- read.xport("data/P_PAQ.xpt")
bmx  <- read.xport("data/P_BMX.xpt")
smq  <- read.xport("data/P_SMQ.xpt")
alq  <- read.xport("data/P_ALQ.xpt")
diq  <- read.xport("data/P_DIQ.xpt")

# Merge datasets by SEQN
nhanes <- demo %>%
  left_join(bpx, by = "SEQN") %>%
  left_join(paq, by = "SEQN") %>%
  left_join(bmx, by = "SEQN") %>%
  left_join(smq, by = "SEQN") %>%
  left_join(alq, by = "SEQN") %>%
  left_join(diq, by = "SEQN")

#4 Restrict to adults (≥18 years)
nhanes <- nhanes %>% filter(RIDAGEYR >= 18)

#5 Create mean SBP/DBP and BP control
nhanes <- nhanes %>%
  mutate(
    mean_sbp = rowMeans(select(., BPXOSY1, BPXOSY2, BPXOSY3), na.rm = TRUE),
    mean_dbp = rowMeans(select(., BPXODI1, BPXODI2, BPXODI3), na.rm = TRUE),
    bp_control = ifelse(mean_sbp < 130 & mean_dbp < 80, 1, 0)
  )

#6 Create Physical Activity variables
nhanes <- nhanes %>%
  mutate(
    mod_min = ifelse(is.na(PAD615), 0, PAD615),
    vig_min = ifelse(is.na(PAD630), 0, PAD630),
    mvpa_min = mod_min + (2 * vig_min),
    meet_pa = ifelse(mvpa_min >= 150, 1, 0)
  )

#7 Survey Design
nhanes_design <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTMECPRP,
  nest = TRUE,
  data = nhanes
)

#8 Weighted Descriptive Statistics
# Continuous variables
weighted_means <- svymean(~mean_sbp + mean_dbp + RIDAGEYR + BMXBMI, nhanes_design, na.rm = TRUE)
weighted_means_df <- as.data.frame(weighted_means)
names(weighted_means_df) <- c("Mean", "SE")
weighted_means_df$Variable <- rownames(weighted_means_df)
weighted_means_df <- weighted_means_df[, c("Variable", "Mean", "SE")]

write.csv(weighted_means_df, "outputs/Weighted_Means.csv", row.names = FALSE)
write.csv(weighted_means_df, "outputs/Weighted_Means.csv", row.names = FALSE)


# Categorical: PA vs BP control
weighted_table <- svytable(~meet_pa + bp_control, nhanes_design)
weighted_table_df <- as.data.frame(weighted_table)
names(weighted_table_df) <- c("meet_pa", "bp_control", "Weighted_Count")

write.csv(weighted_table_df, "outputs/PA_vs_BP_Crosstab.csv", row.names = FALSE)

#9 Weighted Poisson Regression (Main Model)
model_pr <- svyglm(
  bp_control ~ meet_pa + RIDAGEYR + factor(RIAGENDR) + factor(RIDRETH1) +
    BMXBMI + SMQ020 + ALQ111 + DIQ010,
  design = nhanes_design,
  family = quasipoisson(link = "log")
)

pr_estimates <- exp(cbind(PR = coef(model_pr), confint(model_pr)))
write.csv(as.data.frame(pr_estimates), "outputs/BP_Control_MainModel_PR.csv", row.names = TRUE)

#10 Interaction Model (PA × sex, PA × BMI)
model_inter <- svyglm(
  bp_control ~ meet_pa * factor(RIAGENDR) + meet_pa * BMXBMI + RIDAGEYR +
    factor(RIDRETH1) + SMQ020 + ALQ111 + DIQ010,
  design = nhanes_design,
  family = quasipoisson(link = "log")
)

pr_inter_estimates <- exp(cbind(PR = coef(model_inter), confint(model_inter)))
write.csv(as.data.frame(pr_inter_estimates), "outputs/BP_Control_InteractionModel_PR.csv", row.names = TRUE)

getwd()

if(!dir.exists("outputs")) dir.create("outputs")
write.csv(weighted_means_df, "outputs/Weighted_Means.csv", row.names = FALSE)
write.csv(weighted_table_df, "outputs/PA_vs_BP_Crosstab.csv", row.names = FALSE)
write.csv(as.data.frame(pr_estimates), "outputs/BP_Control_MainModel_PR.csv", row.names = TRUE)
write.csv(as.data.frame(pr_inter_estimates), "outputs/BP_Control_InteractionModel_PR.csv", row.names = TRUE)
C:/Users/ryzen7/Documents/ADA-final-project/outputs/
