#heritability estimate
#load packages
options(brms.backend = "cmdstanr")
library(brms)
library(tidybayes)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(nadiv)
library(sjPlot)
library(cmdstanr)
library(tidyverse)
library(data.table)
library(Matrix)
library(MASS)
library(genio)
library(lqmm)
library(ggdist)
library(lme4)

# load data
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")
PrunedPed_bodymass <- read.csv("./datasets/PrunedPed_bodymass_20260503.csv")
Amat <- as.matrix(nadiv::makeA(PrunedPed_bodymass))

df_clean$minutes_s <- scale(df_clean$minutes)
df_clean$animal <- df_clean$BirdID

df_clean_age <- df_clean %>%
  mutate(AgeClass = case_when(age_year < 1 ~ "Juv",
                              age_year >= 1 & age_year <= 5 ~ "Young",
                              age_year > 5 ~ "Old"))
df_clean_age$AgeClass <- factor(df_clean_age$AgeClass, 
                                levels = c("Juv", "Young", "Old"))
get_prior(bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
               minutes_s + avg_invert + group_size +
               (0 + AgeClass | gr(animal, cov = Amat)) + (0 + AgeClass | BirdID) + (0 + AgeClass | BirthYear) +
               (0 + AgeClass | CatchYear) + (0 + AgeClass | Observer),
             sigma ~ AgeClass - 1),
          data    = df_clean_age,
          data2   = list(Amat = Amat),
          family  = gaussian())
# prior
my_priors <- c(
  prior(normal(0, 1),   class = b, coef = age_year),
  prior(normal(0, 0.5), class = b, coef = Iage_yearE2),
  prior(normal(0, 2),   class = b, coef = RightTarsus),
  prior(normal(0, 1),   class = b, coef = SexEstimate),
  prior(normal(0, 1),   class = b, coef = summer),
  prior(normal(0, 1),   class = b, coef = minutes_s),
  prior(normal(0, 1),   class = b, coef = avg_invert),
  prior(normal(0, 1),   class = b, coef = group_size),
  prior(normal(15.5, 3), class = Intercept),
  prior(student_t(3, 0, 1), class = sd),
  prior(lkj(2), class = cor),
  prior(normal(log(2), 0.5), class = b, coef = AgeClassJuv,   dpar = "sigma"),
  prior(normal(log(2), 0.5), class = b, coef = AgeClassYoung, dpar = "sigma"),
  prior(normal(log(2), 0.5), class = b, coef = AgeClassOld,   dpar = "sigma"))

m_h2_age <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
       minutes_s + avg_invert + group_size +
       (0 + AgeClass | gr(animal, cov = Amat)) + (0 + AgeClass | BirdID) + (0 + AgeClass | BirthYear) +
       (0 + AgeClass | CatchYear) + (0 + AgeClass | Observer),
     sigma ~ AgeClass - 1),
  data    = df_clean_age,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = my_priors,
  warmup  = 3000,
  iter    = 15000,
  chains  = 4,
  cores   = 4,
  threads = threading(4),
  seed    = 123,
  file    = "m_brms_h2_age20260507",
  backend = "cmdstanr")

#estimate h2 from genomic data
grm.data <- read_grm("./datasets/GRM/grm_output_gwas_2.9m_20260506.grm")

fam <- read.table("./datasets/GRM/filterd.final.autosome.maf0.05_bodymass_h2_gwas.fam",
                  header = FALSE)
# GRM
G1a <- grm.data$kinship

rownames(G1a) <- fam[,2]
colnames(G1a) <- fam[,2]

model_birds <- unique(df_clean$animal)

cat("Birds in model:", length(model_birds), "\n")
cat("Birds in GRM:", nrow(G1a), "\n")
cat("Overlap:", sum(model_birds %in% rownames(G1a)), "\n")
cat("Missing from GRM:", sum(!model_birds %in% rownames(G1a)), "\n")
cat("Missing from model:", sum(!rownames(G1a) %in% model_birds), "\n")

# Keep only birds present in model and GRM
birds_in_both <- model_birds[model_birds %in% rownames(G1a)]
cat("Birds to use in SNP model:", length(birds_in_both), "\n")

GRM_sub <- G1a[as.character(birds_in_both), as.character(birds_in_both)]
dim(GRM_sub)
# make positive definite
G1a_sub_final <- make.positive.definite(GRM_sub)

#subset data
df_snp <- df_clean_age[df_clean_age$animal %in% birds_in_both, ]

m_h2_age_genomic <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
       minutes_s + avg_invert + group_size +
       (0 + AgeClass | gr(animal, cov = G1a_sub_final)) + (0 + AgeClass | BirdID) + (0 + AgeClass | BirthYear) +
       (0 + AgeClass | CatchYear) + (0 + AgeClass | Observer),
     sigma ~ AgeClass - 1),
  data    = df_snp,
  data2   = list(G1a_sub_final = G1a_sub_final),
  family  = gaussian(),
  prior   = my_priors,
  warmup  = 3000,
  iter    = 15000,
  chains  = 4,
  cores   = 4,
  threads = threading(16),
  seed    = 123,
  file    = "m_brms_h2_age_genomic20260507",
  backend = "cmdstanr")

