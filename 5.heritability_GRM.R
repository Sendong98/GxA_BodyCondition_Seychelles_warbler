#h2 estimates from GRM
#h2 from GRM
#load packages
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

# load data
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")

df_clean$minutes_s <- scale(df_clean$minutes)
df_clean$animal <- df_clean$BirdID

#run GCTA in R by using system() function or run the command line in terminal to generate GRM
#also a tutorial for this part: https://wildanimalmodels.org/docs/univariate/grm/

#system("/your/path/to/gcta64 or gcta64.exe --bfile filterd.final.autosome.maf0.05_bodymass_h2  --make-grm --autosome-num 29 --out filterd.final.autosome.maf0.05")

#load the GRM, if you want to check how to 
grm.data <- read_grm("./datasets/GRM/grm_output_gwas_2.9m_20260506.grm")

#fam contains the birdid and other information
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
df_snp <- df_clean[df_clean$animal %in% birds_in_both, ]

cat("Original observations:", nrow(df_clean), "\n")
cat("Observations after subsetting:", nrow(df_snp), "\n")

# prior
my_priors <- c(prior(normal(0, 1),class = b,coef = age_year),
               prior(normal(0, 0.5),class = b,coef = Iage_yearE2),
               prior(normal(0, 2),class = b,coef = RightTarsus),
               prior(normal(0, 1),class = b,coef = SexEstimate),
               prior(normal(0, 1),class = b,coef = summer),
               prior(normal(0, 1),class = b,coef = minutes_s),
               prior(normal(0, 1),class = b,coef = avg_invert),
               prior(normal(0, 1),class = b,coef = group_size),
               prior(normal(15.5, 3),class = Intercept),
               prior(student_t(3, 0, 1), class = sd),
               prior(student_t(3, 0, 1), class = sigma))

m_h2_genomic <- brm(
  BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
    minutes_s + avg_invert + group_size +
    (1 | gr(animal, cov = G1a_sub_final)) + (1 | BirdID) + (1 | BirthYear) +
    (1 | CatchYear) + (1 | Observer),
  data    = df_snp,
  data2   = list(G1a_sub_final = G1a_sub_final),
  family  = gaussian(),
  prior   = my_priors,
  warmup  = 1500,
  iter    = 21500,
  chains  = 4,
  cores   = 4,
  threads = threading(4),
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  seed    = 123,
  file    = "m_brms_h2_genomics_full3_20260506",
  backend = "cmdstanr")