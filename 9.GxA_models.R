library(brms)
library(tidybayes)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(nadiv)
library(sjPlot)
library(loo)
library(cmdstanr)

# Data extract
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")
PrunedPed_bodymass <- read.csv("./datasets/PrunedPed_bodymass_20260503.csv")
Amat <- as.matrix(nadiv::makeA(PrunedPed_bodymass))

df_clean$minutes_s <- scale(df_clean$minutes)
df_clean$animal <- df_clean$BirdID

df_clean$RsAge <- -1 + (1 - (-1)) / (max(df_clean$age_year) - min(df_clean$age_year)) *
  (df_clean$age_year - min(df_clean$age_year))

nblocks4 <- 4
df_clean$ageclass4 <- as.numeric(
  arules::discretize(df_clean$age_year, breaks = nblocks4, method = "frequency"))
df_clean$ageclass4 <- as.factor(df_clean$ageclass4)

#set up for different gxa orders
gxa_ctrl <- list(adapt_delta = 0.99, max_treedepth = 12)

#Base priors (fixed effects, overall sd, sigma)
gxa_priors_base <- c(
  prior(normal(0, 1),       class = b),
  prior(normal(15.5, 3),    class = Intercept),
  prior(student_t(3, 0, 1), class = sd),
  prior(normal(0, 1),       class = b, dpar = sigma)
)


# Baseline models — intercept-only for one component

#animal intercept only, PE linear
m_GxA_A0_PE1 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 | gr(animal, cov = Amat)) +
       (1 + RsAge | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept, group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,     group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A0_PE1_20260503", backend = "cmdstanr"
)

#animal linear, PE intercept only
m_GxA_A1_PE0 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge | gr(animal, cov = Amat)) +
       (1 | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept, group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,     group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A1_PE0_20260503", backend = "cmdstanr"
)

# A1 models — animal: (1 + RsAge | gr(animal, cov = Amat))

#A1&PE1
m_GxA_A1_PE1 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge | gr(animal, cov = Amat)) +
       (1 + RsAge | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept, group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,     group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept, group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,     group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A1_PE1_20260503", backend = "cmdstanr")

#A1&PE2
m_GxA_A1_PE2 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge | gr(animal, cov = Amat)) +
       (1 + RsAge + I(RsAge^2) | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,      group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,      group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A1_PE2_20260503", backend = "cmdstanr"
)

#  A1  PE3 
m_GxA_A1_PE3 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge | gr(animal, cov = Amat)) +
       (1 + RsAge + I(RsAge^2) + I(RsAge^3) | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,      group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,      group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE3,   group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A1_PE3_20260503", backend = "cmdstanr"
)

# A2 models — animal: (1 + RsAge + I(RsAge^2) | gr(animal, cov = Amat))
#  A2  PE1 
m_GxA_A2_PE1 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge + I(RsAge^2) | gr(animal, cov = Amat)) +
       (1 + RsAge | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,      group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,      group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A2_PE1_20260503", backend = "cmdstanr"
)

#  A2  PE2 
m_GxA_A2_PE2 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge + I(RsAge^2) | gr(animal, cov = Amat)) +
       (1 + RsAge + I(RsAge^2) | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,      group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,      group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A2_PE2_20260503", backend = "cmdstanr"
)

#  A2  PE3 
m_GxA_A2_PE3 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge + I(RsAge^2) | gr(animal, cov = Amat)) +
       (1 + RsAge + I(RsAge^2) + I(RsAge^3) | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,      group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,      group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE3,   group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A2_PE3_20260503", backend = "cmdstanr"
)


# A3 models — animal: (1 + RsAge + I(RsAge^2) + I(RsAge^3) | gr(animal, cov = Amat))
#  A3  PE1 
m_GxA_A3_PE1 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge + I(RsAge^2) + I(RsAge^3) | gr(animal, cov = Amat)) +
       (1 + RsAge | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,      group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE3,   group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,      group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A3_PE1_20260503", backend = "cmdstanr"
)
#A3  PE2 
m_GxA_A3_PE2 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge + I(RsAge^2) + I(RsAge^3) | gr(animal, cov = Amat)) +
       (1 + RsAge + I(RsAge^2) | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,      group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE3,   group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,      group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A3_PE2_20260503", backend = "cmdstanr"
)

#A3  PE3
m_GxA_A3_PE3 <- brm(
  bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
       summer + minutes_s + avg_invert + group_size +
       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
       (1 + RsAge + I(RsAge^2) + I(RsAge^3) | gr(animal, cov = Amat)) +
       (1 + RsAge + I(RsAge^2) + I(RsAge^3) | BirdID),
     sigma ~ 0 + ageclass4),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = c(gxa_priors_base,
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = animal),
              prior(exponential(2),     class = sd, coef = RsAge,      group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = animal),
              prior(exponential(2),     class = sd, coef = IRsAgeE3,   group = animal),
              prior(lkj(4),             class = cor, group = animal),
              prior(student_t(3, 0, 1), class = sd, coef = Intercept,  group = BirdID),
              prior(exponential(2),     class = sd, coef = RsAge,      group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE2,   group = BirdID),
              prior(exponential(2),     class = sd, coef = IRsAgeE3,   group = BirdID),
              prior(lkj(4),             class = cor, group = BirdID)),
  warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
  control = gxa_ctrl,
  file    = "m_GxA_A3_PE3_20260503", backend = "cmdstanr"
)

# LOO-CV model comparison

m_GxA_A0_PE1 <- add_criterion(m_GxA_A0_PE1, "loo")
m_GxA_A1_PE0 <- add_criterion(m_GxA_A1_PE0, "loo")
m_GxA_A1_PE1 <- add_criterion(m_GxA_A1_PE1, "loo")
m_GxA_A1_PE2 <- add_criterion(m_GxA_A1_PE2, "loo")
m_GxA_A1_PE3 <- add_criterion(m_GxA_A1_PE3, "loo")
m_GxA_A2_PE1 <- add_criterion(m_GxA_A2_PE1, "loo")
m_GxA_A2_PE2 <- add_criterion(m_GxA_A2_PE2, "loo")
m_GxA_A2_PE3 <- add_criterion(m_GxA_A2_PE3, "loo")
m_GxA_A3_PE1 <- add_criterion(m_GxA_A3_PE1, "loo")
m_GxA_A3_PE2 <- add_criterion(m_GxA_A3_PE2, "loo")
m_GxA_A3_PE3 <- add_criterion(m_GxA_A3_PE3, "loo")

loo_compare(m_GxA_A0_PE1, m_GxA_A1_PE0,
  m_GxA_A1_PE1, m_GxA_A1_PE2, m_GxA_A1_PE3,
  m_GxA_A2_PE1, m_GxA_A2_PE2, m_GxA_A2_PE3,
  m_GxA_A3_PE1, m_GxA_A3_PE2, m_GxA_A3_PE3)
