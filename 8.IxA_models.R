library(brms)
library(tidybayes)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(nadiv)
library(sjPlot)
library(loo)
library(cmdstanr)

#Data extract
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")
PrunedPed_bodymass <- read.csv("./datasets/PrunedPed_bodymass_20260503.csv")
Amat <- as.matrix(nadiv::makeA(PrunedPed_bodymass))

df_clean$minutes_s <- scale(df_clean$minutes)
df_clean$animal <- df_clean$BirdID

#Re-scaled age, IxA with hetero residuals
df_clean$RsAge <- -1+(1-(-1))/(max(df_clean$age_year) - min(df_clean$age_year)) * (df_clean$age_year - min(df_clean$age_year))
nblocks4 <-4	# number of 'residual blocks'
df_clean$ageclass4 <- as.numeric(arules::discretize(df_clean$age_year,breaks= nblocks4, method='frequency'))
df_clean$ageclass4 <- as.factor(df_clean$ageclass4)

# IxA polynomial order comparison

ixa_priors_base <- c(prior(normal(0, 1),        class = b),
                     prior(normal(15.5, 3),     class = Intercept),
                     prior(student_t(3, 0, 1),  class = sd),
                     prior(normal(0, 1),        class = b, dpar = sigma))

ixa_ctrl <- list(adapt_delta = 0.99, max_treedepth = 12)

# 0th order: (1 | BirdID), intercept only
m_IxA_0 <- brm(bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
                    summer + minutes_s + avg_invert + group_size +
                    (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) + (1 | BirdID),
                  sigma ~ 0 + ageclass4),
               data    = df_clean,
               family  = gaussian(),
               prior   = ixa_priors_base,
               warmup  = 2000, iter = 32000, chains = 4, cores = 4,threads = threading(4), seed = 123,
               control = ixa_ctrl,
               file    = "m_IxA_0_20260506", backend = "cmdstanr")
#mcmc_plot(m_IxA_0, type = "trace", variable = variables(m_IxA_0)[1:18])
summary(m_IxA_0)

# 1st order: (1 + RsAge | BirdID), linear slope
m_IxA_1 <- brm(bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
                    summer + minutes_s + avg_invert + group_size +
                    (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) + (1 + RsAge | BirdID),
                  sigma ~ 0 + ageclass4),
               data    = df_clean,
               family  = gaussian(),
               prior   = c(ixa_priors_base,
                           prior(exponential(2),       class = sd, coef = RsAge,     group = BirdID),
                           prior(student_t(3, 0, 1),   class = sd, coef = Intercept, group = BirdID),
                           prior(lkj(4),               class = cor)),
               warmup  = 2000, iter = 32000, chains = 4, cores = 4,threads = threading(4), seed = 123,
               control = ixa_ctrl,
               file    = "m_IxA_1_20260506", backend = "cmdstanr")

# 2nd order: (1 + RsAge + I(RsAge^2) | BirdID), quadratic slope
m_IxA_2 <- brm(bf(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
                    summer + minutes_s + avg_invert + group_size +
                    (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) + (1 + RsAge + I(RsAge^2) | BirdID),
                  sigma ~ 0 + ageclass4),
               data    = df_clean,
               family  = gaussian(),
               prior   = c(ixa_priors_base,
                           prior(exponential(2),       class = sd, coef = RsAge,        group = BirdID),
                           prior(exponential(2),       class = sd, coef = IRsAgeE2,   group = BirdID),
                           prior(student_t(3, 0, 1),   class = sd, coef = Intercept,    group = BirdID),
                           prior(lkj(4),               class = cor)),
               warmup  = 2000, iter = 32000, chains = 4, cores = 4,threads = threading(4), seed = 123,
               control = ixa_ctrl,
               file    = "m_IxA_2_20260503", backend = "cmdstanr")

# Compare models via LOO
m_IxA_0 <- add_criterion(m_IxA_0, "loo")
m_IxA_1 <- add_criterion(m_IxA_1, "loo")
m_IxA_2 <- add_criterion(m_IxA_2, "loo")
loo_compare(m_IxA_0, m_IxA_1, m_IxA_2)

