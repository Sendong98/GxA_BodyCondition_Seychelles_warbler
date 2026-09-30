library(Reacnorm)
library(brms)
library(tidybayes)
library(bayesplot)
library(ggplot2)
library(patchwork)
library(nadiv)
library(sjPlot)
library(loo)
library(cmdstanr)
library(rstan)
library(posterior)
library(dplyr)
library(purrr)

# Data extract
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")
PrunedPed_bodymass <- read.csv("./datasets/PrunedPed_bodymass_20260503.csv")
Amat <- as.matrix(nadiv::makeA(PrunedPed_bodymass))

df_clean$minutes_s <- scale(df_clean$minutes)
df_clean$animal <- df_clean$BirdID

df_clean$McAge <- scale(df_clean$age_year-mean(df_clean$age_year))
df_clean$McAge_sq <- (I(df_clean$McAge)^2)
df_clean$McAge_cubic <- (I(df_clean$McAge)^3)

#set up for different gxa orders
gxa_ctrl <- list(adapt_delta = 0.99, max_treedepth = 12)

gxa_priors_reactnorm <- c(
  prior(normal(0, 1),       class = b),
  prior(normal(15.5, 3),    class = Intercept),
  prior(student_t(3, 0, 1), class = sd),
  prior(normal(0, 1),       class = b, dpar = sigma),
  prior(student_t(3, 0, 1), class = sd, coef = Intercept, group = animal),
  prior(student_t(3, 0, 1), class = sd, coef = Intercept, group = BirdID),
  prior(lkj(4),             class = cor, group = animal))

m_reactnorm <- brm(bf(BodyMass ~ McAge + McAge_sq + RightTarsus + SexEstimate +
                        summer + minutes_s + avg_invert + group_size +
                        (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
                        (1 + McAge + McAge_sq | gr(animal, cov = Amat)) +
                        (1 | BirdID),
                      sigma ~ 0 + ageclass4),
                   data    = df_clean,
                   data2   = list(Amat = Amat),
                   family  = gaussian(),
                   prior   = gxa_priors_reactnorm,
                   warmup  = 2000, iter = 32000, chains = 4, cores = 6, seed = 4821,
                   control = gxa_ctrl,
                   file    = "m_Reactnorm_20260911", backend = "cmdstanr")


summary(m_reactnorm)
plot(m_reactnorm)

m_reactnorm_cubic <- brm(bf(BodyMass ~ McAge + McAge_sq + McAge_cubic + RightTarsus + SexEstimate +
                              summer + minutes_s + avg_invert + group_size +
                              (1 | BirthYear) + (1 | CatchYear) + (1 | Observer) +
                              (1 + McAge + McAge_sq + McAge_cubic | gr(animal, cov = Amat)) +
                              (1 | BirdID)),
                         data    = df_clean,
                         data2   = list(Amat = Amat),
                         family  = gaussian(),
                         prior   = gxa_priors_reactnorm,
                         warmup  = 3000, iter = 33000, chains = 4, cores = 6, seed = 4821,
                         control = gxa_ctrl,
                         file    = "m_Reactnorm_scaled_long_cubic_20260922", backend = "cmdstanr")


summary(m_reactnorm_cubic)
plot(m_reactnorm_cubic)

#plot the prediction of the model
tbl_mod <- df_clean %>%
  dplyr::mutate(Predict = predict(m_reactnorm, re_formula = NA) %>%
                  dplyr::as_tibble()) %>%
  tidyr::unpack(Predict) %>%
  dplyr::select(McAge,
                Predict = Estimate,
                Predict_Low = Q2.5,
                Predict_Up = Q97.5) %>%
  dplyr::summarise(dplyr::across(dplyr::starts_with("Predict"), mean),
                   .by = McAge)

p_rn <- ggplot(data = tbl_mod,
               mapping = aes(x = McAge, y = Predict))+geom_ribbon(data = tbl_mod,
                                                                  mapping = aes(x = McAge, ymin = Predict_Low, ymax = Predict_Up),
                                                                  alpha = 0.3) +
  geom_smooth(linewidth = 1)

#plot the prediction of the model
tbl_mod2 <- df_clean %>%
  dplyr::mutate(Predict = predict(m_reactnorm_cubic, re_formula = NA) %>%
                  dplyr::as_tibble()) %>%
  tidyr::unpack(Predict) %>%
  dplyr::select(McAge,
                Predict = Estimate,
                Predict_Low = Q2.5,
                Predict_Up = Q97.5) %>%
  dplyr::summarise(dplyr::across(dplyr::starts_with("Predict"), mean),
                   .by = McAge)

p_rn2 <- ggplot(data = tbl_mod2,
                mapping = aes(x = McAge, y = Predict))+geom_ribbon(data = tbl_mod2,
                                                                   mapping = aes(x = McAge, ymin = Predict_Low, ymax = Predict_Up),
                                                                   alpha = 0.3) +
  geom_smooth(linewidth = 1)

#get the point estimation
seq_age <- seq(min(df_clean$McAge), max(df_clean$McAge), length.out = 200)
seq_X <- cbind(1, seq_age, seq_age^2)

theta_sw <- fixef(m_reactnorm, robust = TRUE)[
  c("Intercept", "McAge", "McAge_sq"), "Estimate"]
names(theta_sw) <- c("a", "b", "c")

S_theta_sw <- vcov(m_reactnorm)[
  c("Intercept", "McAge", "McAge_sq"),
  c("Intercept", "McAge", "McAge_sq")]
rownames(S_theta_sw) <- colnames(S_theta_sw) <- c("a", "b", "c")

G_sw <- VarCorr(m_reactnorm, robust = TRUE)[["animal"]][["cov"]][ , "Estimate", ]
rownames(G_sw) <- colnames(G_sw) <- names(theta_sw)

plas_sw <- rn_phi_decomp(theta = theta_sw,
                         X = seq_X,
                         S = S_theta_sw,
                         wt_env = dnorm(seq_age))
plas_sw

plas_sw_pi <- rn_pi_decomp(theta = theta_sw,
                           V_theta = G_sw,
                           env = seq_age,
                           shape = expression(a + b * x + c * x^2),
                           wt_env = dnorm(seq_age))
plas_sw_pi

gen_sw <- rn_gen_decomp(theta = theta_sw,
                        G_theta = G_sw,
                        X = seq_X,
                        wt_env = dnorm(seq_age))
gen_sw

v_tot_sw <- plas_sw[["V_Plas"]] + gen_sw[["V_Add"]] + mean(vr_sw)
v_tot_sw

var_sw <- c(P2 = plas_sw[["V_Plas"]] / v_tot_sw,
            h2_RN = gen_sw[["V_Add"]] / v_tot_sw,
            h2 = gen_sw[["V_A"]] / v_tot_sw,
            h2_I = gen_sw[["V_AxE"]] / v_tot_sw,
            T2 = (plas_sw[["V_Plas"]] + gen_sw[["V_Add"]]) / v_tot_sw)
var_sw

# Posterior of the quadratic params (subset to a, b, c)
theta_post_sw <- fixef(m_reactnorm, summary = FALSE)[, c("Intercept","McAge","McAge_sq")]
colnames(theta_post_sw) <- c("a","b","c")

# Posterior of the G-matrix
G_post_sw <- VarCorr(m_reactnorm, summary = FALSE)[["animal"]][["cov"]] %>%
  apply(1, \(mat_) mat_, simplify = FALSE) %>%
  map(\(mat_) { rownames(mat_) <- colnames(mat_) <- c("a","b","c"); mat_ })

# Assemble draws_df with theta + G list-columns
post_sw <- as_draws_df(theta_post_sw)

post_sw[["theta"]] <- post_sw %>%
  select(a:c) %>%
  apply(1, \(vec_) { vec_ }, simplify = FALSE)

post_sw[["G"]] <- G_post_sw

post_sw <- thin_draws(post_sw, thin = nrow(theta_post_sw) / 1000)

post_sw_info <- select(post_sw, starts_with("."))
post_plas_sw <- map(post_sw[["theta"]],
                    \(th_) { rn_phi_decomp(theta = th_,
                                           X = seq_X,
                                           S = S_theta_sw,
                                           wt_env = dnorm(seq_age)) },
                    .progress = TRUE) %>%
  bind_rows() %>%
  select(where(\(col_) { abs(mean(col_)) > 10^-5 })) %>%
  # Transform into a "draws" object using posterior package
  cbind(post_sw_info) %>%
  as_draws_df()

summarise_draws(post_plas_sw)

mcmc_trace(post_plas_sw)

mcmc_areas(post_plas_sw,
           regex_pars = "^V",
           prob = 0.95,
           area_method = "scaled height") /
  mcmc_areas(post_plas_sw,
             regex_pars = "^[^V]",
             prob = 0.95,
             area_method = "scaled height") +
  plot_layout(heights = c(1, 2))

post_gen_sw <- map2(post_sw[["theta"]], post_sw[["G"]],
                    \(th_, G_) { rn_gen_decomp(theta = th_,
                                               G_theta = G_,
                                               X = seq_X,
                                               wt_env = dnorm(seq_age)) },
                    .progress = TRUE) %>%
  bind_rows() %>%
  select(where(\(col_) { abs(mean(col_)) > 10^-5 })) %>%
  cbind(post_sw_info) %>%
  as_draws_df()

summarise_draws(post_gen_sw)
mcmc_trace(post_gen_sw)
mcmc_areas(post_gen_sw,
           regex_pars = "^V",
           prob = 0.95,
           area_method = "scaled height")

min(extract_variable(post_gen_sw, "V_AxE"))


post_gen_sw %>%
  mutate(prop_GxE = V_AxE / V_Add) %>%
  posterior::summarise_draws("median", ~quantile(.x, c(0.025, 0.975)))

#reaction norm with cubic genetic random slope
#get the point estimation
seq_age <- seq(min(df_clean$McAge), max(df_clean$McAge), length.out = 200)
seq_X <- cbind(1, seq_age, seq_age^2, seq_age^3)

theta_sw_cubic <- fixef(m_reactnorm_cubic, robust = TRUE)[
  c("Intercept", "McAge", "McAge_sq", "McAge_cubic"), "Estimate"]
names(theta_sw_cubic) <- c("a", "b", "c", "d")

S_theta_sw_cubic <- vcov(m_reactnorm_cubic)[
  c("Intercept", "McAge", "McAge_sq","McAge_cubic"),
  c("Intercept", "McAge", "McAge_sq", "McAge_cubic")]
rownames(S_theta_sw_cubic) <- colnames(S_theta_sw_cubic) <- c("a", "b", "c", "d")

G_sw_cubic <- VarCorr(m_reactnorm_cubic, robust = TRUE)[["animal"]][["cov"]][ , "Estimate", ]
rownames(G_sw_cubic) <- colnames(G_sw_cubic) <- c("a", "b", "c","d")

plas_sw <- rn_phi_decomp(theta = theta_sw_cubic,
                         X = seq_X,
                         S = S_theta_sw_cubic,
                         wt_env = dnorm(seq_age))
plas_sw

plas_sw_pi <- rn_pi_decomp(theta = theta_sw_cubic,
                           V_theta = G_sw_cubic,
                           env = seq_age,
                           shape = expression(a + b * x + c * x^2 + d * x^3),
                           wt_env = dnorm(seq_age))
plas_sw_pi

gen_sw <- rn_gen_decomp(theta = theta_sw_cubic,
                        G_theta = G_sw_cubic,
                        X = seq_X,
                        wt_env = dnorm(seq_age))
gen_sw

# Posterior of the quadratic params (subset to a, b, c, d)
theta_post_sw_cubic <- fixef(m_reactnorm_cubic, summary = FALSE)[, c("Intercept","McAge","McAge_sq","McAge_cubic")]
colnames(theta_post_sw_cubic) <- c("a","b","c","d")

# Posterior of the G-matrix
G_post_sw_cubic <- VarCorr(m_reactnorm_cubic, summary = FALSE)[["animal"]][["cov"]] %>%
  apply(1, \(mat_) mat_, simplify = FALSE) %>%
  map(\(mat_) { rownames(mat_) <- colnames(mat_) <- c("a","b","c","d"); mat_ })

# Assemble draws_df with theta plus G list-columns
post_sw_cubic <- as_draws_df(theta_post_sw_cubic)

post_sw_cubic[["theta"]] <- post_sw_cubic %>%
  select(a:d) %>%
  apply(1, \(vec_) { vec_ }, simplify = FALSE)

post_sw_cubic[["G"]] <- G_post_sw_cubic

post_sw_cubic <- thin_draws(post_sw_cubic, thin = nrow(theta_post_sw_cubic) / 1000)

post_sw_cubic_info <- select(post_sw_cubic, starts_with("."))
post_plas_sw_cubic <- map(post_sw_cubic[["theta"]],
                          \(th_) { rn_phi_decomp(theta = th_,
                                                 X = seq_X,
                                                 S = S_theta_sw_cubic,
                                                 wt_env = dnorm(seq_age)) },
                          .progress = TRUE) %>%
  bind_rows() %>%
  select(where(\(col_) { abs(mean(col_)) > 10^-5 })) %>%
  # Transform into a "draws" object using posterior package
  cbind(post_sw_cubic_info) %>%
  as_draws_df()

summarise_draws(post_plas_sw_cubic)

mcmc_trace(post_plas_sw_cubic)

mcmc_areas(post_plas_sw_cubic,
           regex_pars = "^V",
           prob = 0.95,
           area_method = "scaled height") /
  mcmc_areas(post_plas_sw_cubic,
             regex_pars = "^[^V]",
             prob = 0.95,
             area_method = "scaled height") +
  plot_layout(heights = c(1, 2))

post_gen_sw_cubic <- map2(post_sw_cubic[["theta"]], post_sw_cubic[["G"]],
                          \(th_, G_) { rn_gen_decomp(theta = th_,
                                                     G_theta = G_,
                                                     X = seq_X,
                                                     wt_env = dnorm(seq_age)) },
                          .progress = TRUE) %>%
  bind_rows() %>%
  select(where(\(col_) { abs(mean(col_)) > 10^-5 })) %>%
  cbind(post_sw_cubic_info) %>%
  as_draws_df()

summarise_draws(post_gen_sw_cubic)
mcmc_trace(post_gen_sw_cubic)

post_renamed <- post_gen_sw_cubic %>%
  rename_variables(
    V_Add = V_Add,
    V_G   = V_A,
    V_GxA = V_AxE
  )

mcmc_areas(post_renamed,
           pars = c("V_Add", "V_G", "V_GxA"),
           prob = 0.95,
           area_method = "scaled height")

post_gen_sw_cubic %>%
  mutate(prop_GxE = V_AxE / V_Add) %>%
  posterior::summarise_draws("median", ~quantile(.x, c(0.025, 0.975)))

