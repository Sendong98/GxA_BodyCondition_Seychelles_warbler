library(tidyverse)
library(lme4)
library(brms)
library(tidybayes)
library(bayesplot)
library(nadiv)
library(rstan)

# Data
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")
PrunedPed_bodymass <- read.csv("./datasets/PrunedPed_bodymass_20260503.csv")

df_clean$minutes_s <- scale(df_clean$minutes)
df_clean$animal <- df_clean$BirdID

dataset_mass_Ped <- df_clean %>%
  filter(BirdID %in% PrunedPed_bodymass$BirdID)

#Additive-relationship matrix (A matrix)

Ped <- PrunedPed_bodymass %>%
  dplyr::select(id = BirdID, dam = Dam, sire = Sire)

Amat <- as.matrix(nadiv::makeA(Ped))


#Remove fixed effects (remove_fixed)

remove_fixed <- lmer(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate +
                       summer + minutes_s + avg_invert + group_size +
                       (1 | BirthYear) + (1 | CatchYear) + (1 | Observer),
                     data = dataset_mass_Ped)
summary(remove_fixed)


# Prepare character-state dataset

Mass_Long <- dataset_mass_Ped %>%
  mutate(Mass_resid = resid(remove_fixed),
         AgeClass = case_when(age_year < 1 ~ "Juv",
                              age_year >= 1 & age_year <= 5 ~ "Young",
                              age_year > 5 ~ "Old")) %>%
  group_by(BirdID, AgeClass) %>%
  summarise(MeanResidMass = mean(Mass_resid), .groups = "drop")

Mass_Wide <- Mass_Long %>%
  pivot_wider(names_from = AgeClass, values_from = MeanResidMass) %>%
  relocate(BirdID, Juv, Young, Old)

Mass_Ped_CS <- dataset_mass_Ped %>%
  dplyr::select(BirdID, animal) %>%
  distinct(.keep_all = TRUE) %>%
  left_join(Mass_Wide, by = "BirdID")

# Missingness check
Mass_Ped_CS %>%
  summarise(across(c(Juv, Young, Old),
                   list(n = ~sum(!is.na(.)),
                        percentage_missing = ~round(mean(is.na(.)) * 100, 1)))) %>% print()
write.csv(Mass_Ped_CS, file = "./datasets/Mass_Ped_CS_20260505.csv", row.names = FALSE)

# Priors
priors_CS <- c(prior(student_t(3, 0, 1),   class = sd,    resp = "Juv"),
               prior(student_t(3, 0, 1),   class = sd,    resp = "Young"),
               prior(student_t(3, 0, 1),   class = sd,    resp = "Old"),
               prior(student_t(3, 0, 0.5), class = sigma, resp = "Juv"),
               prior(student_t(3, 0, 0.5), class = sigma, resp = "Young"),
               prior(student_t(3, 0, 0.5), class = sigma, resp = "Old"),
               prior(lkj(2), class = cor,   group = "animal"),  # genetic correlations between age classes
               prior(lkj(2), class = rescor)                    # residual correlations between age classes
)

# Character-state model
f_juv   <- bf(Juv   | mi() ~ 1 + (1 | p | gr(animal, cov = Amat)))
f_young <- bf(Young | mi() ~ 1 + (1 | p | gr(animal, cov = Amat)))
f_old   <- bf(Old   | mi() ~ 1 + (1 | p | gr(animal, cov = Amat)))

# character-state model
model_CS <- brm(
  mvbf(f_juv, f_young, f_old, rescor = TRUE),
  data    = Mass_Ped_CS,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = priors_CS,
  warmup  = 2000,
  iter    = 32000,
  chains  = 4,
  cores   = 4,
  threads = threading(16),
  seed    = 4817,
  normalize = FALSE,
  control = list(adapt_delta=0.95,stepsize=0.001, max_treedepth=15),
  file    = "model_CS_longer_20260506",
  backend = "cmdstanr")

# Convergence diagnostics
summary(model_CS)

#mcmc_plot(model_CS, type = "rhat_hist")
#mcmc_plot(model_CS, type = "neff_hist")

#mcmc_plot(model_CS, type = "trace", variable = "^sd_animal", regex = TRUE)
#mcmc_plot(model_CS, type = "trace", variable = "^cor_animal", regex = TRUE)

#mcmc_plot(model_CS, type = "areas", variable = "^sd_", regex = TRUE) +
#  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
#  labs(title = "Random effect SDs", x = "Posterior estimate")

# Additive genetic variances (VA) and h² per age class

age_class <- c("Juv", "Young", "Old")

VA_draws <- model_CS %>%
  spread_draws(
    sd_animal__Juv_Intercept,
    sd_animal__Young_Intercept,
    sd_animal__Old_Intercept,
    sigma_Juv, sigma_Young, sigma_Old
  ) %>%
  mutate(
    VA_Juv   = sd_animal__Juv_Intercept^2,
    VA_Young = sd_animal__Young_Intercept^2,
    VA_Old   = sd_animal__Old_Intercept^2,
    Ve_Juv   = sigma_Juv^2,
    Ve_Young = sigma_Young^2,
    Ve_Old   = sigma_Old^2,
    Vp_Juv   = VA_Juv + Ve_Juv,
    Vp_Young = VA_Young + Ve_Young,
    Vp_Old   = VA_Old + Ve_Old,
    h2_Juv   = VA_Juv   / Vp_Juv,
    h2_Young = VA_Young / Vp_Young,
    h2_Old   = VA_Old   / Vp_Old
  )

VA_summary <- VA_draws %>%
  select(.draw, starts_with("VA_"), starts_with("h2_")) %>%
  pivot_longer(-.draw, names_to = "param", values_to = "value") %>%
  separate(param, into = c("metric", "AgeClass"), sep = "_(?=[^_]+$)") %>%
  group_by(metric, AgeClass) %>%
  mean_hdi(value) %>%
  mutate(AgeClass = factor(AgeClass, levels = age_class))
print(VA_summary)

p_va_h2 <- ggplot(VA_summary, aes(x = AgeClass, y = value)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.15) +
  facet_wrap(~ metric, scales = "free_y") +
  labs(
    x = "Age class",
    y = "Estimate",
    title = "Additive genetic variance (VA) and heritability (h²) across age classes"
  ) +
  theme_classic(base_size = 12)

ggsave("./p_va_h2_cs_20260502.pdf", plot = p_va_h2, width = 300, height = 150, units = "mm", dpi = 300)

# Genetic correlations between age classes

gen_cor_draws <- model_CS %>%
  spread_draws(
    cor_animal__Juv_Intercept__Young_Intercept,
    cor_animal__Juv_Intercept__Old_Intercept,
    cor_animal__Young_Intercept__Old_Intercept
  )

gen_cor_summary <- gen_cor_draws %>%
  select(.draw, starts_with("cor_animal")) %>%
  pivot_longer(-.draw, names_to = "parameter", values_to = "r") %>%
  group_by(parameter) %>%
  mean_hdi(r) %>%
  mutate(parameter = gsub("cor_animal__|_Intercept", "", parameter))
print(gen_cor_summary)


# Genetic correlation heatmap

age_pairs   <- list(c("Juv", "Young"), c("Juv", "Old"), c("Young", "Old"))

cor_var_names <- c(
  "cor_animal__Juv_Intercept__Young_Intercept",
  "cor_animal__Juv_Intercept__Old_Intercept",
  "cor_animal__Young_Intercept__Old_Intercept"
)

r_mat <- matrix(1, nrow = 3, ncol = 3, dimnames = list(age_class,age_class))

for (i in seq_along(age_pairs)) {
  r_med <- median(gen_cor_draws[[cor_var_names[i]]])
  r_mat[age_pairs[[i]][1], age_pairs[[i]][2]] <- r_med
  r_mat[age_pairs[[i]][2], age_pairs[[i]][1]] <- r_med
}

p_heatmap <- as.data.frame(r_mat) %>%
  rownames_to_column("Age1") %>%
  pivot_longer(-Age1, names_to = "Age2", values_to = "r") %>%
  mutate(
    Age1 = factor(Age1, levels = age_class),
    Age2 = factor(Age2, levels = age_class)
  ) %>%
  ggplot(aes(x = Age2, y = Age1, fill = r)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(r, 2)), size = 3.5) +
  scale_fill_gradient2(low = "#003f5c", mid = "white", high = "#ffa600",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Genetic\ncorrelation") +
  labs(x = NULL, y = NULL,
       title = "Genetic correlations across age classes (model_CS)") +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(size = 11))

ggsave("./p_heatmap_cs_20260502.pdf", plot = p_heatmap, width = 150, height = 150, units = "mm", dpi = 300)

