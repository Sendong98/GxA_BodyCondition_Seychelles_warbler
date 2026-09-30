#load packages
library(brms)
library(loo)
library(tidybayes)
library(tidyverse)
library(patchwork)
options(brms.backend = "cmdstanr")

#IxA
m_ixa0 <- read_rds("m_IxA_0_20260506.rds")
m_ixa1 <- read_rds("m_IxA_1_20260506.rds")
m_ixa2 <- read_rds("m_IxA_2_20260503.rds")

check_hmc_diagnostics(m_ixa0)
summary(m_ixa0)
check_hmc_diagnostics(m_ixa1)
summary(m_ixa1)
check_hmc_diagnostics(m_ixa2)
summary(m_ixa2)

m_ixa0 <- add_criterion(m_ixa0, "loo")
m_ixa1 <- add_criterion(m_ixa1, "loo")
m_ixa2 <- add_criterion(m_ixa2, "loo")
loo_compare(m_ixa0, m_ixa1, m_ixa2)

m_ixa0 <- add_criterion(m_ixa0, "waic")
m_ixa1 <- add_criterion(m_ixa1, "waic")
m_ixa2 <- add_criterion(m_ixa2, "waic")
loo_compare(m_ixa0, m_ixa1, m_ixa2, criterion = "waic")

#plot ixa2
# Sample draws and birds
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")
df_clean$RsAge <- -1+(1-(-1))/(max(df_clean$age_year) - min(df_clean$age_year)) * (df_clean$age_year - min(df_clean$age_year))

# Unique age grid from df_clean
age_grid <- df_clean %>%
  distinct(RsAge, age_year) %>%
  arrange(age_year)

local_Phi <- cbind(
  1/sqrt(2),
  sqrt(3/2) * age_grid$RsAge,
  sqrt(5/2) * 0.5 * (3 * age_grid$RsAge^2 - 1)
)
re_draws <- spread_draws(m_ixa2, r_animal[BirdID, term])
unique(re_draws$term)
n_draws_sample <- 200
sampled_draws  <- sample(unique(re_draws$.draw), n_draws_sample)
sampled_birds  <- sample(unique(re_draws$BirdID), 300)

# Filter and pivot once
re_wide <- re_draws %>%
  filter(.draw %in% sampled_draws, BirdID %in% sampled_birds) %>%
  tidyr::pivot_wider(names_from = term, values_from = r_animal)

# Single matrix multiply
coef_mat <- as.matrix(re_wide[, c("Intercept", "RsAge", "IRsAgeE2")])
BV_mat   <- coef_mat %*% t(local_Phi)

#Flatten to long
bv_traj <- re_wide %>%
  dplyr::select(BirdID, .draw) %>%
  bind_cols(as.data.frame(BV_mat)) %>%
  tidyr::pivot_longer(
    cols      = -c(BirdID, .draw),
    names_to  = "age_idx",
    values_to = "BV"
  ) %>%
  mutate(Age = age_grid$age_year[as.integer(sub("V", "", age_idx))]) %>%
  dplyr::select(BirdID, .draw, Age, BV)

# Summarise across draws per bird×age
bv_summary <- bv_traj %>%
  group_by(BirdID, Age) %>%
  summarise(
    BV_mean  = mean(BV),
    BV_lower = quantile(BV, 0.025),
    BV_upper = quantile(BV, 0.975),
    .groups  = "drop"
  )
#Extract posterior means of fixed effects for age
fix_draws <- spread_draws(m_ixa2, b_age_year, b_Iage_yearE2)

fix_sampled <- fix_draws %>%
  dplyr::select(.draw, b_age_year, b_Iage_yearE2)

# Merge fixed + random effects
bv_traj_full <- bv_traj %>%
  left_join(fix_sampled, by = ".draw") %>%
  left_join(
    df_clean %>% distinct(age_year, RsAge),
    by = c("Age" = "age_year")
  ) %>%
  mutate(
    # Full trajectory = population mean age curve + individual deviation
    Traj = BV + b_age_year * Age + b_Iage_yearE2 * Age^2
  )

# Summarise
traj_summary <- bv_traj_full %>%
  group_by(BirdID, Age) %>%
  summarise(
    Traj_mean  = mean(Traj),
    Traj_lower = quantile(Traj, 0.025),
    Traj_upper = quantile(Traj, 0.975),
    .groups    = "drop"
  )

# Population mean trajectory across all sampled birds
pop_mean <- traj_summary %>%
  group_by(Age) %>%
  summarise(Pop_mean = mean(Traj_mean))

# Plot IxA trajectories (reaction norm)
p_ixa <- ggplot(traj_summary,
                aes(x = Age, y = Traj_mean, group = BirdID)) +
  geom_line(alpha = 0.2, linewidth = 0.3, colour = "#808000") +
  # Population mean
  geom_line(data = pop_mean,
            aes(x = Age, y = Pop_mean, group = 1),
            colour = "black", linewidth = 1.2, inherit.aes = FALSE) +
  labs(x = "Age (years)", y = "Body condition trajectory") +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.line         = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_line(linewidth = 0.4, colour = "black"),
    axis.text.x       = element_text(size = 10, colour = "black"),
    axis.text.y       = element_text(size = 14, colour = "black"),
    axis.title        = element_text(size = 16, colour = "black"),
    legend.position   = "none",
    panel.border      = element_blank())

p_ixa
ggsave("p_ixa20260509.pdf",plot = p_ixa, width = 200, height = 150, units = "mm", dpi = 300)
ggsave("p_ixa20260509.jpg",plot = p_ixa, width = 200, height = 150, units = "mm", dpi = 300)

#GxA
model_files <- list(
  m_gxa00 = "m_brms_h2_1_20260503.rds",
  m_gxa01 = "m_IxA_1_20260506.rds",
  m_gxa10 = "m_GxA_A1_PE0_20260506_99_longer.rds",
  m_gxa11 = "m_GxA_A1_PE1_20260506_99_longer.rds",
  m_gxa20 = "m_GxA_A2_PE0_20260506_99_longer.rds",
  m_gxa21 = "m_GxA_A2_PE1_20260506_99_longer.rds",
  m_gxa22 = "m_GxA_A2_PE2_20260506_99_longer.rds",
  m_gxa30 = "m_GxA_A3_PE0_20260504_longer.rds",
  m_gxa31 = "m_GxA_A3_PE1_20260506_99_longer.rds",
  m_gxa32 = "m_GxA_A3_PE2_20260506_99_longer.rds",
  m_gxa33 = "m_GxA_A3_PE3_20260504_longer.rds"
)

# Compute and save LOO/WAIC
for (name in names(model_files)) {
  cat("Processing", name, "\n")
  
  m <- read_rds(model_files[[name]])
  m <- add_criterion(m, c("loo", "waic"))
  
  saveRDS(m$criteria, file = paste0("criteria_", name, ".rds"))
  
  rm(m)
  gc()
}

# Load criteria back
criteria_list <- lapply(names(model_files), function(name) {
  readRDS(paste0("criteria_", name, ".rds"))
})
names(criteria_list) <- names(model_files)

# Extract LOO objects
loo_list <- lapply(criteria_list, `[[`, "loo")

# Compare
loo_compare(loo_list)

# WAIC - not accurate
waic_list <- lapply(criteria_list, `[[`, "waic")
loo_compare(waic_list)

#plot gxa30
m_gxa30 <- read_rds("m_GxA_A3_PE0_20260504_longer.rds")
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")
df_clean$RsAge <- -1+(1-(-1))/(max(df_clean$age_year) - min(df_clean$age_year)) * (df_clean$age_year - min(df_clean$age_year))
re_draws <- spread_draws(m_gxa30, r_animal[animal, term])
unique(re_draws$term)

# Age grid
age_grid <- df_clean %>%
  distinct(RsAge, age_year) %>%
  arrange(age_year)

# Phi matrix
local_Phi <- cbind(
  1/sqrt(2),
  sqrt(3/2) * age_grid$RsAge,
  sqrt(5/2) * 0.5 * (3 * age_grid$RsAge^2 - 1),
  sqrt(7/2) * 0.5 * (5 * age_grid$RsAge^3 - 3 * age_grid$RsAge)
)

# Extract sd and cor draws
vcv_draws <- as_draws_df(m_gxa30) %>%
  dplyr::select(
    .draw,
    sd_animal__Intercept, sd_animal__RsAge, sd_animal__IRsAgeE2, sd_animal__IRsAgeE3,
    cor_animal__Intercept__RsAge, cor_animal__Intercept__IRsAgeE2,
    cor_animal__Intercept__IRsAgeE3, cor_animal__RsAge__IRsAgeE2,
    cor_animal__RsAge__IRsAgeE3, cor_animal__IRsAgeE2__IRsAgeE3
  )

# Sample draws for speed
set.seed(42)
sampled_idx <- sample(1:nrow(vcv_draws), 500)
vcv_sample  <- vcv_draws[sampled_idx, ]

# Function to build G from sd + cor draws
build_G <- function(row) {
  sd  <- c(row$sd_animal__Intercept, row$sd_animal__RsAge,
           row$sd_animal__IRsAgeE2,  row$sd_animal__IRsAgeE3)
  
  # Correlation matrix
  R <- diag(4)
  R[1,2] <- R[2,1] <- row$cor_animal__Intercept__RsAge
  R[1,3] <- R[3,1] <- row$cor_animal__Intercept__IRsAgeE2
  R[1,4] <- R[4,1] <- row$cor_animal__Intercept__IRsAgeE3
  R[2,3] <- R[3,2] <- row$cor_animal__RsAge__IRsAgeE2
  R[2,4] <- R[4,2] <- row$cor_animal__RsAge__IRsAgeE3
  R[3,4] <- R[4,3] <- row$cor_animal__IRsAgeE2__IRsAgeE3
  
  diag(sd) %*% R %*% diag(sd)
}

# Compute Vg per draw per age
Vg_draws <- lapply(1:nrow(vcv_sample), function(d) {
  G_d <- build_G(vcv_sample[d, ])
  vg  <- sapply(1:nrow(age_grid), function(i) {
    phi <- local_Phi[i, ]
    as.numeric(phi %*% G_d %*% phi)
  })
  data.frame(draw = d, Age = age_grid$age_year, Vg = vg)
})

Vg_all <- bind_rows(Vg_draws)

# Summarise
Vg_df <- Vg_all %>%
  group_by(Age) %>%
  summarise(
    GeneticVariance = mean(Vg),
    Vg_lower        = quantile(Vg, 0.025),
    Vg_upper        = quantile(Vg, 0.975),
    .groups         = "drop"
  )

# Plot pA
pA <- ggplot(Vg_df, aes(x = Age, y = GeneticVariance)) +
  geom_ribbon(aes(ymin = Vg_lower, ymax = Vg_upper),
              fill = "#808000", alpha = 0.2) +
  geom_line(linewidth = 1.2, colour = "#808000") +
  labs(x = "Age (yr)", y = "Additive genetic variance") +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.line         = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_line(linewidth = 0.4, colour = "black"),
    axis.text.x       = element_text(size = 10, colour = "black"),
    axis.text.y       = element_text(size = 14, colour = "black"),
    axis.title        = element_text(size = 16, colour = "black"),
    legend.position   = "none",
    panel.border      = element_blank())
pA

# Sample draws and birds
set.seed(42)
n_draws_sample <- 200
sampled_draws  <- sample(unique(re_draws$.draw), n_draws_sample)
sampled_birds  <- sample(unique(re_draws$animal), 300)

# Filter and pivot
re_wide <- re_draws %>%
  filter(.draw %in% sampled_draws, animal %in% sampled_birds) %>%
  tidyr::pivot_wider(names_from = term, values_from = r_animal)

# Matrix multiply
coef_mat <- as.matrix(re_wide[, c("Intercept", "RsAge", "IRsAgeE2", "IRsAgeE3")])
BV_mat   <- coef_mat %*% t(local_Phi)

# Flatten to long
bv_traj <- re_wide %>%
  dplyr::select(animal, .draw) %>%
  bind_cols(as.data.frame(BV_mat)) %>%
  tidyr::pivot_longer(
    cols      = -c(animal, .draw),
    names_to  = "age_idx",
    values_to = "BV"
  ) %>%
  mutate(Age = age_grid$age_year[as.integer(sub("V", "", age_idx))]) %>%
  dplyr::select(animal, .draw, Age, BV)

# Summarise per animal x age across draws
bv_summary <- bv_traj %>%
  group_by(animal, Age) %>%
  summarise(
    BV_mean  = mean(BV),
    BV_lower = quantile(BV, 0.025),
    BV_upper = quantile(BV, 0.975),
    .groups  = "drop"
  )

# Plot
pB <- ggplot(bv_summary,
             aes(x = Age, y = BV_mean, group = animal)) +
  geom_line(alpha = 0.15, linewidth = 0.3, colour = "#808000") +
  labs(x = "Age (yr)", y = "Breeding value") +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.line         = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.x      = element_blank(),
    axis.ticks.y      = element_line(linewidth = 0.4, colour = "black"),
    axis.text.x       = element_text(size = 10, colour = "black"),
    axis.text.y       = element_text(size = 14, colour = "black"),
    axis.title        = element_text(size = 16, colour = "black"),
    legend.position   = "none",
    panel.border      = element_blank())

# Combined figure
GxA_figure <- pA / pB +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 14, face = "bold")))

GxA_figure
ggsave("GxA_figure.pdf",plot = GxA_figure, width = 150, height = 200, units = "mm", dpi = 300)
ggsave("GxA_figure.jpg",plot = GxA_figure, width = 150, height = 200, units = "mm", dpi = 300)


