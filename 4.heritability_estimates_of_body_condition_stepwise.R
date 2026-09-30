#heritability estimate - stepwise and using pedigree data
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

#set up the cmdstan path on my laptop
set_cmdstan_path("/Users/sendong/Downloads/cmdstan-2.38.0")
Sys.setenv(PATH = "/usr/bin:/bin:/usr/sbin:/sbin")

# load data
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")
PrunedPed_bodymass <- read.csv("./datasets/PrunedPed_bodymass_20260503.csv")
Amat <- as.matrix(nadiv::makeA(PrunedPed_bodymass))

df_clean$minutes_s <- scale(df_clean$minutes)
df_clean$animal <- df_clean$BirdID

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

# S1 additive genetic only
m_h2_s1 <- brm(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
                 minutes_s + avg_invert + group_size +
                 (1 | gr(animal, cov = Amat)),
               data    = df_clean,
               data2   = list(Amat = Amat),
               family  = gaussian(),
               prior   = my_priors,
               warmup  = 1500,
               iter    = 11500,
               chains  = 4,
               cores   = 16,
               seed    = 123,
               file    = "m_h2_s1_20260505",
               backend = "cmdstanr")

# S2 plus permanent environment (BirdID)
m_h2_s2 <- brm(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
                 minutes_s + avg_invert + group_size +
                 (1 | gr(animal, cov = Amat)) + (1 | BirdID),
               data    = df_clean,
               data2   = list(Amat = Amat),
               family  = gaussian(),
               prior   = my_priors,
               warmup  = 1500,
               iter    = 11500,
               chains  = 4,
               cores   = 6,
               seed    = 123,
               file    = "m_h2_s2_20260505",
               backend = "cmdstanr")

# S3 plus birth year
m_h2_s3 <- brm(BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
                 minutes_s + avg_invert + group_size +
                 (1 | gr(animal, cov = Amat)) + (1 | BirdID) + (1 | BirthYear),
               data    = df_clean,
               data2   = list(Amat = Amat),
               family  = gaussian(),
               prior   = my_priors,
               warmup  = 1500,
               iter    = 11500,
               chains  = 4,
               cores   = 6,
               seed    = 123,
               file    = "m_h2_s3_20260505",
               backend = "cmdstanr")

# S4 plus catch year
m_h2_s4 <- brm(
  BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
    minutes_s + avg_invert + group_size +
    (1 | gr(animal, cov = Amat)) + (1 | BirdID) + (1 | BirthYear) +
    (1 | CatchYear),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = my_priors,
  warmup  = 1500,
  iter    = 11500,
  chains  = 4,
  cores   = 4,
  seed    = 123,
  file    = "m_h2_s4_20260505",
  backend = "cmdstanr")

# plus observer (full model)
m_h2_s5 <- brm(
  BodyMass ~ age_year + I(age_year^2) + RightTarsus + SexEstimate + summer +
    minutes_s + avg_invert + group_size +
    (1 | gr(animal, cov = Amat)) + (1 | BirdID) + (1 | BirthYear) +
    (1 | CatchYear) + (1 | Observer),
  data    = df_clean,
  data2   = list(Amat = Amat),
  family  = gaussian(),
  prior   = my_priors,
  warmup  = 1500,
  iter    = 11500,
  chains  = 4,
  cores   = 4,
  seed    = 123,
  file    = "m_brms_h2_1_20260503",
  backend = "cmdstanr")

#compare and plot the components
step_models <- list(
  "S1\n(+animal)"         = m_h2_s1,
  "S2\n(+BirdID)"         = m_h2_s2,
  "S3\n(+BirthYear)"      = m_h2_s3,
  "S4\n(+CatchYear)"      = m_h2_s4,
  "S5\n(+Observer)\n[full]" = m_h2_s5
)

vc_map <- tribble(~param,~label,
                  "sd_animal__Intercept",    "V[g]~(additive~genetic)",
                  "sd_BirdID__Intercept",    "V[pe]~(perm.~env.)",
                  "sd_BirthYear__Intercept", "V[by]~(birth~year)",
                  "sd_CatchYear__Intercept", "V[cy]~(catch~year)",
                  "sd_Observer__Intercept",  "V[obs]~(observer)",
                  "sigma",                   "V[e]~(residual)")

#fixed effect
get_vf <- function(model) {
  mu <- posterior_linpred(model, re_formula = NA)
  apply(mu, 1, var)
}

# Extract raw variances and per-draw h2
extract_step <- function(model, step_label) {
  Vf_vec  <- get_vf(model)
  raw     <- as_draws_df(model)
  present <- vc_map %>% filter(param %in% names(raw))
  
  vc_mat <- map_dfc(
    set_names(present$param, present$label),
    \(p) raw[[p]]^2
  ) %>%
    mutate(Vf = Vf_vec[seq_len(n())])
  
  vc_cols <- c(present$label, "Vf")
  
  vc_mat %>%
    mutate(
      Vp = rowSums(across(all_of(vc_cols))),
      h2 = `V[g]~(additive~genetic)` / Vp
    ) %>%
    select(all_of(vc_cols), Vp, h2) %>%
    mutate(step = step_label)
}

# Run (posterior_linpred)
step_draws <- imap_dfr(step_models, extract_step)

comp_levels <- c(
  "Vf",
  "V[e]~(residual)",
  "V[obs]~(observer)",
  "V[cy]~(catch~year)",
  "V[by]~(birth~year)",
  "V[pe]~(perm.~env.)",
  "V[g]~(additive~genetic)"
)

# Stacked bar: mean(Vi) / mean(Vp)
vc_cols_present <- intersect(comp_levels, names(step_draws))

step_summary <- step_draws %>%
  group_by(step) %>%
  summarise(across(all_of(c(vc_cols_present, "Vp")), mean), .groups = "drop") %>%
  mutate(across(all_of(vc_cols_present), \(x) x / Vp)) %>%
  select(step, all_of(vc_cols_present)) %>%
  pivot_longer(-step, names_to = "component", values_to = "proportion") %>%
  mutate(
    step      = factor(step, levels = names(step_models)),
    component = factor(component,
                       levels = intersect(comp_levels, unique(component)))
  )

# h2: mean_hdi() on per-draw h2 for error bar
h2_summary <- step_draws %>%
  group_by(step) %>%
  mean_hdi(h2) %>%
  rename(mean = h2, lo95 = .lower, hi95 = .upper) %>%
  select(step, mean, lo95, hi95) %>%
  mutate(step = factor(step, levels = names(step_models)))

vc_colors <- c(
  "V[g]~(additive~genetic)" = "#7B7B7B",
  "V[pe]~(perm.~env.)"      = "#D4B8E0",
  "V[by]~(birth~year)"      = "#5DA0D0",
  "V[cy]~(catch~year)"      = "#2B5FA5",
  "V[obs]~(observer)"       = "#3BAF8A",
  "V[e]~(residual)"         = "#CCCCCC",
  "Vf"                      = "#F4A261"
)

# h2 per step point
p_h2 <- ggplot(h2_summary, aes(x = step, y = mean)) +
  geom_pointrange(aes(ymin = lo95, ymax = hi95),
                  size = 0.6, linewidth = 0.8, color = "black") +
  geom_text(aes(label = round(mean, 3)),
            vjust = -2, size = 4, color = "black") +
  scale_y_continuous(
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.2)),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(x = NULL,
       y = expression(h^2 ~ "(posterior mean ± 95% HDI)")) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(colour = "black", size = 10))

# stacked bar
p_stack <- step_summary %>%
  ggplot(aes(x = step, y = proportion, fill = component)) +
  geom_col(position = "stack", width = 0.65, color = "white", linewidth = 0.3) +
  scale_fill_manual(
    values = vc_colors,
    name   = "Variance\ncomponent",
    labels = function(x) parse(text = x)
  ) +
  scale_y_continuous(
    expand = c(0, 0),
    labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = NULL,
    y     = "Proportion of phenotypic variance",
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.line       = element_line(linewidth = 0.5),
    axis.text       = element_text(colour = "black", size = 10),
    legend.key.size = unit(0.45, "cm")
  )

#plot

p_h2_stepwise <- p_h2 / p_stack +
  plot_annotation(tag_levels = "A")

ggsave("./figures/p_h2_stepwise_20260505.pdf",
       plot   = p_h2_stepwise,
       width  = 220,
       height = 200,
       units  = "mm",
       dpi    = 300)
