library(brms)
library(tidyverse)
library(tidybayes)
library(ggdist)
library(patchwork)
library(lme4)

# Load model
m_h2_age <- readRDS("m_brms_h2_age20260507.rds")
#m_h2_age <- readRDS("m_brms_h2_age_genomic20260507.rds")

age_classes <- c("Juv", "Young", "Old")

# Variance component map
vc_map_age <- tribble(~param_prefix, ~label,
                      "sd_animal__AgeClass",    "V[a]~(additive~genetic)",
                      "sd_BirdID__AgeClass",    "V[pe]~(perm.~env.)",
                      "sd_BirthYear__AgeClass", "V[by]~(birth~year)",
                      "sd_CatchYear__AgeClass", "V[cy]~(catch~year)",
                      "sd_Observer__AgeClass",  "V[obs]~(observer)"
)

# Vf
get_vf_safe <- function(model, draw_ids) {
  all_draws <- as_draws_df(model)
  b_cols    <- grep("^b_(?!sigma)", colnames(all_draws), perl = TRUE, value = TRUE)
  B         <- as.matrix(as.data.frame(all_draws[draw_ids, b_cols]))
  fe_formula <- nobars(model$formula$formula)
  X          <- model.matrix(fe_formula, data = model$data)
  apply(B, 1, \(b) var(as.numeric(X %*% b)))
}

# Extractor
age_means <- df_clean_age %>%
  group_by(AgeClass) %>%
  summarise(mean_mass = mean(BodyMass, na.rm = TRUE))

extract_age <- function(model, age_class, ndraws = 10000, seed = 42) {
  
  message("Extracting age class: ", age_class)
  
  # Age-class-specific mean body mass for evolvability
  mean_mass <- df_clean_age %>%
    filter(AgeClass == age_class) %>%
    pull(BodyMass) %>%
    mean(na.rm = TRUE)
  
  all_draws <- as_draws_df(model)
  n_total   <- nrow(all_draws)
  
  set.seed(seed)
  draw_ids <- sample(seq_len(n_total), size = ndraws)
  
  Vf_vec <- get_vf_safe(model, draw_ids)
  raw    <- as.data.frame(all_draws[draw_ids, ])
  
  present <- vc_map_age %>%
    filter(map_lgl(param_prefix,
                   \(p) paste0(p, age_class) %in% names(raw)))
  
  vc_list <- setNames(
    lapply(paste0(present$param_prefix, age_class), \(p) raw[[p]]^2),
    present$label
  )
  
  sigma_col <- paste0("b_sigma_AgeClass", age_class)
  vc_list[["V[e]~(residual)"]] <- exp(raw[[sigma_col]])^2
  vc_list[["Vf"]]              <- Vf_vec
  
  df         <- as.data.frame(vc_list, check.names = FALSE)
  vc_cols    <- names(vc_list)
  df$Vp      <- rowSums(df[, vc_cols, drop = FALSE])
  df$h2_with    <- df[["V[a]~(additive~genetic)"]] / df$Vp
  df$h2_without <- df[["V[a]~(additive~genetic)"]] / (df$Vp - df[["V[obs]~(observer)"]])
  df$rep        <- (df[["V[a]~(additive~genetic)"]] + df[["V[pe]~(perm.~env.)"]]) /df$Vp
  df$CVa2       <- 100 * (sqrt(df[["V[a]~(additive~genetic)"]]) / mean_mass)
  df$AgeClass   <- age_class
  df}

# Run for all age classes
age_draws <- map_dfr(age_classes, \(ac) extract_age(m_h2_age, ac))
age_draws_tabel <- age_draws %>%
  group_by(AgeClass) %>%
  summarise(
    across(
      c(`V[a]~(additive~genetic)`,
        `V[pe]~(perm.~env.)`,
        `V[by]~(birth~year)`,
        `V[cy]~(catch~year)`,
        `V[obs]~(observer)`,
        `V[e]~(residual)`,
        Vf, Vp,
        h2_with, h2_without,
        rep, CVa2),
      list(mean  = mean,
           lower = ~quantile(.x, 0.025),
           upper = ~quantile(.x, 0.975)),
      .names = "{.col}__{.fn}"
    ),
    .groups = "drop"
  )

View(age_draws_tabel)

# Component order
comp_levels_age <- c(
  "Vf",
  "V[e]~(residual)",
  "V[obs]~(observer)",
  "V[cy]~(catch~year)",
  "V[by]~(birth~year)",
  "V[pe]~(perm.~env.)",
  "V[a]~(additive~genetic)"
)

# Stacked bar summary
vc_cols_age <- intersect(comp_levels_age, names(age_draws))

age_summary <- age_draws %>%
  group_by(AgeClass) %>%
  summarise(across(all_of(c(vc_cols_age, "Vp")), mean), .groups = "drop") %>%
  mutate(across(all_of(vc_cols_age), \(x) x / Vp)) %>%
  pivot_longer(all_of(vc_cols_age), names_to = "component", values_to = "proportion") %>%
  mutate(
    AgeClass  = factor(AgeClass,  levels = age_classes),
    component = factor(component, levels = intersect(comp_levels_age, unique(component))))

# h2 summary: per-draw HDI
h2_summary_age <- age_draws %>%
  group_by(AgeClass) %>%
  mean_hdi(h2) %>%
  rename(mean = h2, lo95 = .lower, hi95 = .upper) %>%
  mutate(AgeClass = factor(AgeClass, levels = age_classes))

# Colours and labels
comp_colors_age <- c(
  "V[a]~(additive~genetic)" = "#1a0a3d",
  "V[pe]~(perm.~env.)"      = "#2d1b69",
  "V[by]~(birth~year)"      = "#5c2d91",
  "V[cy]~(catch~year)"      = "#8b44b8",
  "V[obs]~(observer)"       = "#e07070",
  "V[e]~(residual)"         = "#f0a060",
  "Vf"                      = "#f5e090"
)

age_labels <- c("Juv" = "Juvenile", "Young" = "Young", "Old" = "Old")

# Stacked bar
p_stack_age <- ggplot(age_summary,
                      aes(x = AgeClass, y = proportion, fill = component)) +
  geom_bar(stat = "identity", width = 0.65) +
  scale_fill_manual(
    values = comp_colors_age,
    labels = \(x) parse(text = x),
    name   = "Variance\ncomponents"
  ) +
  scale_x_discrete(labels = age_labels) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                     breaks = seq(0, 1, 0.25)) +
  labs(x = NULL,
       y = "Proportion of variance explained") +
  theme_bw(base_size = 15, base_family = "Helvetica") +
  theme(
    panel.grid        = element_blank(),
    plot.title        = element_text(hjust = 0.5),
    legend.text.align = 0)

# Plot 2: h2 posterior with HDI
p_h2_age <- ggplot() +
  stat_halfeye(
    data = age_draws %>% mutate(AgeClass = factor(AgeClass, levels = age_classes)),
    aes(x = AgeClass, y = h2),
    .width         = c(0.66, 0.95),
    point_interval = median_hdi,
    slab_alpha     = 0.6,
    normalize      = "panels"
  ) +
  #geom_pointrange(
  #  data = h2_summary_age,
  # aes(x = AgeClass, y = mean, ymin = lo95, ymax = hi95, colour = AgeClass),
  #  size        = 0.8,
  #  linewidth   = 1.1,
  #  show.legend = FALSE
  #) +
  #scale_fill_manual(
  #  values = c("Juv" = "#1a0a3d", "Young" = "#8b44b8", "Old" = "#e07070"),
  #  guide  = "none"
  #) +
  #scale_colour_manual(
  #  values = c("Juv" = "#1a0a3d", "Young" = "#8b44b8", "Old" = "#e07070")
  #) +
  scale_x_discrete(labels = age_labels) +
  scale_y_continuous(limits = c(0, 0.3), breaks = seq(0, 1, 0.1)) +
  labs(
    x = "Age class",
    y = expression(italic(h)^2)) +
  theme_bw(base_size = 15, base_family = "Helvetica") +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position    = "none")

# Combine and save
p_stack_age + p_h2_age +
  plot_layout(widths = c(1.2, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("h2_age_classes_pedigree.pdf", width = 11, height = 5, dpi = 300)
ggsave("h2_age_classes_pedigree.png", width = 11, height = 5, dpi = 300)
#ggsave("h2_age_classes_genomic.pdf", width = 11, height = 5, dpi = 300)
#ggsave("h2_age_classes_genomic.png", width = 11, height = 5, dpi = 300)
