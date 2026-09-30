#load packages
library(dplyr)
library(ggplot2)
library(patchwork)

#load dataset
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")

dat <- df_clean %>%
  group_by(BirdID) %>%
  mutate(age_mean = mean(age_year, na.rm = TRUE),age_dev  = age_year - age_mean) %>%
  ungroup()

#Plot the distribution of age and bodymass
p_age <- ggplot(dat, aes(x = age_year)) +
  geom_histogram(binwidth = 1, fill = "grey70", color = "white") +
  stat_bin(binwidth = 1,geom = "text",
           aes(label = after_stat(count)),
           vjust = -0.3,
           size = 3.5) +
  labs(x = "Age (yr)",
       y = "Number of observations",
       title = "Age distribution across observations") +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 16, colour = "black"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)) +ylim(0, 1150)
p_mass <- ggplot(dat, aes(x = BodyMass)) +
  geom_histogram(binwidth = 1, fill = "grey70", color = "white") +
  stat_bin(binwidth = 1,geom = "text",
           aes(label = after_stat(count)),
           vjust = -0.3,
           size = 3.5) +
  labs(x = "Body mass (g)",
       y = "Number of observations",
       title = "Body mass distribution across observations") +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 16, colour = "black"),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16)) +ylim(0, 1150)

p_age_mass <- (p_age | p_mass) + plot_annotation(tag_levels = "A")
ggsave("./figures/p_age_mass20260502.pdf", plot = p_age_mass, width = 300, height = 150, units = "mm", dpi = 300)

#Plot the distribution of bird records
counts_per_bird <- dat %>%
  group_by(BirdID) %>%
  summarise(n_measurements = n()) %>%
  ungroup()

#Plot the distribution
ggplot(counts_per_bird, aes(x = n_measurements)) +
  geom_histogram(binwidth = 1, fill = "grey70", color = "white") +
  labs(x = "Number of body mass measurements per bird",
       y = "Number of individuals",
       title = "Distribution of body mass measurements per individual") +
  theme_classic()+
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_blank())

accum_data <- counts_per_bird %>%
  count(n_measurements) %>%
  arrange(n_measurements) %>%
  mutate(cumulative_individuals = cumsum(n),
         total_individuals = sum(n),
         percent = cumulative_individuals / total_individuals * 100)

p_measurements <- ggplot(accum_data, aes(x = n_measurements)) +
  geom_col(aes(y = n), fill = "grey70", color = "white") +
  geom_text(aes(y = n, label = n),vjust = -0.3, size = 3.5) +
  geom_line(aes(y = percent * max(n) / 100), color = "black", size = 1,alpha=0.4) +
  geom_point(aes(y = percent * max(n) / 100), color = "black", size = 2,alpha=0.4) +
  scale_y_continuous(name = "Number of individuals",
                     sec.axis = sec_axis(~ . * 100 / max(accum_data$n), name = "Cumulative % of individuals")) +
  labs(x = "Number of body mass measurements per bird",title = "Histogram and Accumulation of Body Mass Measurements") +
  theme_classic()+
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_blank())
ggsave("./figures/p_measurements20260502.pdf",plot = p_measurements,width = 200, height = 150, units = "mm",dpi = 300)

# Compute stats
ci_age  <- t.test(dat$age_year)$conf.int
ci_mass <- t.test(dat$BodyMass)$conf.int

age_lbl  <- sprintf("Mean = %.2f yr (95%% CI: %.2f\u2013%.2f)",
                    mean(dat$age_year), ci_age[1], ci_age[2])
mass_lbl <- sprintf("Mean = %.2f g (95%% CI: %.2f\u2013%.2f)",
                    mean(dat$BodyMass), ci_mass[1], ci_mass[2])

#Sex differences in body mass
sex_stats <- dat %>%
  mutate(Sex = ifelse(SexEstimate == 1, "Male", "Female")) %>%
  group_by(Sex) %>%
  summarise(
    mean  = mean(BodyMass),
    lower = t.test(BodyMass)$conf.int[1],
    upper = t.test(BodyMass)$conf.int[2],
    .groups = "drop"
  )

mass_sex_lbl <- sprintf(
  "Female: %.2f g (95%% CI: %.2f\u2013%.2f) \n Male: %.2f g (95%% CI: %.2f\u2013%.2f)",
  sex_stats$mean[1], sex_stats$lower[1], sex_stats$upper[1],
  sex_stats$mean[2], sex_stats$lower[2], sex_stats$upper[2]
)

sex_colors <- c("Female" = "#C47E5A", "Male" = "#4D9AA8")

p_mass_sex <- dat %>%
  mutate(Sex = ifelse(SexEstimate == 1, "Male", "Female")) %>%
  ggplot(aes(x = BodyMass, fill = Sex)) +
  geom_histogram(binwidth = 1, color = "white",
                 position = "identity", alpha = 0.65) +
  geom_vline(data = sex_stats,
             aes(xintercept = mean, color = Sex),
             linetype = "dashed", linewidth = 0.8) +
  scale_fill_manual(values  = sex_colors) +
  scale_color_manual(values = sex_colors, guide = "none") +
  labs(x        = "Body mass (g)",
       y        = "Number of observations",
       title    = "Body mass distribution by sex",
       fill     = NULL) +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(legend.position      = c(0.88, 0.88),
        legend.text          = element_text(size = 15),
        legend.background    = element_blank(),
        legend.key.size      = unit(0.5, "cm"))

ggsave("./figures/p_mass_sex20260502.pdf",plot = p_mass_sex, width = 300, height = 150, units = "mm", dpi = 300)

p_mass


ggsave("./figures/p_sex_diff_mass20260502.pdf",plot = p_mass_sex,width = 200, height = 150, units = "mm",dpi = 300)

#Check for sex biases in tarsus control
plot(dat$RightTarsus, dat$BodyMass, col = dat$SexEstimate + 1, pch = 1)
legend("topleft", legend = c("Female", "Male"),col = c(1, 2), pch = 16)


