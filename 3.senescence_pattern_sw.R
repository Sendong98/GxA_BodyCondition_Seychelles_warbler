library(dplyr)
library(lmerTest)
library(DHARMa)
library(ggplot2)
library(ggeffects)
library(sjPlot)
library(car)
library(lubridate)
library(patchwork)

#load dataset
setwd("~/GxA_bodycondition_brms")
df_clean <- read.csv("./datasets/dataset_mass_20260503.csv")

dat <- df_clean %>%
  group_by(BirdID) %>%
  mutate(age_mean = mean(age_year, na.rm = TRUE), # age_i
         age_dev  = age_year - age_mean, # within-individual delta age
         age2_dev = age_year^2 - age_mean^2 # quadratic within  
  ) %>%ungroup()

age_at_last_seen <- readxl::read_xlsx("./datasets/sys_LastSeenwithGenParentandDeathRecord.xlsx")
dat <- dplyr::left_join(dat, age_at_last_seen, by = "BirdID")
dat$age_last_seen <- time_length(interval(dat$BirthDate, dat$LastSeen), "years")

#scale for modelling
dat$minutes_s <- scale(dat$minutes)

#basic model1 (log didnt improve much)
m1 <- lmer(formula = BodyMass ~ age_year + I(age_year^2) + SexEstimate + RightTarsus + summer + minutes_s + avg_invert + group_size +
             (1|BirthYear) +(1|CatchYear) + (1|Observer) + (1|BirdID), REML = F ,data = dat)
simulationOutput <- simulateResiduals(fittedModel = m1, plot = T)
qqnorm(resid(m1))
qqline(resid(m1))
vif(m1)
summary(m1)
tab_model(m1)

#plot ageing pattern
dat_m1 <- predict_response(m1, terms = c("age_year [all]"))
plot_data <- as.data.frame(dat_m1)

p_m1 <- ggplot(plot_data, aes(x = x, y = predicted)) +
  geom_point(data = dat, aes(x = age_year, y = BodyMass),
             alpha = 0.2, size = 2, colour = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              fill = "#808000", alpha = 0.2) +
  geom_line(colour = "#808000", linewidth = 1) +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(axis.line         = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks        = element_line(linewidth = 0.4, colour = "black"),
        axis.text         = element_text(size = 16, colour = "black"),
        axis.title        = element_text(size = 16, colour = "black"),
        legend.title      = element_text(size = 16, face = "bold"),
        legend.text       = element_text(size = 16),
        legend.key.size   = unit(0.45, "cm"),
        legend.background = element_blank(),
        plot.title        = element_blank(),
        panel.border      = element_blank()) +
  xlab("Age (yr)") + ylab("Body Condition")

ggsave("./figures/p_mass_plot_all_ageing_20260502.pdf", plot = p_m1, width = 200, height = 150, units = "mm", dpi = 300)

#post-peak ageing
fix_m1 <- fixef(m1)
a <- fix_m1["I(age_year^2)"]
b <- fix_m1["age_year"]

# calculate peak
peak_age <- -b / (2 * a)

post_peak <- dat %>% mutate(age_cat = factor(age_year >= peak_age, 
                                             labels = c("Pre", "Post")))

mass_age_mod_pre <- lmer(BodyMass ~ age_year + RightTarsus + SexEstimate  + summer + minutes_s + avg_invert + group_size + (1|BirthYear) +(1|CatchYear) + (1|Observer) + (1|BirdID),
                         data = post_peak[post_peak$age_cat=="Pre",],control = lmerControl(optimizer = "bobyqa",
                                                                                           optCtrl   = list(maxfun = 2e5)))
allFit(mass_age_mod_pre)
simulationOutput <- simulateResiduals(fittedModel = mass_age_mod_pre, plot = T)
qqnorm(resid(mass_age_mod_pre))
qqline(resid(mass_age_mod_pre))
vif(mass_age_mod_pre)
summary(mass_age_mod_pre)
tab_model(mass_age_mod_pre)

mass_age_mod_post <- lmer(BodyMass ~ age_year + RightTarsus + SexEstimate  + summer + minutes_s + avg_invert + group_size + (1|BirthYear) +(1|CatchYear) + (1|Observer) + (1|BirdID), 
                          data = post_peak[post_peak$age_cat=="Post",])
simulationOutput <- simulateResiduals(fittedModel = mass_age_mod_post, plot = T)
qqnorm(resid(mass_age_mod_post))
qqline(resid(mass_age_mod_post))
vif(mass_age_mod_post)
summary(mass_age_mod_post)
tab_model(mass_age_mod_post)

# age_year p-values from both models
coef_pre  <- summary(mass_age_mod_pre)$coefficients
coef_post <- summary(mass_age_mod_post)$coefficients

p_pre  <- coef_pre["age_year", "Pr(>|t|)"]
b_pre  <- coef_pre["age_year", "Estimate"]
p_post <- coef_post["age_year", "Pr(>|t|)"]
b_post <- coef_post["age_year", "Estimate"]

cat(sprintf("Pre-peak:  β = %.3f, p = %.4f\n", b_pre,  p_pre))
cat(sprintf("Post-peak: β = %.3f, p = %.4f\n", b_post, p_post))
cat(sprintf("Peak age: %.2f years\n", peak_age))

mass_prediction_pre <- ggpredict(mass_age_mod_pre, terms = c("age_year")) %>% 
  as_tibble() %>% mutate(age = x ,AgeClass="PrePeak") %>% filter(age<=peak_age)

mass_prediction_post <- ggpredict(mass_age_mod_post, terms = c("age_year")) %>%
  as_tibble() %>% mutate(age = x, AgeClass="PostPeak") %>% filter(age>peak_age)

prediction_all <- rbind(mass_prediction_pre,mass_prediction_post)

mass_plot <- ggplot(prediction_all, aes(age, predicted)) +
  geom_point(data = post_peak,aes(x = age_year, y = BodyMass),alpha = 0.2,size = 2) +
  geom_ribbon(data = mass_prediction_pre,aes(x = age, ymin = conf.low, ymax = conf.high),alpha = 0.2,fill = "#c8c87a") +
  geom_line(data = mass_prediction_pre,aes(x = age, y = predicted),linewidth = 1,color="#c8c87a")  +
  geom_ribbon(data = mass_prediction_post,aes(x = age, ymin = conf.low, ymax = conf.high),alpha = 0.2,fill = "#3d3d06") +
  geom_line(data = mass_prediction_post,aes(x = age, y = predicted),linewidth = 1,color="#3d3d06")  +
  geom_vline(xintercept = peak_age, linetype = "dashed") +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(axis.line         = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks        = element_line(linewidth = 0.4, colour = "black"),
        axis.text         = element_text(size = 16, colour = "black"),
        axis.title        = element_text(size = 16, colour = "black"),
        legend.title      = element_text(size = 16, face = "bold"),
        legend.text       = element_text(size = 16),
        legend.key.size   = unit(0.45, "cm"),
        legend.background = element_blank(),
        plot.title        = element_blank(),
        panel.border      = element_blank()
  ) + xlab("Age (yr)") + ylab("Body Condition")

format_p <- function(p) {
  if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
}

pre_label  <- sprintf("β = %.3f, %s",  b_pre,  format_p(p_pre))
post_label <- sprintf("β = %.3f, %s", b_post, format_p(p_post))
peak_label <- sprintf("Peak age: %.1f yr", peak_age)

# x positions for annotations — just inside each segment
x_pre  <- peak_age / 2
x_post <- peak_age + (max(post_peak$age_year, na.rm = TRUE) - peak_age) / 2

mass_plot_post_peak <- mass_plot + annotate("text", x = x_pre,  y = Inf, label = pre_label,
                                            hjust = 0.5, vjust = 1.5, size = 5, colour = "black") +
  annotate("text", x = x_post, y = Inf, label = post_label,
           hjust = 0.5, vjust = 1.5, size = 5, colour = "black") +
  annotate("text", x = peak_age, y = -Inf, label = peak_label,
           hjust = -0.08, vjust = -0.6, size = 5, fontface = "italic")
ggsave("./figures/p_mass_plot_post_peak_20260502.pdf", plot = mass_plot_post_peak, width = 200, height = 150, units = "mm", dpi = 300,device = cairo_pdf)

#combine together
overall_ageing <- (p_m1 | mass_plot_post_peak)+ plot_annotation(tag_levels = "A")
ggsave("./figures/p_overall_ageing_20260502.pdf", plot = overall_ageing, width = 300, height = 150, units = "mm", dpi = 300,device = cairo_pdf)

#test for minutes^2, not significate
m2 <- lmer(formula = BodyMass ~ age_year + I(age_year^2) + SexEstimate + RightTarsus + summer + minutes_s + I(minutes_s^2)+ avg_invert + group_size +
             (1|BirthYear) +(1|CatchYear) + (1|Observer) + (1|BirdID), REML = F ,data = dat)
simulationOutput <- simulateResiduals(fittedModel = m2, plot = T)
qqnorm(resid(m2))
qqline(resid(m2))
vif(m2)
summary(m2)
tab_model(m2)

#Sex interaction
#Add sex*age_year interaction
m3 <- lmer(formula = BodyMass ~ age_year*SexEstimate + I(age_year^2)*SexEstimate + RightTarsus + summer + minutes_s + I(minutes_s^2)+ avg_invert + group_size +
             (1|BirthYear) +(1|CatchYear) + (1|Observer) + (1|BirdID), REML = F ,data = dat)
simulationOutput <- simulateResiduals(fittedModel = m3, plot = T)
qqnorm(resid(m3))
qqline(resid(m3))
vif(m3)
summary(m3)
tab_model(m3)

m3.1 <- lmer(formula = BodyMass ~ age_year*SexEstimate + RightTarsus + summer + minutes_s + avg_invert + group_size +
               (1|BirthYear) +(1|CatchYear) + (1|Observer) + (1|BirdID), REML = F ,data = dat)
simulationOutput <- simulateResiduals(fittedModel = m3, plot = T)
qqnorm(resid(m3.1))
qqline(resid(m3.1))
vif(m3.1)
summary(m3.1)
tab_model(m3.1)

#test selective disappearance / within–between individual ageing process

# Mean-centering within and among 
m4 <- lmer(BodyMass ~ age_dev + age2_dev + age_mean + SexEstimate + RightTarsus + summer + minutes_s + avg_invert + group_size +
             (1|BirthYear) + (1|CatchYear) + (1|Observer) + (1|BirdID),data = dat, REML = FALSE)
simulationOutput <- simulateResiduals(fittedModel = m4, plot = T)
qqnorm(resid(m4))
qqline(resid(m4))
vif(m4)
summary(m4)
tab_model(m4)

#additive last seen
m5 <- lmer(BodyMass ~ age_year + I(age_year^2) + age_last_seen + SexEstimate + RightTarsus + summer + minutes_s + avg_invert + group_size +
             (1|BirthYear) + (1|CatchYear) + (1|Observer) + (1|BirdID), data = dat, REML = FALSE)
simulationOutput <- simulateResiduals(fittedModel = m5, plot = T)
qqnorm(resid(m5))
qqline(resid(m5))
vif(m5)
summary(m5)
tab_model(m5)

# age at last seen as interaction
m6 <- lmer(BodyMass ~ age_year * age_last_seen + I(age_year^2) * age_last_seen + SexEstimate + RightTarsus + summer + minutes_s + avg_invert + group_size +
             (1|BirthYear) + (1|CatchYear) + (1|Observer) + (1|BirdID), data = dat, REML = FALSE)
allFit(m6)
simulationOutput <- simulateResiduals(fittedModel = m6, plot = T)
qqnorm(resid(m6))
qqline(resid(m6))
vif(m6)
summary(m6)
tab_model(m6)

AIC(m1, m4, m5, m6) #lowest AIC:m6 
anova(m1, m4)
anova(m4, m6)
anova(m5, m6)

#Plot ageing pattern
#Plot within individual pattern
b <- fixef(m4)
# Fixed-effects fitted values using each obs's actual covariate values
X_actual <- model.matrix(~ age_dev + age2_dev + age_mean + SexEstimate + RightTarsus +
                           summer + minutes_s + avg_invert + group_size, data = dat)
fitted_fixef_actual <- as.numeric(X_actual %*% b)

# Fixed-effects fitted values using each obs's actual age_dev/age2_dev
# but all other covariates held at the same means used in newdat
X_age_at_mean <- model.matrix(~ age_dev + age2_dev + age_mean + SexEstimate + RightTarsus +
                                summer + minutes_s + avg_invert + group_size,
                              data = data.frame(
                                age_dev      = dat$age_dev,
                                age2_dev     = dat$age2_dev,
                                age_mean     = mean(dat$age_year),
                                SexEstimate  = 0,
                                RightTarsus  = mean(dat$RightTarsus),
                                summer       = 0,
                                minutes_s    = mean(dat$minutes_s),
                                avg_invert   = mean(dat$avg_invert),
                                group_size   = mean(dat$group_size)))

fitted_age_at_mean <- as.numeric(X_age_at_mean %*% b)

#remove all non-age fixed-effects variation, keep age contribution and noise
dat$partial_resid <- dat$BodyMass - fitted_fixef_actual + fitted_age_at_mean

# Build prediction grid for m4: vary age_dev across observed range,
# hold all other covariates at their means (Sex is female and at winter)
# age2_dev = age_year^2 - age_mean^2 = age_dev^2 + 2*age_dev*age_mean,
# so the cross term must be included to match the model's parameterisation

age_dev_seq  <- seq(min(dat$age_dev), max(dat$age_dev), length.out = 100)
age_mean_ref <- mean(dat$age_year)
newdat <- data.frame(
  age_dev     = age_dev_seq,
  age2_dev    = age_dev_seq^2 + 2 * age_dev_seq * age_mean_ref,
  age_mean    = age_mean_ref,
  SexEstimate = 0,
  RightTarsus = mean(dat$RightTarsus, na.rm = TRUE),
  summer      = 0,
  minutes_s   = 0,
  avg_invert  = mean(dat$avg_invert, na.rm = TRUE),
  group_size  = mean(dat$group_size, na.rm = TRUE)
)

mm_newdat    <- model.matrix(~ age_dev + age2_dev + age_mean + SexEstimate +
                               RightTarsus + summer + minutes_s + avg_invert +
                               group_size, data = newdat)
se_newdat    <- sqrt(diag(mm_newdat %*% vcov(m4) %*% t(mm_newdat)))
newdat$fit   <- as.numeric(mm_newdat %*% b)
newdat$lower <- newdat$fit - 1.96 * se_newdat
newdat$upper <- newdat$fit + 1.96 * se_newdat

head(newdat[, c("age_dev", "age2_dev", "fit", "lower", "upper")])

p_m4 <- ggplot(newdat, aes(x = age_dev, y = fit)) +
  geom_point(data = dat, aes(x = age_dev, y = partial_resid), alpha = 0.15, size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2,fill="#808000") +
  geom_line(linewidth = 1.2,colour = "#808000",alpha = 0.8) +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(axis.line = element_line(linewidth = 0.4), axis.ticks = element_line(linewidth = 0.4),
        axis.text = element_text(size = 16), axis.title = element_text(size = 16),
        legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.45, "cm"),
        legend.background = element_blank())+
  labs(x = "Within-individual age deviation", y = "Predicted body condition (partial residuals)")+
  coord_cartesian(ylim = c(9, 20))

p_m4

#Plot selective disappearance
pred_m6 <- ggpredict(m6, terms = c("age_year[all]", "age_last_seen [quart]"))
pred_df <- as.data.frame(pred_m6)
pred_df$group <- factor(pred_df$group,labels = c("Q1 (low lifespan)","Q2","Q3","Q4","Q5 (high lifespan)"))

olive_pal <- colorRampPalette(c("#c8c87a", "#808000", "#3d3d06"))(5)

p_m6 <- ggplot(pred_df, aes(x = x, y = predicted, colour = group, fill = group)) +
  geom_point(data = dat, aes(x = age_year, y = BodyMass),
             alpha = 0.05, size = 0.5,
             inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1) +
  labs(
    x = "Age (yr)",
    y = "Predicted body condition",
    colour = "Lifespan class",
    fill = "Lifespan class"
  ) +
  theme_classic(base_size = 16, base_family = "Helvetica") +
  theme(
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4),
    axis.text    = element_text(size = 16),
    axis.title   = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text  = element_text(size = 16),
    legend.key.size = unit(0.45, "cm"),
    legend.background = element_blank(),
    legend.position = c(0.3, 0.2)
  ) +
  scale_colour_manual(values = olive_pal) +
  scale_fill_manual(values = olive_pal)

p_m6

p_ageing_pattern <- (p_m4 | p_m6) + plot_annotation(tag_levels = "A")
p_ageing_pattern

ggsave("./figures/p_ageing_pattern_20260502.pdf", plot = p_ageing_pattern, width = 300, height = 150, units = "mm", dpi = 300)

