####
##
##  Santtu Pentikäinen (Helsinki Urban Rat Project) 2026
##
##  Code for Helsinki rat dynamics study
##
##  bayesian model, diagnostics and plots
##
####


# require packages (and install if necessary)
#####
pkg_list <- c(
  "tidyverse",
  "brms",
  "cmdstanr",
  "bayesplot",
  "patchwork",
  "DHARMa",
  "loo",
  "performance"
)

# install
missing_pkgs <- pkg_list[!(pkg_list %in% installed.packages()[, "Package"])]
if (length(missing_pkgs)) install.packages(missing_pkgs)
rm(missing_pkgs)

# require
lapply(pkg_list, library, character.only = TRUE)

#####

# data and model
#####

# download data from Zenodo
path <- "path-to-downloaded-files/"
helsinkirats_data <- readRDS(paste0(path,"helsinkirats_data.rds"))


## Option 1: load the fitted model (recommended)

# download helsinkirats_model.rds from Zenodo
helsinkirats_model <- readRDS(paste0(path,"helsinkirats_model.rds"))


## Option 2: reproduce the model from the published data
## This may take several hours

theformula <- bf(rottasaalis ~ log_effort * asemamalli_faktori + 
                   kohdetyyppi_faktori + 
                   elintarviketoimijat_log_skaalattu + 
                   stat_sosioekon_indeksi_skaalattu + 
                   vih_osuus_rak_log_skaalattu * vuodenaika + 
                   interval_tempmean_allweek_skaalattu + 
                   interval_tempmean_daysbelow10_6months_log_skaalattu + 
                   interval_mpi_90days_mean_skaalattu + 
                   (1 | vuosi_faktori) + (1 | kohde_faktori) + 
                   asukastiheys_km2_kokoala_log_skaalattu + 
                   edellinen_cpue_log_skaalattu + 
                   edellinen_cpue_log_skaalattu:edellinen_effort_log_skaalattu)


# refit the model by uncommenting:
# helsinkirats_model_reproduced <- brm(
#   formula = theformula,
#   data = helsinkirats_data,
#   family = negbinomial(),
#   chains = 3,
#   cores = 3,
#   iter = 4000,
#   warmup = 2000,
#   threads = threading(2),
#   backend = "cmdstanr",
#   save_pars = save_pars(all = TRUE),
#   control = list(adapt_delta = 0.95)
# )
# # continue the analysis using the refitted model, uncomment:
# helsinkirats_model <- helsinkirats_model_reproduced

#####

# diagnostics
#####

# R2-values
the_R2cond <- bayes_R2(helsinkirats_model)
the_R2marg <- bayes_R2(helsinkirats_model, re_formula = NA)

# LOO & Moment Match
if (!"loo" %in% names(helsinkirats_model$criteria)) {
  helsinkirats_model <- add_criterion(helsinkirats_model, 
                                      "loo", moment_match = TRUE)
}

# DHARMa-residuals
the_res <- createDHARMa(
  simulatedResponse = posterior_predict(helsinkirats_model, 
                                        ndraws = 500) %>% t(),
  observedResponse = helsinkirats_model$data$rottasaalis,
  fittedPredictedResponse = apply(posterior_epred(helsinkirats_model, ndraws = 500), 
                                  2, mean),
  integerResponse = TRUE)


# diagnostic summaries
###
##
# model summary
print(summary(helsinkirats_model))
# loo-cv
print(helsinkirats_model$criteria$loo)
# bayesian R2
round(mean(the_R2marg), 3)
round(mean(the_R2cond), 3)
# dharma residuals
print(testDispersion(the_res, plot = FALSE))
# outliers
print(testOutliers(the_res, plot = FALSE))
# zeroinf
print(testZeroInflation(the_res, plot = FALSE))
# temporal autocorrelation
print(testTemporalAutocorrelation(the_res, time = helsinkirats_model$data$pvm, plot = FALSE))
###

# diagnostic plots
###
##
#

# convergence of chains
plot(helsinkirats_model, ask = FALSE)

# posterior predictive check
pp_check(helsinkirats_model, ndraws = 100) + 
  scale_x_continuous(trans = "log1p") 

# zeros
prop_zero <- function(x) mean(x == 0)
pp_check(helsinkirats_model, 
         type = "stat", 
         stat = "prop_zero")

# dharma residuals
plot(the_res)

# acf
res_map <- aggregate(the_res$scaledResiduals, 
                     list(helsinkirats_data$pvm), mean)
acf(res_map$x)

# effects
plot(conditional_effects(helsinkirats_model), ask = FALSE)
###
#####

# table1
#####
# include rhat and ess
tab1 <- summary(helsinkirats_model)$fixed |>
  as.data.frame() |>
  tibble::rownames_to_column("Variable")
print(tab1)

# write.csv(tab1,
#           "table1_model_coefficients.csv",
#           row.names = FALSE)

#####

# plot weather and socioeconomics (fig2)
#####

# Recover scaling parameters for back-transformation
get_scaling_params <- function(scaled_val, 
                               orig_val, 
                               is_log = FALSE) {
  if (is_log) {
    x_input <- log(abs(orig_val) + 1)
  } else {
    x_input <- abs(orig_val)
  }
  fit <- lm(scaled_val ~ x_input)
  slope <- as.numeric(coef(fit)[2])
  intercept <- as.numeric(coef(fit)[1])
  return(list(
    sd = 1 / slope,
    mu = -intercept / slope
  ))
}

p_temp   <- get_scaling_params(helsinkirats_data$interval_tempmean_allweek_skaalattu, 
                               helsinkirats_data$original_interval_tempmean_allweek)
p_winter <- get_scaling_params(helsinkirats_data$interval_tempmean_daysbelow10_6months_log_skaalattu, 
                               helsinkirats_data$original_interval_tempmean_daysbelow10_6months, 
                               is_log = TRUE)
p_mpi    <- get_scaling_params(helsinkirats_data$interval_mpi_90days_mean_skaalattu, 
                               helsinkirats_data$original_interval_mpi_90days_mean)
p_socio  <- get_scaling_params(helsinkirats_data$stat_sosioekon_indeksi_skaalattu, 
                               helsinkirats_data$original_stat_sosioekon_indeksi)

weather_vars <- c("interval_tempmean_allweek_skaalattu", 
                  "interval_tempmean_daysbelow10_6months_log_skaalattu", 
                  "interval_mpi_90days_mean_skaalattu")
ce_weather <- conditional_effects(helsinkirats_model, 
                                  effects = weather_vars, 
                                  re_formula = NA)
ce_socio <- conditional_effects(helsinkirats_model, 
                                effects = "stat_sosioekon_indeksi_skaalattu", 
                                re_formula = NA)

df_A <- ce_weather[["interval_tempmean_allweek_skaalattu"]] %>%
  transmute(
    x_orig = effect1__ * p_temp$sd + p_temp$mu,
    estimate_rel = (estimate__ / mean(estimate__)) * 100,
    lower_rel = (lower__ / mean(estimate__)) * 100,
    upper_rel = (upper__ / mean(estimate__)) * 100
  )

df_B <- ce_weather[["interval_tempmean_daysbelow10_6months_log_skaalattu"]] %>%
  transmute(
    x_orig = exp(effect1__ * p_winter$sd + p_winter$mu) - 1,
    estimate_rel = (estimate__ / mean(estimate__)) * 100,
    lower_rel = (lower__ / mean(estimate__)) * 100,
    upper_rel = (upper__ / mean(estimate__)) * 100
  )

df_C <- ce_weather[["interval_mpi_90days_mean_skaalattu"]] %>%
  transmute(
    x_orig = effect1__ * p_mpi$sd + p_mpi$mu,
    estimate_rel = (estimate__ / mean(estimate__)) * 100,
    lower_rel = (lower__ / mean(estimate__)) * 100,
    upper_rel = (upper__ / mean(estimate__)) * 100
  )

df_D <- ce_socio[[1]] %>%
  transmute(
    x_orig = stat_sosioekon_indeksi_skaalattu * p_socio$sd + p_socio$mu,
    estimate_rel = (estimate__ / mean(estimate__)) * 100,
    lower_rel = (lower__ / mean(estimate__)) * 100,
    upper_rel = (upper__ / mean(estimate__)) * 100
  )

all_lowers <- c(df_A$lower_rel, df_B$lower_rel, df_C$lower_rel, df_D$lower_rel)
all_uppers <- c(df_A$upper_rel, df_B$upper_rel, df_C$upper_rel, df_D$upper_rel)
y_min <- floor(min(all_lowers) / 10) * 10
y_max <- ceiling(max(all_uppers) / 10) * 10

# plot function
base_plot <- function(df, x_label, tag_letter, show_y_axis = TRUE) {
  p <- ggplot(df, aes(x = x_orig, y = estimate_rel)) +
    geom_ribbon(aes(ymin = lower_rel, ymax = upper_rel), alpha = 0.2, fill = "steelblue") +
    geom_line(linewidth = 1, color = "darkblue") +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = x_label, title = tag_letter) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank()
    )
  
  if (show_y_axis) {
    p <- p + labs(y = "% of mean abundance")
  } else {
    p <- p + theme(
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
  }
  return(p)
}

plot_A <- base_plot(df_A, "Temperature (7-day mean, °C)", 
                    "A", show_y_axis = TRUE)
plot_B <- base_plot(df_B, "Prolonged cold (6-month sum of days ≤ -10 °C)", 
                    "B", show_y_axis = FALSE)
plot_C <- base_plot(df_C, "Prolonged drought (90-day mean wildfire index)", 
                    "C", show_y_axis = TRUE)
plot_D <- base_plot(df_D, "Regional Socioeconomic Index (RSI)", 
                    "D", show_y_axis = FALSE)

fig_2 <- (plot_A + plot_B) / (plot_C + plot_D)
fig_2
# 
# ggsave(
#   filename = "fig2_rat_abundance_effects.png",
#   plot = fig_2,
#   device = "png",
#   width = 180,
#   height = 160,
#   units = "mm",
#   dpi = 300
# )
# 
# ggsave(
#   filename = "fig2_rat_abundance_effects.tiff",
#   plot = fig_2,
#   device = "tiff",
#   width = 180,
#   height = 160,
#   units = "mm",
#   dpi = 300,
#   compression = "lzw"
# )
#####

# plot greenspace effects (fig3)
#####

# Calculate effect for each season
effects_df <- as.data.frame(helsinkirats_model) %>% 
  transmute(
    Winter = b_vih_osuus_rak_log_skaalattu,
    Spring = b_vih_osuus_rak_log_skaalattu + `b_vih_osuus_rak_log_skaalattu:vuodenaikab_kevät`,
    Summer = b_vih_osuus_rak_log_skaalattu + `b_vih_osuus_rak_log_skaalattu:vuodenaikac_kesä`,
    Autumn = b_vih_osuus_rak_log_skaalattu + `b_vih_osuus_rak_log_skaalattu:vuodenaikad_syksy`
  ) %>% 
  pivot_longer(cols = everything(), names_to = "Season", values_to = "Effect") %>% 
  mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn")))

# Posterior mean and 95% credible interval
effects_summary <- effects_df %>% 
  group_by(Season) %>% 
  summarise(
    estimate = mean(Effect),
    lower = quantile(Effect, 0.025),
    upper = quantile(Effect, 0.975)
  )

# Plot
fig_3 <- ggplot(effects_summary, aes(x = Season, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") + # zeroline for comparison
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15, linewidth = 1.2) +
  geom_point(size = 4) +
  labs(
    x = "Season",
    y = "Effect of green area on rat abundance"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")
fig_3


#####

# Remove temporary objects
#####
rm(
  pkg_list,
  theformula,
  prop_zero,
  get_scaling_params,
  base_plot,
  p_temp,
  p_winter,
  p_mpi,
  p_socio,
  weather_vars,
  ce_weather,
  ce_socio,
  df_A,
  df_B,
  df_C,
  df_D,
  all_lowers,
  all_uppers,
  y_min,
  y_max,
  effects_df,
  effects_summary
)
#####
