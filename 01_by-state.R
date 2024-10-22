library(tidyverse)
library(brms)
library(lme4)
library(cmdstanr) # remotes::install_github("stan-dev/cmdstanr")
cmdstanr::install_cmdstan()

# example data: y = repeal the ACA 
dat <- read_csv("data/cces_by-cd.csv.gz")


# try a sample so each model runs in a minute ---
set.seed(06510)
samp_dat <- dat |> slice_sample(prop = 0.4, by = st, replace = FALSE)

## Variable to use
grpvar = ~st
grp_ff = update(grpvar, y ~ (1 | . ))

# Direct estimates
dir_est <- samp_dat |> 
  summarize(
    pct_raw = mean(y),
    n  = n(),
    .by = all_of(all.vars(grpvar))
  )


# states (poststrat)
st_dat <- distinct(samp_dat, .data[[all.vars(grpvar)]])


# lmer (empirical bayes?) 
## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2009q4/002984.html
fit_lmer <- glmer(grp_ff, samp_dat, family = binomial)
pred_lmer <- st_dat |> 
  mutate(yhat = predict(fit_lmer, st_dat, type = "response"), 
         model = "lmer")

#' Function to run brms random effects
fit_pool <- function(
    prior = prior_string("normal(0, 1)", class = "sd"), 
    form = grp_ff, 
    data = samp_dat,
    type = "no") {
  brms::brm(
    form, 
    data, 
    family = bernoulli, 
    iter = 1000,
    prior = prior, 
    cores = 4,
    sample_prior = type,
    backend = "cmdstanr")
}

# different priors
ssfit0 <- fit_pool(prior = NULL)
ssfit1 <- fit_pool(prior_string("cauchy(0, 0.1)", class = "sd"))
ssfit2 <- fit_pool(prior_string("cauchy(0, 10)", class = "sd"))
ssfit3 <- fit_pool(c(prior_string("normal(1, 0.1)", class = "sd")))
ssfit4 <- fit_pool(c(prior_string("normal(1, 1)", class = "sd")))
ssfit5 <- fit_pool(c(prior_string("normal(5, 1)", class = "sd")))
ssfitf1 <- fit_pool(form = update(grpvar, y ~ .), prior = prior_string("normal(0, 1)", class = "b"))

# create estimates
ests <- list(
  `0: sigma ~ half studentt_3` = ssfit0,
  `1: sigma ~ half cauchy(0, 0.1)` = ssfit1,
  `2: sigma ~ half cauchy(0, 10)` = ssfit2,
  `3: sigma ~ half normal(1, 0.1)` = ssfit3,
  `4: sigma ~ half normal(1, 1)` = ssfit4,
  `5: sigma ~ half normal(5, 1)` = ssfit5,
  `7: FE` = ssfitf1
) |> 
  map(
    .f = \(x) {
      spred <- x |> 
        posterior_epred(newdata = st_dat)
      colnames(spred) <- st_dat[[1]]
      
      sdraws <- spred |> 
        as_tibble() |> 
        pivot_longer(cols = everything(), names_to = "st", values_to = "yhat") |> 
        summarize(
          sd = sd(yhat), 
          yhat = mean(yhat), 
          n_iter = n(), 
          .by = all_of(all.vars(grpvar)))
    }
  ) |> 
  list_rbind(names_to = "model")


# visualize -----
ests |> 
  bind_rows(pred_lmer) |> 
  left_join(dir_est, by = c("st")) |> 
  ggplot(aes(x = pct_raw, y = yhat)) +
  # geom_point(size = 0.5) + 
  geom_text(aes_(label = grpvar), size = 1) + 
  lemon::facet_rep_wrap(
    ~ model, labeller = label_wrap_gen(),
    nrow = 2) +
  geom_abline(linewidth = 0.1) +
  theme_classic() +
  coord_equal() +
  labs(x = "Direct estimate", y = "Modeled Estimate") +
  labs(caption = glue::glue("RE models are {as.character(grp_ff)}. FE models are y ~ state. n = {nrow(samp_dat)}"))
