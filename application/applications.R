library(tidyverse)

##### Setup #####

scale_outcome <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

run_analysis <- function(
  data,
  X,
  A,
  Y,
  dataset_name,
  thresholds = 10^seq(-5, log10(0.05), 0.2),
  smoothness = c(0.001)
) {
  cat("Analyzing: ", dataset_name, "\n")

  set.seed(10016)
  bounds <- ate_bounds(
    data[, c(X, A, Y)],
    X,
    A,
    Y,
    smoothness = smoothness,
    thresholds = thresholds,
    learners_trt = c("SL.glm"),
    learners_outcome = c("SL.glm")
  )

  pi_hat <- bounds$nuisance$pi_hat

  calc_non_overlap_size <- function(threshold) {
    below <- pi_hat <= threshold
    above <- (1 - pi_hat) <= threshold
    mean(below | above)
  }
  cs <- c(seq(0, 0.5, length.out = 100), thresholds)
  non_overlap_size <- map_dbl(cs, calc_non_overlap_size)
  min_pscore <- summary(bounds)$propensity_score_range[1]
  onestep_width <- summary(bounds)$onestep$upper - summary(bounds)$onestep$lower

  min_observed_pscore <- min(c(
    bounds$nuisance$pi_hat[data[[A]] == 1],
    (1 - bounds$nuisance$pi_hat[data[[A]] == 0])
  ))

  plot_pscore_distribution <- ggplot(tibble(pi = pi_hat), aes(x = pi)) +
    geom_histogram(color = "white", bins = 50, alpha = 0.7) +
    labs(
      title = paste0(dataset_name, ": Estimated Propensity Score Distribution"),
      x = expression(hat(pi)(X)),
      y = "Count"
    ) +
    theme_minimal()

  plot_log_pscore_distribution <- ggplot(
    tibble(pi = pi_hat),
    aes(x = 1 / pi)
  ) +
    geom_histogram(color = "white", bins = 100, alpha = 0.7) +
    labs(
      title = paste0(dataset_name, ": Estimated Propensity Score Distribution"),
      x = expression(hat(pi)(X)),
      y = "Count"
    ) +
    scale_x_log10() +
    theme_minimal()

  plot_nonoverlap_size <- ggplot(
    tibble(c = cs, size = non_overlap_size),
    aes(x = 1 / c, y = size)
  ) +
    geom_line(linewidth = 0.8) +
    labs(
      title = paste0(dataset_name),
      x = expression(log[10](1 / c)),
      y = expression(P(hat(pi) <= c) + P(hat(pi) >= 1 - c))
    ) +
    scale_x_log10(
      limits = c(1, 1 / min_observed_pscore)
    ) +
    scale_y_continuous(
      labels = scales::percent_format()
    ) +
    theme_minimal()

  uniform_width <- min(
    summary(bounds)$shortest_intervals$upper_uniform -
      summary(bounds)$shortest_intervals$lower_uniform
  )

  list(
    dataset_name = dataset_name,
    n = nrow(data),
    n_covariates = length(X),
    prevalence = mean(data[[A]]),
    bounds = bounds,
    cs = cs,
    pi_hat = pi_hat,
    thresholds = thresholds,
    non_overlap_size = non_overlap_size,
    min_pscore = min_pscore,
    min_observed_pscore = min_observed_pscore,
    onestep_width = onestep_width,
    uniform_width = uniform_width,
    plot_pscore_distribution = plot_pscore_distribution,
    plot_log_pscore_distribution = plot_log_pscore_distribution,
    plot_nonoverlap_size = plot_nonoverlap_size
  )
}

##### Right-heart catheterization #####
# Package: ATbounds
# Treatment: right-heart catheterization
# Outcome: mortality (binary)
load_rhc <- function() {
  data("RHC", package = "ATbounds")

  X <- setdiff(colnames(RHC), c("survival", "RHC"))
  A <- "RHC"
  Y <- "survival"

  list(
    data = as_tibble(RHC),
    X = X,
    A = A,
    Y = Y,
    name = "Right heart catheterization"
  )
}

##### Politicians #####
# Package: causaldata
# Treatment: legislator race
# Outcome: response
load_politicians <- function() {
  data("black_politicians", package = "causaldata")
  df <- as_tibble(black_politicians)
  Y <- "responded"
  A <- "leg_black"
  X <- setdiff(colnames(df), c(A, Y))

  list(
    data = df,
    X = X,
    A = A,
    Y = Y,
    name = "Black Politicians"
  )
}

##### LaLonde / NSW (Deheijia-Wahba Sample) #####
# Package: causaldata
# Treatment: job program enrollment
# Outcome: wages
load_lalonde <- function() {
  data("nsw_mixtape", package = "causaldata")

  df <- as_tibble(nsw_mixtape)
  df$Y <- scale_outcome(df$re78)
  X <- c(
    "age",
    "educ",
    "black",
    "hisp",
    "marr",
    "nodegree",
    "re74",
    "re75"
  )

  list(
    data = df,
    X = X,
    A = "treat",
    Y = "Y",
    name = "LaLonde (Dehejia-Wahba)"
  )
}

##### NHEFS (National Health and Nutrition Examination Survey) #####
# Package: causaldata
# Treatment: smoking cessation
# Outcome: weight change
load_nhefs <- function() {
  data("nhefs", package = "causaldata")
  nhefs <- as_tibble(nhefs[!is.na(nhefs$wt82_71), ])
  nhefs$Y <- scale_outcome(nhefs$wt82_71)

  X <- c(
    "sex",
    "age",
    "race",
    "education",
    "ht",
    "smokeintensity",
    "smokeyrs",
    "exercise",
    "active",
    "wt71"
  )

  list(data = nhefs, X = X, A = "qsmk", Y = "Y", name = "NHEFS")
}

##### STAR #####
# Package: AER
# Treatment: small class vs. regular class
# Outcome: math scores
load_star <- function() {
  data("STAR", package = "AER")

  star <- as_tibble(STAR) |>
    filter(stark %in% c("small", "regular")) |>
    filter(!is.na(mathk)) |>
    mutate(treat = as.integer(stark == "small"), Y = scale_outcome(mathk))

  star <- na.omit(star[, c(
    "gender",
    "ethnicity",
    "birth",
    "lunchk",
    "treat",
    "Y"
  )])

  star <- as_tibble(model.matrix(
    ~ -1 + Y + gender + ethnicity + birth + lunchk + treat,
    data = star
  ))

  X <- c(
    "gendermale",
    "genderfemale",
    "ethnicityafam",
    "ethnicityasian",
    "ethnicityhispanic",
    "ethnicityamindian",
    "ethnicityother",
    "birth",
    "lunchkfree"
  )

  list(
    data = star[, c(X, "treat", "Y")],
    X = X,
    A = "treat",
    Y = "Y",
    name = "STAR (Kindergarten Math)"
  )
}

##### Crown Court Sentencing #####
load_ccdrug <- function() {
  data("ccdrug", package = "causaldata")
  df <- as_tibble(ccdrug)

  Y <- "custody"
  A <- "male"
  X <- setdiff(colnames(df), c(A, Y))

  list(
    data = df,
    X = X,
    A = A,
    Y = Y,
    name = "Crown Court Sentencing"
  )
}

##### Results #####

results <- tibble(
  loader = list(
    load_ccdrug,
    load_politicians,
    load_rhc,
    load_lalonde,
    load_nhefs,
    load_star
  ),
) |>
  mutate(
    data = map(loader, \(f) f())
  )

results <- results |>
  mutate(
    results = map(data, \(data) {
      run_analysis(data$data, data$X, data$A, data$Y, data$name)
    })
  )

results |>
  mutate(
    name = map_chr(data, `[[`, "name"),
    n = map_dbl(results, `[[`, "n"),
    p = map_dbl(results, `[[`, "n_covariates"),
    prev = map_dbl(results, `[[`, "prevalence"),
    uniform_width = pmin(2, map_dbl(results, `[[`, "uniform_width")),
    onestep_width = pmin(2, map_dbl(results, `[[`, "onestep_width"))
  ) |>
  dplyr::select(-loader, -data, -results) |>
  arrange(n) |>
  mutate_at(vars(ends_with("_width")), scales::number_format(0.01)) |>
  mutate_at(vars(prev), scales::percent_format(accuracy = 1))

plot(results$results[[1]]$bounds)

results$results[[1]]$plot_nonoverlap_size

summary(results$results[[1]]$bounds)


combined_data <- map_dfr(
  results$results,
  \(results) {
    tibble(
      Dataset = results$dataset_name,
      c = results$cs,
      size = results$non_overlap_size
    )
  }
)

ggplot(combined_data, aes(x = c, y = size, color = Dataset)) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Non-overlap Subpopulation Size Across Application Datasets",
    x = expression(c),
    y = expression(P(hat(pi) <= c) + P(hat(pi) >= 1 - c))
  ) +
  theme_minimal()
ggsave("application/propensity_score_distributions.pdf", width = 6, height = 6)

pdf("application/non_overlap_size.pdf", width = 8, height = 6)
do.call(
  gridExtra::grid.arrange,
  map(results$results, `[[`, "plot_nonoverlap_size")
)
dev.off()


pdf(
  "application/rhc_non_overlap_size_pscore_distribution.pdf",
  width = 8,
  height = 4
)
gridExtra::grid.arrange(
  results$results[[3]]$plot_pscore_distribution +
    ggtitle("Estimated Propensity Score Distribution"),
  results$results[[3]]$plot_nonoverlap_size +
    ggtitle("Non-overlap subpopulation size"),
  ncol = 2
)
dev.off()
