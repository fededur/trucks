library(jsonlite)
library(ggplot2)
config <- fromJSON("config.json", simplifyVector = TRUE)
`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0L || is.na(value[1])) fallback else value
}

# Spatial Cascading State-Transition Monte Carlo Engine ----------------------
# Generic A -> B -> C -> D state engine. Domain meaning, durations, spatial
# extent, production value and protective attributes are supplied by JSON.
# The simulation uses jsonlite for configuration, ggplot2 for presentation,
# and otherwise relies on base R.

# ---- Required parameter mapping --------------------------------------------
MONTE_CARLO_RUNS <- config$simulation_controls$monte_carlo_iterations
POPULATION_SIZE <- config$model_parameters$population_size
SIM_DURATION <- config$model_parameters$simulation_days
LAMBDA <- config$model_parameters$spatial_decay_factor
COOLDOWN_PHASE_B <- config$model_parameters$incubation_period_days
COOLDOWN_PHASE_C <- config$model_parameters$active_period_days
BASE_YIELD <- config$production_metrics$average_bird_yield
ALLOC_RATE <- config$housing_attributes$barn_allocation_rate
PROT_EFFICIENCY <- config$housing_attributes$barn_protection_efficiency

# ---- Dynamic scenario labels -----------------------------------------------
SCENARIO_NAME <- config$scenario_metadata$name
PHASE_A_LABEL <- config$scenario_metadata$phase_a_label
PHASE_B_LABEL <- config$scenario_metadata$phase_b_label
PHASE_C_LABEL <- config$scenario_metadata$phase_c_label
PHASE_D_LABEL <- config$scenario_metadata$phase_d_label

# Optional labels have generic fallbacks, keeping older configuration files
# valid if presentation metadata is not supplied.
VALUE_AXIS_LABEL <- config$scenario_metadata$value_axis_label %||%
  paste("Total", PHASE_D_LABEL)
SURVIVAL_AXIS_LABEL <- config$scenario_metadata$survival_axis_label %||%
  paste("Average final", PHASE_A_LABEL, "rate (%)")
SHELTERED_LABEL <- config$scenario_metadata$sheltered_label %||% "Sheltered"
UNSHELTERED_LABEL <- config$scenario_metadata$unsheltered_label %||% "Unsheltered"

# ---- Optional configuration with safe defaults -----------------------------
RANDOM_SEED <- config$simulation_controls$random_seed %||% 2026L
INITIAL_PHASE_C <- config$simulation_controls$initial_phase_c_entities %||% 3L
MIN_LON <- config$spatial_bounds$minimum_longitude %||% 166
MAX_LON <- config$spatial_bounds$maximum_longitude %||% 179
MIN_LAT <- config$spatial_bounds$minimum_latitude %||% -47.5
MAX_LAT <- config$spatial_bounds$maximum_latitude %||% -34

# ---- Configuration validation ----------------------------------------------
assert_scalar <- function(value, name, lower = -Inf, upper = Inf,
                          integer = FALSE) {
  valid <- is.numeric(value) && length(value) == 1L && is.finite(value) &&
    value >= lower && value <= upper
  if (integer) valid <- valid && value == as.integer(value)
  if (!valid) {
    stop(sprintf(
      "Invalid config value '%s'. Expected %s in [%s, %s].",
      name, if (integer) "one integer" else "one number", lower, upper
    ), call. = FALSE)
  }
  invisible(TRUE)
}

assert_scalar(MONTE_CARLO_RUNS, "monte_carlo_iterations", 1, Inf, TRUE)
assert_scalar(POPULATION_SIZE, "population_size", 3, Inf, TRUE)
assert_scalar(SIM_DURATION, "simulation_days", 1, Inf, TRUE)
assert_scalar(LAMBDA, "spatial_decay_factor", 0, Inf)
assert_scalar(COOLDOWN_PHASE_B, "incubation_period_days", 1, Inf, TRUE)
assert_scalar(COOLDOWN_PHASE_C, "active_period_days", 1, Inf, TRUE)
assert_scalar(BASE_YIELD, "average_bird_yield", 0, Inf)
assert_scalar(ALLOC_RATE, "barn_allocation_rate", 0, 1)
assert_scalar(PROT_EFFICIENCY, "barn_protection_efficiency", 0, 1)
assert_scalar(INITIAL_PHASE_C, "initial_phase_c_entities", 1,
              POPULATION_SIZE, TRUE)
assert_scalar(MIN_LON, "minimum_longitude", -180, 180)
assert_scalar(MAX_LON, "maximum_longitude", -180, 180)
assert_scalar(MIN_LAT, "minimum_latitude", -90, 90)
assert_scalar(MAX_LAT, "maximum_latitude", -90, 90)
if (MIN_LON >= MAX_LON || MIN_LAT >= MAX_LAT) {
  stop("Spatial minimum bounds must be below maximum bounds.", call. = FALSE)
}

required_labels <- c(PHASE_A_LABEL, PHASE_B_LABEL, PHASE_C_LABEL,
                     PHASE_D_LABEL, SCENARIO_NAME)
if (any(!nzchar(required_labels))) {
  stop("Scenario name and all four phase labels must be non-empty.",
       call. = FALSE)
}

# ---- Pure geographic distance function -------------------------------------
# Returns great-circle distance in kilometres. Arguments may be scalars or
# equal-length vectors; no external state is read or changed.
haversine_distance <- function(lon1, lat1, lon2, lat2) {
  earth_radius_km <- 6371.0088
  radians <- pi / 180
  phi1 <- lat1 * radians
  phi2 <- lat2 * radians
  delta_phi <- (lat2 - lat1) * radians
  delta_lambda <- (lon2 - lon1) * radians
  a <- sin(delta_phi / 2)^2 +
    cos(phi1) * cos(phi2) * sin(delta_lambda / 2)^2
  2 * earth_radius_km * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

# ---- Monte Carlo state engine ----------------------------------------------
set.seed(as.integer(RANDOM_SEED))
production_loss_log <- numeric(MONTE_CARLO_RUNS)
sheltered_survival_log <- numeric(MONTE_CARLO_RUNS)
unsheltered_survival_log <- numeric(MONTE_CARLO_RUNS)

for (sim in seq_len(MONTE_CARLO_RUNS)) {
  # A brand-new spatial population is generated for every Monte Carlo run.
  # This represents uncertainty in the geographic arrangement of agents.
  sheltered_count <- as.integer(round(POPULATION_SIZE * ALLOC_RATE))
  sheltered_assignment <- rep(FALSE, POPULATION_SIZE)
  if (sheltered_count > 0L) {
    sheltered_assignment[sample.int(POPULATION_SIZE, sheltered_count)] <- TRUE
  }

  nodes <- data.frame(
    id = seq_len(POPULATION_SIZE),
    lon = runif(POPULATION_SIZE, MIN_LON, MAX_LON),
    lat = runif(POPULATION_SIZE, MIN_LAT, MAX_LAT),
    current_state = rep("A", POPULATION_SIZE),
    timer_phase_b = integer(POPULATION_SIZE),
    timer_phase_c = integer(POPULATION_SIZE),
    is_sheltered = sheltered_assignment,
    production_yield = pmax(
      0,
      rnorm(POPULATION_SIZE, mean = BASE_YIELD, sd = BASE_YIELD * 0.15)
    ),
    stringsAsFactors = FALSE
  )

  # Geography is fixed within a run, so compute all pairwise distances once.
  # Matrix rows are targets and columns are potential active sources. Reusing
  # this matrix avoids expensive geographic recalculation on every day.
  spatial_distance_km <- matrix(
    haversine_distance(
      lon1 = rep(nodes$lon, times = POPULATION_SIZE),
      lat1 = rep(nodes$lat, times = POPULATION_SIZE),
      lon2 = rep(nodes$lon, each = POPULATION_SIZE),
      lat2 = rep(nodes$lat, each = POPULATION_SIZE)
    ),
    nrow = POPULATION_SIZE,
    ncol = POPULATION_SIZE
  )

  # Exactly three initial active entities by default, configurable if needed.
  seed_ids <- sample.int(POPULATION_SIZE, INITIAL_PHASE_C)
  nodes$current_state[seed_ids] <- "C"

  for (day in seq_len(SIM_DURATION)) {
    active_ids <- which(nodes$current_state == "C")
    target_ids <- which(nodes$current_state == "A")
    phase_b_at_start_of_timing <- which(nodes$current_state == "B")
    phase_c_at_start_of_timing <- active_ids

    if (length(active_ids) > 0L && length(target_ids) > 0L) {
      newly_exposed <- logical(length(target_ids))

      for (target_position in seq_along(target_ids)) {
        target_id <- target_ids[target_position]
        distances <- spatial_distance_km[target_id, active_ids]

        raw_probability <- exp(-LAMBDA * distances)
        protection_multiplier <- if (nodes$is_sheltered[target_id]) {
          1 - PROT_EFFICIENCY
        } else {
          1
        }
        source_probabilities <- pmin(
          1, pmax(0, raw_probability * protection_multiplier)
        )

        # Independent exposure routes are combined into the probability that
        # at least one active node reaches this target during the current day.
        final_probability <- 1 - prod(1 - source_probabilities)
        newly_exposed[target_position] <-
          runif(1) <= final_probability
      }

      entering_phase_b <- target_ids[newly_exposed]
      if (length(entering_phase_b) > 0L) {
        nodes$current_state[entering_phase_b] <- "B"
        nodes$timer_phase_b[entering_phase_b] <- 0L
      }
    }

    # Only nodes already in a timed phase at the start of this day advance their
    # clocks. A newly exposed node therefore enters B at timer zero.
    if (length(phase_b_at_start_of_timing) > 0L) {
      nodes$timer_phase_b[phase_b_at_start_of_timing] <-
        nodes$timer_phase_b[phase_b_at_start_of_timing] + 1L
      entering_phase_c <- phase_b_at_start_of_timing[
        nodes$timer_phase_b[phase_b_at_start_of_timing] >= COOLDOWN_PHASE_B
      ]
      if (length(entering_phase_c) > 0L) {
        nodes$current_state[entering_phase_c] <- "C"
        nodes$timer_phase_c[entering_phase_c] <- 0L
      }
    }

    if (length(phase_c_at_start_of_timing) > 0L) {
      nodes$timer_phase_c[phase_c_at_start_of_timing] <-
        nodes$timer_phase_c[phase_c_at_start_of_timing] + 1L
      entering_phase_d <- phase_c_at_start_of_timing[
        nodes$timer_phase_c[phase_c_at_start_of_timing] >= COOLDOWN_PHASE_C
      ]
      if (length(entering_phase_d) > 0L) {
        nodes$current_state[entering_phase_d] <- "D"
      }
    }
  }

  production_loss_log[sim] <- sum(
    nodes$production_yield[nodes$current_state == "D"]
  )

  sheltered_total <- sum(nodes$is_sheltered)
  unsheltered_total <- sum(!nodes$is_sheltered)
  sheltered_survival_log[sim] <- if (sheltered_total > 0L) {
    100 * sum(nodes$current_state == "A" & nodes$is_sheltered) /
      sheltered_total
  } else {
    NA_real_
  }
  unsheltered_survival_log[sim] <- if (unsheltered_total > 0L) {
    100 * sum(nodes$current_state == "A" & !nodes$is_sheltered) /
      unsheltered_total
  } else {
    NA_real_
  }
}

# ---- Results and ggplot2 visualisation --------------------------------------
average_survival <- c(
  mean(sheltered_survival_log, na.rm = TRUE),
  mean(unsheltered_survival_log, na.rm = TRUE)
)
names(average_survival) <- c(SHELTERED_LABEL, UNSHELTERED_LABEL)

simulation_results <- list(
  scenario = SCENARIO_NAME,
  configuration = config,
  production_loss_log = production_loss_log,
  sheltered_survival_log = sheltered_survival_log,
  unsheltered_survival_log = unsheltered_survival_log,
  average_survival = average_survival,
  phase_labels = c(A = PHASE_A_LABEL, B = PHASE_B_LABEL,
                   C = PHASE_C_LABEL, D = PHASE_D_LABEL)
)

assets_dir <- "assets"
if (!dir.exists(assets_dir)) dir.create(assets_dir, recursive = TRUE)
saveRDS(simulation_results,
        file.path(assets_dir, "spatial_cascade_results.rds"))

loss_data <- data.frame(production_loss = production_loss_log)
survival_data <- data.frame(
  protection_group = factor(
    names(average_survival), levels = names(average_survival)
  ),
  survival_rate = as.numeric(average_survival)
)

professional_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", colour = "#173B57"),
    plot.subtitle = element_text(colour = "#526575"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

loss_plot <- ggplot(loss_data, aes(x = production_loss)) +
  geom_histogram(
    bins = max(12L, as.integer(round(sqrt(MONTE_CARLO_RUNS)))),
    fill = "#2A788E", colour = "white", linewidth = 0.25
  ) +
  geom_vline(
    xintercept = mean(production_loss_log),
    colour = "#E76F51", linewidth = 1, linetype = "dashed"
  ) +
  labs(
    title = paste(SCENARIO_NAME, "loss distribution"),
    subtitle = sprintf(
      "%s Monte Carlo runs; dashed line shows the average outcome",
      format(MONTE_CARLO_RUNS, big.mark = ",")
    ),
    x = VALUE_AXIS_LABEL,
    y = "Number of simulation runs"
  ) +
  scale_x_continuous(labels = function(x) format(x, big.mark = ",",
                                                 scientific = FALSE)) +
  professional_theme

survival_plot <- ggplot(
  survival_data,
  aes(x = protection_group, y = survival_rate, fill = protection_group)
) +
  geom_col(width = 0.62) +
  geom_text(
    aes(label = sprintf("%.1f%%", survival_rate)),
    vjust = -0.5, fontface = "bold", colour = "#173B57"
  ) +
  scale_fill_manual(values = c("#3B9AB2", "#F28E2B")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = paste("Average final", PHASE_A_LABEL, "rate"),
    subtitle = paste("Comparison by protective attribute across",
                     format(MONTE_CARLO_RUNS, big.mark = ","), "runs"),
    x = NULL,
    y = SURVIVAL_AXIS_LABEL
  ) +
  professional_theme

ggsave(
  file.path(assets_dir, "spatial_cascade_loss_distribution.svg"),
  loss_plot, device = grDevices::svg, width = 9, height = 6, units = "in",
  bg = "white"
)
ggsave(
  file.path(assets_dir, "spatial_cascade_survival_comparison.svg"),
  survival_plot, device = grDevices::svg, width = 8, height = 6, units = "in",
  bg = "white"
)

summary_table <- data.frame(
  scenario = SCENARIO_NAME,
  monte_carlo_runs = MONTE_CARLO_RUNS,
  mean_production_loss = mean(production_loss_log),
  median_production_loss = median(production_loss_log),
  loss_p05 = unname(quantile(production_loss_log, 0.05)),
  loss_p95 = unname(quantile(production_loss_log, 0.95)),
  sheltered_survival_pct = unname(average_survival[1]),
  unsheltered_survival_pct = unname(average_survival[2]),
  protection_benefit_percentage_points =
    unname(average_survival[1] - average_survival[2]),
  stringsAsFactors = FALSE
)
write.csv(summary_table,
          file.path(assets_dir, "spatial_cascade_summary.csv"),
          row.names = FALSE)

if (interactive()) {
  print(loss_plot)
  print(survival_plot)
}

invisible(simulation_results)
