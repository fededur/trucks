library(jsonlite)
library(ggplot2)
config <- fromJSON("config.json", simplifyVector = TRUE)

`%||%` <- function(x, fallback) if (is.null(x) || !length(x) || is.na(x[1])) fallback else x

# Fixed-Population Spatial Cascade Monte Carlo Engine ------------------------
# Farms and attributes remain fixed. Uncertainty comes from event origins and
# transmission rolls, not from inventing a new population for every run.

MONTE_CARLO_RUNS <- config$simulation_controls$monte_carlo_iterations
POPULATION_SIZE <- config$model_parameters$population_size
SIM_DURATION <- config$model_parameters$simulation_days
LAMBDA <- config$model_parameters$spatial_decay_factor
CROSS_BARRIER_MULTIPLIER <-
  config$model_parameters$cross_barrier_transmission_multiplier %||% 0
COOLDOWN_PHASE_B <- config$model_parameters$incubation_period_days
COOLDOWN_PHASE_C <- config$model_parameters$active_period_days
BASE_YIELD <- config$production_metrics$average_bird_yield
ALLOC_RATE <- config$housing_attributes$barn_allocation_rate
PROT_EFFICIENCY <- config$housing_attributes$barn_protection_efficiency

meta <- config$scenario_metadata
SCENARIO_NAME <- meta$name
PHASE_A_LABEL <- meta$phase_a_label
PHASE_B_LABEL <- meta$phase_b_label
PHASE_C_LABEL <- meta$phase_c_label
PHASE_D_LABEL <- meta$phase_d_label
VALUE_AXIS_LABEL <- meta$value_axis_label %||% paste("Total", PHASE_D_LABEL)
SURVIVAL_AXIS_LABEL <- meta$survival_axis_label %||% paste("Final", PHASE_A_LABEL, "rate (%)")
SHELTERED_LABEL <- meta$sheltered_label %||% "Sheltered"
UNSHELTERED_LABEL <- meta$unsheltered_label %||% "Unsheltered"
PRODUCTION_TYPE_LABEL <- meta$production_type_label %||% "Production type"

RANDOM_SEED <- config$simulation_controls$random_seed %||% 2026L
INITIAL_PHASE_C <- config$simulation_controls$initial_phase_c_entities %||% 3L
START_MODE <- config$simulation_controls$starting_location_mode %||% "random"
ITERATIONS_PER_START <- config$simulation_controls$iterations_per_starting_farm %||% 25L
INPUT_FILE <- config$input_data$file %||% "data/farms.csv"
COLUMN_MAP <- unlist(config$input_data$column_mappings, use.names = TRUE)

assert_number <- function(x, name, lower = -Inf, integer = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lower
  if (integer) ok <- ok && x == as.integer(x)
  if (!ok) stop("Invalid configuration value: ", name, call. = FALSE)
}
assert_number(MONTE_CARLO_RUNS, "monte_carlo_iterations", 1, TRUE)
assert_number(SIM_DURATION, "simulation_days", 1, TRUE)
assert_number(LAMBDA, "spatial_decay_factor", 0)
assert_number(CROSS_BARRIER_MULTIPLIER,
              "cross_barrier_transmission_multiplier", 0)
if (CROSS_BARRIER_MULTIPLIER > 1) {
  stop("cross_barrier_transmission_multiplier must be between 0 and 1.")
}
assert_number(COOLDOWN_PHASE_B, "incubation_period_days", 1, TRUE)
assert_number(COOLDOWN_PHASE_C, "active_period_days", 1, TRUE)
assert_number(PROT_EFFICIENCY, "barn_protection_efficiency", 0)
if (PROT_EFFICIENCY > 1) stop("Protection efficiency must be at most 1.")
if (!START_MODE %in% c("random", "each_farm")) stop("starting_location_mode must be random or each_farm.")

haversine_distance <- function(lon1, lat1, lon2, lat2) {
  r <- 6371.0088; rad <- pi / 180
  phi1 <- lat1 * rad; phi2 <- lat2 * rad
  dphi <- (lat2 - lat1) * rad; dlambda <- (lon2 - lon1) * rad
  a <- sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlambda / 2)^2
  2 * r * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

# Load and standardise the fixed farm table.
required_fields <- c("id", "lon", "lat", "barrier_group", "region",
                     "production_type", "is_sheltered", "production_yield")
if (!file.exists(INPUT_FILE)) stop("Farm input file not found: ", INPUT_FILE)
if (!all(required_fields %in% names(COLUMN_MAP))) stop("Incomplete input column_mappings.")
raw_farms <- read.csv(INPUT_FILE, stringsAsFactors = FALSE, check.names = FALSE)
missing_columns <- setdiff(unname(COLUMN_MAP[required_fields]), names(raw_farms))
if (length(missing_columns)) stop("Farm input is missing: ", paste(missing_columns, collapse = ", "))
farms <- raw_farms[, unname(COLUMN_MAP[required_fields]), drop = FALSE]
names(farms) <- required_fields
farms$id <- as.character(farms$id)
farms$barrier_group <- as.character(farms$barrier_group)
farms$region <- as.character(farms$region)
farms$production_type <- as.character(farms$production_type)
farms$lon <- as.numeric(farms$lon); farms$lat <- as.numeric(farms$lat)
farms$production_yield <- as.numeric(farms$production_yield)
farms$is_sheltered <- as.logical(farms$is_sheltered)
if (anyDuplicated(farms$id)) stop("Farm IDs must be unique.")
if (anyNA(farms) || any(farms$production_yield < 0)) stop("Farm data contain invalid or missing values.")
POPULATION_SIZE <- nrow(farms)
if (POPULATION_SIZE < 3L || INITIAL_PHASE_C > POPULATION_SIZE) stop("Invalid fixed farm population size.")

# Fixed distances are computed once. Barrier pairs receive the configurable
# multiplier: 0 blocks crossings, values between 0 and 1 penalise them, and 1
# removes the barrier effect.
distance_km <- matrix(haversine_distance(
  rep(farms$lon, times = POPULATION_SIZE), rep(farms$lat, times = POPULATION_SIZE),
  rep(farms$lon, each = POPULATION_SIZE), rep(farms$lat, each = POPULATION_SIZE)
), nrow = POPULATION_SIZE)
same_group <- outer(farms$barrier_group, farms$barrier_group, `==`)
barrier_multiplier <- ifelse(same_group, 1, CROSS_BARRIER_MULTIPLIER)

set.seed(as.integer(RANDOM_SEED))
if (START_MODE == "random") {
  start_sets <- replicate(MONTE_CARLO_RUNS, sample.int(POPULATION_SIZE, INITIAL_PHASE_C), simplify = FALSE)
} else {
  assert_number(ITERATIONS_PER_START, "iterations_per_starting_farm", 1, TRUE)
  start_sets <- rep(as.list(seq_len(POPULATION_SIZE)), each = ITERATIONS_PER_START)
  MONTE_CARLO_RUNS <- length(start_sets)
}

types <- sort(unique(farms$production_type))
groups <- c(SHELTERED_LABEL, UNSHELTERED_LABEL)
production_loss_log <- numeric(MONTE_CARLO_RUNS)
type_loss_log <- matrix(0, MONTE_CARLO_RUNS, length(types), dimnames = list(NULL, types))
type_lost_count_log <- type_loss_log
type_survival_log <- array(NA_real_, c(MONTE_CARLO_RUNS, length(types), 2L),
                           dimnames = list(NULL, types, groups))
group_survival_log <- matrix(NA_real_, MONTE_CARLO_RUNS, 2L, dimnames = list(NULL, groups))
farm_start_count <- farm_affected_count <- farm_lost_count <- integer(POPULATION_SIZE)

for (sim in seq_len(MONTE_CARLO_RUNS)) {
  state <- rep("A", POPULATION_SIZE); timer_b <- timer_c <- integer(POPULATION_SIZE)
  seeds <- start_sets[[sim]]; state[seeds] <- "C"
  farm_start_count[seeds] <- farm_start_count[seeds] + 1L

  for (day in seq_len(SIM_DURATION)) {
    active <- which(state == "C"); targets <- which(state == "A")
    existing_b <- which(state == "B"); existing_c <- active
    if (length(active) && length(targets)) {
      new_b <- integer()
      for (target in targets) {
        sources <- active
        raw <- exp(-LAMBDA * distance_km[target, sources])
        protection <- if (farms$is_sheltered[target]) 1 - PROT_EFFICIENCY else 1
        adjusted <- raw * protection * barrier_multiplier[target, sources]
        final_probability <- 1 - prod(1 - pmin(1, pmax(0, adjusted)))
        if (runif(1) <= final_probability) new_b <- c(new_b, target)
      }
      if (length(new_b)) { state[new_b] <- "B"; timer_b[new_b] <- 0L }
    }
    if (length(existing_b)) {
      timer_b[existing_b] <- timer_b[existing_b] + 1L
      new_c <- existing_b[timer_b[existing_b] >= COOLDOWN_PHASE_B]
      if (length(new_c)) { state[new_c] <- "C"; timer_c[new_c] <- 0L }
    }
    if (length(existing_c)) {
      timer_c[existing_c] <- timer_c[existing_c] + 1L
      new_d <- existing_c[timer_c[existing_c] >= COOLDOWN_PHASE_C]
      if (length(new_d)) state[new_d] <- "D"
    }
  }

  affected <- state != "A"; lost <- state == "D"
  farm_affected_count <- farm_affected_count + affected
  farm_lost_count <- farm_lost_count + lost
  production_loss_log[sim] <- sum(farms$production_yield[lost])
  for (type in types) {
    type_rows <- farms$production_type == type
    type_loss_log[sim, type] <- sum(farms$production_yield[type_rows & lost])
    type_lost_count_log[sim, type] <- sum(type_rows & lost)
    for (group in groups) {
      sheltered <- identical(group, SHELTERED_LABEL)
      rows <- type_rows & farms$is_sheltered == sheltered
      if (sum(rows)) type_survival_log[sim, type, group] <- 100 * sum(rows & state == "A") / sum(rows)
    }
  }
  for (group in groups) {
    rows <- farms$is_sheltered == identical(group, SHELTERED_LABEL)
    group_survival_log[sim, group] <- 100 * sum(rows & state == "A") / sum(rows)
  }
}

average_type_survival <- apply(type_survival_log, c(2, 3), mean, na.rm = TRUE)
average_survival <- colMeans(group_survival_log, na.rm = TRUE)
assets_dir <- "assets"; if (!dir.exists(assets_dir)) dir.create(assets_dir, recursive = TRUE)

farm_risk <- transform(farms,
  times_selected_as_start = farm_start_count,
  probability_affected_pct = 100 * farm_affected_count / MONTE_CARLO_RUNS,
  probability_lost_pct = 100 * farm_lost_count / MONTE_CARLO_RUNS,
  expected_loss_contribution = production_yield * farm_lost_count / MONTE_CARLO_RUNS)
write.csv(farm_risk, file.path(assets_dir, "spatial_cascade_farm_risk.csv"), row.names = FALSE)

type_summary <- do.call(rbind, lapply(types, function(type) {
  losses <- type_loss_log[, type]; lost_counts <- type_lost_count_log[, type]
  population <- sum(farms$production_type == type)
  data.frame(production_type = type, farms_in_population = population,
    farms_in_population_pct = 100 * population / POPULATION_SIZE,
    average_farms_lost = mean(lost_counts), average_farms_lost_pct = 100 * mean(lost_counts) / population,
    mean_total_production_loss = mean(losses), sd_total_production_loss = sd(losses),
    total_loss_p25 = unname(quantile(losses, .25)), median_total_production_loss = median(losses),
    total_loss_p75 = unname(quantile(losses, .75)), total_loss_p05 = unname(quantile(losses, .05)),
    total_loss_p95 = unname(quantile(losses, .95)),
    sheltered_survival_pct = average_type_survival[type, SHELTERED_LABEL],
    unsheltered_survival_pct = average_type_survival[type, UNSHELTERED_LABEL])
}))
write.csv(type_summary, file.path(assets_dir, "spatial_cascade_production_type_summary.csv"), row.names = FALSE)

summary_table <- data.frame(scenario = SCENARIO_NAME, monte_carlo_runs = MONTE_CARLO_RUNS,
  population_size = POPULATION_SIZE, population_size_source = paste("fixed table", INPUT_FILE),
  starting_location_mode = START_MODE, mean_production_loss = mean(production_loss_log),
  median_production_loss = median(production_loss_log), loss_p05 = unname(quantile(production_loss_log, .05)),
  loss_p95 = unname(quantile(production_loss_log, .95)),
  sheltered_survival_pct = average_survival[SHELTERED_LABEL],
  unsheltered_survival_pct = average_survival[UNSHELTERED_LABEL],
  protection_benefit_percentage_points = average_survival[SHELTERED_LABEL] - average_survival[UNSHELTERED_LABEL])
write.csv(summary_table, file.path(assets_dir, "spatial_cascade_summary.csv"), row.names = FALSE)

results <- list(scenario = SCENARIO_NAME, configuration = config, farms = farms,
  population_size = POPULATION_SIZE, population_size_source = paste("fixed table", INPUT_FILE),
  starting_location_mode = START_MODE, monte_carlo_runs = MONTE_CARLO_RUNS,
  production_loss_log = production_loss_log, production_loss_by_type_log = type_loss_log,
  lost_count_by_type_log = type_lost_count_log, survival_by_type_and_shelter_log = type_survival_log,
  average_survival = average_survival, farm_risk_summary = farm_risk)
saveRDS(results, file.path(assets_dir, "spatial_cascade_results.rds"))

report_theme <- theme_minimal(base_size = 12) + theme(
  plot.title = element_text(face = "bold", colour = "#173B57"),
  plot.subtitle = element_text(colour = "#526575"), panel.grid.minor = element_blank(), legend.position = "none")
loss_plot <- ggplot(data.frame(loss = production_loss_log), aes(loss)) +
  geom_histogram(bins = max(12, round(sqrt(MONTE_CARLO_RUNS))), fill = "#2A788E", colour = "white") +
  geom_vline(xintercept = mean(production_loss_log), colour = "#E76F51", linewidth = 1, linetype = "dashed") +
  labs(title = paste(SCENARIO_NAME, "loss distribution"), subtitle = paste(MONTE_CARLO_RUNS, "runs over one fixed farm population"), x = VALUE_AXIS_LABEL, y = "Simulation runs") + report_theme
type_loss_data <- as.data.frame(as.table(type_loss_log)); names(type_loss_data) <- c("run", "production_type", "loss")
type_loss_plot <- ggplot(type_loss_data, aes(production_type, loss, fill = production_type)) +
  geom_boxplot(width = .62, outlier.alpha = .2) + scale_fill_manual(values = setNames(hcl.colors(length(types), "Dark 3"), types)) +
  labs(title = "Total loss by production type", subtitle = "Production type reports outcomes; it does not alter spread", x = PRODUCTION_TYPE_LABEL, y = VALUE_AXIS_LABEL) + report_theme
survival_data <- as.data.frame(as.table(average_type_survival)); names(survival_data) <- c("production_type", "protection_group", "survival")
survival_plot <- ggplot(survival_data, aes(production_type, survival, fill = protection_group)) +
  geom_col(position = position_dodge(.72), width = .66) + geom_text(aes(label = sprintf("%.1f%%", survival)), position = position_dodge(.72), vjust = -.4) +
  scale_fill_manual(values = c("#3B9AB2", "#F28E2B"), name = "Protection") + scale_y_continuous(limits = c(0, 100), expand = expansion(c(0, .04))) +
  labs(title = paste("Average final", PHASE_A_LABEL, "rate"), subtitle = "Fixed farms, split by production type and shelter", x = PRODUCTION_TYPE_LABEL, y = SURVIVAL_AXIS_LABEL) + report_theme + theme(legend.position = "top")
farm_risk_plot <- ggplot(farm_risk, aes(lon, lat, colour = probability_lost_pct, size = production_yield)) +
  geom_point(alpha = .82) + coord_quickmap() +
  scale_colour_viridis_c(option = "C", name = "Probability lost (%)") +
  scale_size_area(max_size = 6, name = "Production") + labs(title = "Farm-level probability of reaching the lost state", subtitle = "Conditional on the configured starting-location mode", x = "Longitude", y = "Latitude") + report_theme + theme(legend.position = "right")

plots <- list(spatial_cascade_loss_distribution = loss_plot,
  spatial_cascade_loss_by_production_type = type_loss_plot,
  spatial_cascade_survival_comparison = survival_plot,
  spatial_cascade_farm_risk_map = farm_risk_plot)
for (nm in names(plots)) ggsave(file.path(assets_dir, paste0(nm, ".svg")), plots[[nm]], device = grDevices::svg, width = 9, height = 6, bg = "white")
if (interactive()) invisible(lapply(plots, print))
invisible(results)
