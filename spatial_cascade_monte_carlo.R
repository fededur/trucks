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
CULL_COMPARISON <- config$cull_response$compare_with_no_control %||% TRUE
CULL_COVERAGE <- config$cull_response$coverage %||% 0.8
CULL_RESPONSE_DELAY <- config$cull_response$response_delay_days %||% 1L

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
BOUNDARY_FILE <- config$map_options$boundary_file %||%
  "data/nz_regional_council_2025.geojson"
SHOW_REGION_BOUNDARIES <- config$map_options$show_region_boundaries %||% FALSE
if (!is.logical(SHOW_REGION_BOUNDARIES) || length(SHOW_REGION_BOUNDARIES) != 1L) {
  stop("map_options.show_region_boundaries must be true or false.")
}

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
if (!is.logical(CULL_COMPARISON) || length(CULL_COMPARISON) != 1L) {
  stop("cull_response.compare_with_no_control must be true or false.")
}
assert_number(CULL_COVERAGE, "cull_response.coverage", 0)
if (CULL_COVERAGE > 1) stop("cull_response.coverage must be at most 1.")
assert_number(CULL_RESPONSE_DELAY, "cull_response.response_delay_days", 0, TRUE)
if (!START_MODE %in% c("random", "each_farm")) stop("starting_location_mode must be random or each_farm.")

haversine_distance <- function(lon1, lat1, lon2, lat2) {
  r <- 6371.0088; rad <- pi / 180
  phi1 <- lat1 * rad; phi2 <- lat2 * rad
  dphi <- (lat2 - lat1) * rad; dlambda <- (lon2 - lon1) * rad
  a <- sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlambda / 2)^2
  2 * r * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

# Convert Polygon and MultiPolygon GeoJSON into a plain data frame for
# ggplot2. This keeps map drawing free of additional spatial dependencies.
read_geojson_polygons <- function(path) {
  if (!file.exists(path)) stop("Map boundary file not found: ", path)
  geo <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  rows <- list(); group_id <- 0L
  add_polygon <- function(polygon, region_name) {
    for (ring in polygon) {
      group_id <<- group_id + 1L
      xy <- do.call(rbind, ring)
      rows[[length(rows) + 1L]] <<- data.frame(
        lon = as.numeric(xy[, 1]), lat = as.numeric(xy[, 2]),
        group = group_id, region = region_name, stringsAsFactors = FALSE)
    }
  }
  for (feature in geo$features) {
    region_name <- unlist(feature$properties, use.names = FALSE)[1] %||% "Region"
    geometry <- feature$geometry
    if (identical(geometry$type, "Polygon")) {
      add_polygon(geometry$coordinates, region_name)
    } else if (identical(geometry$type, "MultiPolygon")) {
      for (polygon in geometry$coordinates) add_polygon(polygon, region_name)
    }
  }
  if (!length(rows)) stop("No polygon geometry found in: ", path)
  do.call(rbind, rows)
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
cull_type_loss_log <- matrix(NA_real_, MONTE_CARLO_RUNS, length(types),
                             dimnames = list(NULL, types))
cull_type_lost_count_log <- cull_type_loss_log
type_survival_log <- array(NA_real_, c(MONTE_CARLO_RUNS, length(types), 2L),
                           dimnames = list(NULL, types, groups))
cull_type_survival_log <- array(NA_real_, c(MONTE_CARLO_RUNS, length(types), 2L),
                                dimnames = list(NULL, types, groups))
group_survival_log <- matrix(NA_real_, MONTE_CARLO_RUNS, 2L, dimnames = list(NULL, groups))
cull_group_survival_log <- matrix(NA_real_, MONTE_CARLO_RUNS, 2L,
                                  dimnames = list(NULL, groups))
farm_start_count <- farm_affected_count <- farm_lost_count <- integer(POPULATION_SIZE)
cull_farm_affected_count <- cull_farm_lost_count <-
  cull_farm_culled_count <- integer(POPULATION_SIZE)
baseline_affected_count_log <- baseline_lost_count_log <- numeric(MONTE_CARLO_RUNS)
cull_total_loss_log <- cull_direct_loss_log <- cull_affected_count_log <-
  cull_lost_count_log <- cull_farm_count_log <- net_loss_avoided_log <-
  rep(NA_real_, MONTE_CARLO_RUNS)

simulate_run <- function(seeds, transmission_rolls, coverage_rolls,
                         apply_culling = FALSE) {
  state <- rep("A", POPULATION_SIZE)
  timer_b <- timer_c <- integer(POPULATION_SIZE)
  selected_for_cull <- culled <- rep(FALSE, POPULATION_SIZE)
  cull_due_day <- rep(NA_integer_, POPULATION_SIZE)
  state[seeds] <- "C"

  if (apply_culling) {
    selected_for_cull[seeds] <- coverage_rolls[seeds] <= CULL_COVERAGE
    selected_seeds <- seeds[selected_for_cull[seeds]]
    if (length(selected_seeds)) cull_due_day[selected_seeds] <- CULL_RESPONSE_DELAY
  }

  for (day in seq_len(SIM_DURATION)) {
    if (apply_culling) {
      due <- which(state == "C" & selected_for_cull & !culled &
                     !is.na(cull_due_day) & cull_due_day <= day)
      if (length(due)) {
        due <- due[order(cull_due_day[due], farms$id[due])]
        to_cull <- due
        culled[to_cull] <- TRUE
        state[to_cull] <- "D"
      }
    }

    active <- which(state == "C")
    targets <- which(state == "A")
    existing_b <- which(state == "B")
    existing_c <- active

    if (length(active) && length(targets)) {
      new_b <- integer()
      for (target in targets) {
        raw <- exp(-LAMBDA * distance_km[target, active])
        protection <- if (farms$is_sheltered[target]) 1 - PROT_EFFICIENCY else 1
        adjusted <- raw * protection * barrier_multiplier[target, active]
        final_probability <- 1 - prod(1 - pmin(1, pmax(0, adjusted)))
        if (transmission_rolls[day, target] <= final_probability) {
          new_b <- c(new_b, target)
        }
      }
      if (length(new_b)) {
        state[new_b] <- "B"
        timer_b[new_b] <- 0L
      }
    }

    if (length(existing_b)) {
      timer_b[existing_b] <- timer_b[existing_b] + 1L
      new_c <- existing_b[timer_b[existing_b] >= COOLDOWN_PHASE_B]
      if (length(new_c)) {
        state[new_c] <- "C"
        timer_c[new_c] <- 0L
        if (apply_culling) {
          selected_for_cull[new_c] <- coverage_rolls[new_c] <= CULL_COVERAGE
          selected <- new_c[selected_for_cull[new_c]]
          if (length(selected)) cull_due_day[selected] <- day + CULL_RESPONSE_DELAY
        }
      }
    }

    if (length(existing_c)) {
      timer_c[existing_c] <- timer_c[existing_c] + 1L
      new_d <- existing_c[timer_c[existing_c] >= COOLDOWN_PHASE_C]
      if (length(new_d)) state[new_d] <- "D"
    }
  }

  list(state = state, affected = state != "A", lost = state == "D",
       culled = culled)
}

for (sim in seq_len(MONTE_CARLO_RUNS)) {
  seeds <- start_sets[[sim]]
  farm_start_count[seeds] <- farm_start_count[seeds] + 1L
  transmission_rolls <- matrix(runif(SIM_DURATION * POPULATION_SIZE),
                               nrow = SIM_DURATION)
  coverage_rolls <- runif(POPULATION_SIZE)
  baseline <- simulate_run(seeds, transmission_rolls, coverage_rolls, FALSE)
  state <- baseline$state
  affected <- baseline$affected
  lost <- baseline$lost
  farm_affected_count <- farm_affected_count + affected
  farm_lost_count <- farm_lost_count + lost
  production_loss_log[sim] <- sum(farms$production_yield[lost])
  baseline_affected_count_log[sim] <- sum(affected)
  baseline_lost_count_log[sim] <- sum(lost)
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

  if (CULL_COMPARISON) {
    controlled <- simulate_run(seeds, transmission_rolls, coverage_rolls, TRUE)
    cull_farm_affected_count <- cull_farm_affected_count + controlled$affected
    cull_farm_lost_count <- cull_farm_lost_count + controlled$lost
    cull_farm_culled_count <- cull_farm_culled_count + controlled$culled
    cull_total_loss_log[sim] <- sum(farms$production_yield[controlled$lost])
    cull_direct_loss_log[sim] <- sum(farms$production_yield[controlled$culled])
    cull_affected_count_log[sim] <- sum(controlled$affected)
    cull_lost_count_log[sim] <- sum(controlled$lost)
    cull_farm_count_log[sim] <- sum(controlled$culled)
    net_loss_avoided_log[sim] <- production_loss_log[sim] - cull_total_loss_log[sim]
    for (type in types) {
      type_rows <- farms$production_type == type
      cull_type_loss_log[sim, type] <-
        sum(farms$production_yield[type_rows & controlled$lost])
      cull_type_lost_count_log[sim, type] <-
        sum(type_rows & controlled$lost)
      for (group in groups) {
        sheltered <- identical(group, SHELTERED_LABEL)
        rows <- type_rows & farms$is_sheltered == sheltered
        if (sum(rows)) {
          cull_type_survival_log[sim, type, group] <-
            100 * sum(rows & controlled$state == "A") / sum(rows)
        }
      }
    }
    for (group in groups) {
      rows <- farms$is_sheltered == identical(group, SHELTERED_LABEL)
      cull_group_survival_log[sim, group] <-
        100 * sum(rows & controlled$state == "A") / sum(rows)
    }
  }
}

average_type_survival <- apply(type_survival_log, c(2, 3), mean, na.rm = TRUE)
average_survival <- colMeans(group_survival_log, na.rm = TRUE)
cull_average_type_survival <- apply(cull_type_survival_log, c(2, 3), mean,
                                    na.rm = TRUE)
cull_average_survival <- colMeans(cull_group_survival_log, na.rm = TRUE)
assets_dir <- "assets"; if (!dir.exists(assets_dir)) dir.create(assets_dir, recursive = TRUE)

farm_risk <- transform(farms,
  times_selected_as_start = farm_start_count,
  probability_affected_pct = 100 * farm_affected_count / MONTE_CARLO_RUNS,
  probability_lost_pct = 100 * farm_lost_count / MONTE_CARLO_RUNS,
  expected_loss_contribution = production_yield * farm_lost_count / MONTE_CARLO_RUNS)
write.csv(farm_risk, file.path(assets_dir, "spatial_cascade_farm_risk.csv"), row.names = FALSE)

cull_farm_risk <- NULL
if (CULL_COMPARISON) {
  cull_farm_risk <- transform(farms,
    times_selected_as_start = farm_start_count,
    probability_affected_pct = 100 * cull_farm_affected_count / MONTE_CARLO_RUNS,
    probability_lost_pct = 100 * cull_farm_lost_count / MONTE_CARLO_RUNS,
    probability_culled_pct = 100 * cull_farm_culled_count / MONTE_CARLO_RUNS,
    expected_loss_contribution = production_yield * cull_farm_lost_count /
      MONTE_CARLO_RUNS)
  write.csv(cull_farm_risk,
            file.path(assets_dir, "spatial_cascade_cull_farm_risk.csv"),
            row.names = FALSE)
}

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

cull_comparison_summary <- NULL
if (CULL_COMPARISON) {
  cull_comparison_runs <- data.frame(
    run = seq_len(MONTE_CARLO_RUNS),
    baseline_total_loss = production_loss_log,
    cull_total_loss = cull_total_loss_log,
    direct_cull_loss = cull_direct_loss_log,
    net_loss_avoided = net_loss_avoided_log,
    baseline_farms_affected = baseline_affected_count_log,
    cull_farms_affected = cull_affected_count_log,
    baseline_farms_lost = baseline_lost_count_log,
    cull_farms_lost = cull_lost_count_log,
    farms_culled = cull_farm_count_log)
  write.csv(cull_comparison_runs,
            file.path(assets_dir, "spatial_cascade_cull_comparison.csv"),
            row.names = FALSE)

  cull_type_summary <- do.call(rbind, lapply(types, function(type) {
    data.frame(
      production_type = type,
      farms_in_population = sum(farms$production_type == type),
      farms_in_population_pct = 100 * mean(farms$production_type == type),
      average_farms_lost = mean(cull_type_lost_count_log[, type]),
      average_farms_lost_pct = 100 * mean(cull_type_lost_count_log[, type]) /
        sum(farms$production_type == type),
      mean_total_production_loss = mean(cull_type_loss_log[, type]),
      total_loss_p25 = unname(quantile(cull_type_loss_log[, type], .25)),
      total_loss_p75 = unname(quantile(cull_type_loss_log[, type], .75)))
  }))
  write.csv(cull_type_summary,
            file.path(assets_dir, "spatial_cascade_cull_production_type_summary.csv"),
            row.names = FALSE)

  cull_comparison_summary <- data.frame(
    runs = MONTE_CARLO_RUNS,
    coverage_pct = 100 * CULL_COVERAGE,
    response_delay_days = CULL_RESPONSE_DELAY,
    mean_baseline_total_loss = mean(production_loss_log),
    mean_cull_total_loss = mean(cull_total_loss_log),
    mean_direct_cull_loss = mean(cull_direct_loss_log),
    mean_net_loss_avoided = mean(net_loss_avoided_log),
    net_loss_reduction_pct = 100 * mean(net_loss_avoided_log) /
      mean(production_loss_log),
    runs_with_lower_total_loss_pct = 100 * mean(cull_total_loss_log < production_loss_log),
    mean_baseline_farms_affected = mean(baseline_affected_count_log),
    mean_cull_farms_affected = mean(cull_affected_count_log),
    mean_farms_culled = mean(cull_farm_count_log),
    net_loss_avoided_p05 = unname(quantile(net_loss_avoided_log, .05)),
    net_loss_avoided_p95 = unname(quantile(net_loss_avoided_log, .95)))
  write.csv(cull_comparison_summary,
            file.path(assets_dir, "spatial_cascade_cull_comparison_summary.csv"),
            row.names = FALSE)
}

results <- list(scenario = SCENARIO_NAME, configuration = config, farms = farms,
  population_size = POPULATION_SIZE, population_size_source = paste("fixed table", INPUT_FILE),
  starting_location_mode = START_MODE, monte_carlo_runs = MONTE_CARLO_RUNS,
  production_loss_log = production_loss_log, production_loss_by_type_log = type_loss_log,
  lost_count_by_type_log = type_lost_count_log, survival_by_type_and_shelter_log = type_survival_log,
  average_survival = average_survival, farm_risk_summary = farm_risk,
  cull_farm_risk_summary = cull_farm_risk,
  cull_comparison_summary = cull_comparison_summary,
  cull_total_loss_log = if (CULL_COMPARISON) cull_total_loss_log else NULL,
  cull_production_loss_by_type_log = if (CULL_COMPARISON) cull_type_loss_log else NULL,
  cull_lost_count_by_type_log = if (CULL_COMPARISON) cull_type_lost_count_log else NULL,
  cull_survival_by_type_and_shelter_log = if (CULL_COMPARISON) cull_type_survival_log else NULL,
  cull_average_survival = if (CULL_COMPARISON) cull_average_survival else NULL,
  direct_cull_loss_log = if (CULL_COMPARISON) cull_direct_loss_log else NULL,
  net_loss_avoided_log = if (CULL_COMPARISON) net_loss_avoided_log else NULL)
saveRDS(results, file.path(assets_dir, "spatial_cascade_results.rds"))

report_theme <- theme_minimal(base_size = 12) + theme(
  plot.title = element_text(face = "bold", colour = "#173B57"),
  plot.subtitle = element_text(colour = "#526575"), panel.grid.minor = element_blank(), legend.position = "none")
total_population_capacity <- sum(farms$production_yield)
make_loss_plot_rows <- function(loss_values, response_label) {
  ordered_loss <- sort(loss_values / total_population_capacity)
  rbind(
    data.frame(response = response_label,
               scenario_order = seq_along(ordered_loss),
               outcome = "Production remaining", share = 1 - ordered_loss),
    data.frame(response = response_label,
               scenario_order = seq_along(ordered_loss),
               outcome = "Production lost", share = ordered_loss))
}
loss_plot_data <- make_loss_plot_rows(production_loss_log, "No culling")
if (CULL_COMPARISON) {
  loss_plot_data <- rbind(loss_plot_data,
    make_loss_plot_rows(cull_total_loss_log, "Culling"))
}
loss_plot_data$response <- factor(loss_plot_data$response,
  levels = c("No culling", "Culling"))
loss_plot_data$outcome <- factor(loss_plot_data$outcome,
  levels = c("Production lost", "Production remaining"))
loss_plot <- ggplot(loss_plot_data,
                    aes(scenario_order, share, fill = outcome)) +
  geom_col(width = 1, colour = NA) +
  scale_fill_manual(values = c("Production lost" = "#E76F51",
                               "Production remaining" = "#D8E1E6"),
                    name = NULL) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .25),
                     labels = function(x) paste0(round(100 * x), "%"),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
  facet_wrap(~response, ncol = if (CULL_COMPARISON) 2 else 1) +
  labs(title = paste(SCENARIO_NAME, "production outcome across repeated scenarios"),
       subtitle = paste(MONTE_CARLO_RUNS, "matched scenarios per response;",
         "each panel is ordered from lower to higher percentage loss"),
       x = "Complete scenarios (ordered by loss, not time)",
       y = "Share of national production capacity") +
  report_theme + theme(legend.position = "top",
    strip.text = element_text(face = "bold", colour = "#173B57", size = 12),
    panel.spacing.x = grid::unit(1, "lines"))
type_loss_data <- as.data.frame(as.table(type_loss_log))
names(type_loss_data) <- c("run", "production_type", "loss")
type_loss_data$response <- "No culling"
if (CULL_COMPARISON) {
  cull_type_loss_data <- as.data.frame(as.table(cull_type_loss_log))
  names(cull_type_loss_data) <- c("run", "production_type", "loss")
  cull_type_loss_data$response <- "Culling"
  type_loss_data <- rbind(type_loss_data, cull_type_loss_data)
}
type_loss_data$response <- factor(type_loss_data$response,
                                  levels = c("No culling", "Culling"))
type_loss_plot <- ggplot(type_loss_data, aes(production_type, loss, fill = production_type)) +
  geom_boxplot(width = .62, outlier.alpha = .2) +
  facet_wrap(~response, nrow = 1) +
  scale_fill_manual(values = setNames(hcl.colors(length(types), "Dark 3"), types)) +
  labs(title = "Total loss by production type and response",
       subtitle = "Matched management scenarios; production type does not alter spread",
       x = PRODUCTION_TYPE_LABEL, y = VALUE_AXIS_LABEL) + report_theme +
  theme(strip.text = element_text(face = "bold", colour = "#173B57", size = 12))
survival_data <- as.data.frame(as.table(average_type_survival))
names(survival_data) <- c("production_type", "protection_group", "survival")
survival_data$response <- "No culling"
if (CULL_COMPARISON) {
  cull_survival_data <- as.data.frame(as.table(cull_average_type_survival))
  names(cull_survival_data) <- c("production_type", "protection_group", "survival")
  cull_survival_data$response <- "Culling"
  survival_data <- rbind(survival_data, cull_survival_data)
}
survival_data$response <- factor(survival_data$response,
                                 levels = c("No culling", "Culling"))
survival_plot <- ggplot(survival_data, aes(production_type, survival, fill = protection_group)) +
  geom_col(position = position_dodge(.72), width = .66) + geom_text(aes(label = sprintf("%.1f%%", survival)), position = position_dodge(.72), vjust = -.4) +
  facet_wrap(~response, nrow = 1) +
  scale_fill_manual(values = c("#3B9AB2", "#F28E2B"), name = "Protection") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(c(0, .04)),
                     labels = function(x) paste0(round(x), "%")) +
  labs(title = paste("Average final", PHASE_A_LABEL, "rate by response"),
       subtitle = "Matched management scenarios, split by production type and shelter",
       x = PRODUCTION_TYPE_LABEL, y = SURVIVAL_AXIS_LABEL) + report_theme +
  theme(legend.position = "top",
        strip.text = element_text(face = "bold", colour = "#173B57", size = 12))
nz_boundaries <- read_geojson_polygons(BOUNDARY_FILE)
boundary_line_colour <- if (SHOW_REGION_BOUNDARIES) "#87929A" else NA
map_risk_data <- if (CULL_COMPARISON) {
  rbind(transform(farm_risk, probability_culled_pct = 0,
                  response = "No control"),
        transform(cull_farm_risk, response = "Cull"))
} else {
  transform(farm_risk, probability_culled_pct = 0,
            response = "No control")
}
map_risk_data$response <- factor(map_risk_data$response,
                                 levels = c("No control", "Cull"))
farm_risk_plot <- ggplot() +
  geom_polygon(data = nz_boundaries, aes(lon, lat, group = group),
    fill = "#E5E7E9", colour = boundary_line_colour, linewidth = .25) +
  geom_point(data = map_risk_data,
    aes(lon, lat, colour = probability_lost_pct, size = production_yield),
    alpha = .82) +
  facet_wrap(~response, nrow = 1) +
  coord_quickmap(xlim = c(config$spatial_bounds$minimum_longitude,
                         config$spatial_bounds$maximum_longitude),
                 ylim = c(config$spatial_bounds$minimum_latitude,
                          config$spatial_bounds$maximum_latitude),
                 expand = FALSE) +
  scale_colour_viridis_c(option = "C", name = "Probability lost (%)") +
  scale_size_area(max_size = 6, name = "Production") +
  labs(title = "Farm exposure with and without culling",
       subtitle = "Matched scenarios and one shared probability scale",
       x = NULL, y = NULL) +
  report_theme +
  theme(legend.position = "right", axis.text = element_blank(),
        axis.ticks = element_blank())

if (CULL_COMPARISON) {
  cull_plot_data <- rbind(
    data.frame(response = "No culling", total_loss = production_loss_log),
    data.frame(response = "Culling", total_loss = cull_total_loss_log))
  cull_plot_data$response <- factor(cull_plot_data$response,
                                    levels = c("No culling", "Culling"))
  cull_comparison_plot <- ggplot(cull_plot_data,
      aes(x = "", y = total_loss, fill = response)) +
    geom_boxplot(width = .45, outlier.alpha = .18) +
    facet_wrap(~response, nrow = 1) +
    scale_fill_manual(values = c("No culling" = "#7B8790",
                                 "Culling" = "#167D8D")) +
    labs(title = "Total production loss by response",
         subtitle = paste0(CULL_COVERAGE * 100, "% coverage, ",
           CULL_RESPONSE_DELAY, "-day response delay"),
         x = NULL, y = VALUE_AXIS_LABEL) + report_theme +
    theme(strip.text = element_text(face = "bold", colour = "#173B57", size = 12),
          axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plots <- list(spatial_cascade_loss_distribution = loss_plot,
  spatial_cascade_loss_by_production_type = type_loss_plot,
  spatial_cascade_survival_comparison = survival_plot,
  spatial_cascade_farm_risk_map = farm_risk_plot)
if (CULL_COMPARISON) {
  plots$spatial_cascade_cull_comparison <- cull_comparison_plot
}
for (nm in names(plots)) {
  plot_width <- if (identical(nm, "spatial_cascade_farm_risk_map")) 13 else 9
  ggsave(file.path(assets_dir, paste0(nm, ".svg")), plots[[nm]],
         device = grDevices::svg, width = plot_width, height = 6, bg = "white")
}

# Appendix concept figure: timing controls when a farm can spread, while the
# distance-decay curve controls how likely it is to reach another farm.
timeline_end <- COOLDOWN_PHASE_B + COOLDOWN_PHASE_C +
  max(2, ceiling(COOLDOWN_PHASE_C * .3))
timeline_data <- rbind(
  data.frame(response = "No control", phase = "Unaffected", start = -2, end = 0),
  data.frame(response = "No control", phase = PHASE_B_LABEL, start = 0,
             end = COOLDOWN_PHASE_B),
  data.frame(response = "No control", phase = PHASE_C_LABEL,
             start = COOLDOWN_PHASE_B,
             end = COOLDOWN_PHASE_B + COOLDOWN_PHASE_C),
  data.frame(response = "No control", phase = PHASE_D_LABEL,
             start = COOLDOWN_PHASE_B + COOLDOWN_PHASE_C, end = timeline_end),
  data.frame(response = "Covered farm with culling", phase = "Unaffected",
             start = -2, end = 0),
  data.frame(response = "Covered farm with culling", phase = PHASE_B_LABEL,
             start = 0, end = COOLDOWN_PHASE_B),
  data.frame(response = "Covered farm with culling", phase = PHASE_C_LABEL,
             start = COOLDOWN_PHASE_B,
             end = COOLDOWN_PHASE_B + CULL_RESPONSE_DELAY),
  data.frame(response = "Covered farm with culling", phase = "Culled",
             start = COOLDOWN_PHASE_B + CULL_RESPONSE_DELAY,
             end = timeline_end))
timeline_data <- timeline_data[timeline_data$end > timeline_data$start, ]
timeline_data$mid <- (timeline_data$start + timeline_data$end) / 2
timeline_data$label <- ifelse(timeline_data$end - timeline_data$start >= 1.5,
                              timeline_data$phase, "")
timeline_data$label[timeline_data$phase == PHASE_D_LABEL] <- "Lost"
phase_colours <- setNames(
  c("#D8DEE2", "#F2B134", "#E76F51", "#7B8790", "#167D8D"),
  c("Unaffected", PHASE_B_LABEL, PHASE_C_LABEL, PHASE_D_LABEL, "Culled"))
timing_plot <- ggplot(timeline_data) +
  geom_rect(aes(xmin = start, xmax = end,
                ymin = as.numeric(factor(response)) - .34,
                ymax = as.numeric(factor(response)) + .34, fill = phase),
            colour = "white", linewidth = .5) +
  geom_text(aes(x = mid, y = as.numeric(factor(response)), label = label),
            size = 3, colour = "#203040") +
  scale_fill_manual(values = phase_colours, name = NULL) +
  scale_y_continuous(
    breaks = seq_along(levels(factor(timeline_data$response))),
    labels = levels(factor(timeline_data$response)), expand = expansion(.18)) +
  scale_x_continuous(breaks = unique(c(-2, 0, COOLDOWN_PHASE_B,
    COOLDOWN_PHASE_B + CULL_RESPONSE_DELAY,
    COOLDOWN_PHASE_B + COOLDOWN_PHASE_C, timeline_end))) +
  labs(title = "A. When a reached farm can spread",
       subtitle = paste0("Incubation: ", COOLDOWN_PHASE_B,
         " days; active period: ", COOLDOWN_PHASE_C,
         " days; culling delay: ", CULL_RESPONSE_DELAY, " day(s)"),
       x = "Days since the farm was reached", y = NULL) +
  report_theme + theme(legend.position = "none")

decay_distance <- seq(0, 300, by = 1)
decay_data <- rbind(
  data.frame(distance_km = decay_distance,
             probability = exp(-LAMBDA * decay_distance),
             condition = "Same island, not sheltered"),
  data.frame(distance_km = decay_distance,
             probability = exp(-LAMBDA * decay_distance) * (1 - PROT_EFFICIENCY),
             condition = "Same island, sheltered"),
  data.frame(distance_km = decay_distance,
             probability = exp(-LAMBDA * decay_distance) *
               CROSS_BARRIER_MULTIPLIER,
             condition = "Different island, not sheltered"))
decay_plot <- ggplot(decay_data,
    aes(distance_km, probability, colour = condition, linetype = condition)) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = c("Same island, not sheltered" = "#E76F51",
    "Same island, sheltered" = "#167D8D",
    "Different island, not sheltered" = "#7B8790"), name = NULL) +
  scale_linetype_manual(values = c("Same island, not sheltered" = "solid",
    "Same island, sheltered" = "dashed",
    "Different island, not sheltered" = "dotted"), name = NULL) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, 1)) +
  labs(title = "B. How distance reduces one transmission opportunity",
       subtitle = paste0("Spatial decay factor: ", LAMBDA,
         "; shelter protection: ", round(100 * PROT_EFFICIENCY),
         "%; cross-island multiplier: ", CROSS_BARRIER_MULTIPLIER),
       x = "Distance between farms (km)",
       y = "Modelled transmission probability") +
  report_theme + theme(legend.position = "bottom")

concept_path <- file.path(assets_dir, "spatial_cascade_timing_distance_concept.svg")
grDevices::svg(concept_path, width = 13, height = 5.8, bg = "white")
grid::grid.newpage()
concept_layout <- grid::grid.layout(nrow = 1, ncol = 2,
                                    widths = grid::unit(c(1, 1), "null"))
grid::pushViewport(grid::viewport(layout = concept_layout))
print(timing_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(decay_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
grid::popViewport()
grDevices::dev.off()
if (interactive()) invisible(lapply(plots, print))
invisible(results)
