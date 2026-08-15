# Standalone cluster analysis for the truck population -----------------------
#
# Run the whole file from the project root:
# source("truck_cluster_analysis.R")
#
# Required packages:
# install.packages(c("ggplot2", "dplyr", "cluster", "maps"))
#
# SVG figures and CSV results are written to the outputs/ folder.

required_packages <- c("ggplot2", "dplyr", "cluster", "maps")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Install required packages first: ",
       paste(missing_packages, collapse = ", "), call. = FALSE)
}

library(ggplot2)
library(dplyr)

# ---- User-adjustable analysis settings -------------------------------------
RANDOM_SEED <- 2026

# Choosing the cluster-count range -------------------------------------------
# For a fleet of roughly 100-200 trucks, testing 2-12 clusters is a sensible
# broad range. However, this analysis runs separately for each island, so the
# useful maximum should depend on how many trucks are present on that island.
#
# The settings below aim for an AVERAGE of at least 8 trucks per cluster and
# never test more than 12 clusters. Examples:
#   25 trucks  -> test k = 2:3
#   50 trucks  -> test k = 2:6
#   80 trucks  -> test k = 2:10
#   100-200 trucks -> test up to k = 12
#
# This is an operational safeguard against producing many tiny fuel-risk zones.
# It does not guarantee every final cluster has 8 trucks because geographic
# groups can be uneven. Always review cluster counts in the saved summaries.
MIN_CLUSTERS <- 2
MAX_CLUSTERS <- 12
MIN_TRUCKS_PER_CLUSTER <- 8
K_RANGE <- MIN_CLUSTERS:MAX_CLUSTERS

island_k_range <- function(number_of_trucks) {
  upper_k <- min(MAX_CLUSTERS,
                 floor(number_of_trucks / MIN_TRUCKS_PER_CLUSTER))
  if (upper_k < MIN_CLUSTERS) {
    stop(
      "An island needs at least ",
      MIN_CLUSTERS * MIN_TRUCKS_PER_CLUSTER,
      " trucks for the configured cluster-size rule. Reduce ",
      "MIN_TRUCKS_PER_CLUSTER if smaller islands must be analysed.",
      call. = FALSE
    )
  }
  seq.int(MIN_CLUSTERS, upper_k)
}

# The risk zones are geographic, so location is the default and cargo does not
# decide cluster membership. The number of geographic zones is still selected
# automatically. Set this to "all" only for a separate sensitivity test.
FEATURE_MODE <- "geography"

if (!FEATURE_MODE %in% c("automatic", "all", "geography")) {
  stop("FEATURE_MODE must be 'automatic', 'all', or 'geography'.")
}

# The data-generating function is shared with app.R, ensuring that this script
# analyses exactly the same reproducible truck population as the Shiny app.
source("truck_data.R", local = TRUE)
trucks <- load_truck_source()

output_dir <- "outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
data_dir <- "data"
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# ---- Automatic model selection ---------------------------------------------
# Standardisation prevents cargo tonnes or coordinate units from dominating
# merely because their numerical scales differ.
candidate_feature_sets <- switch(
  FEATURE_MODE,
  automatic = list(
    geography = c("latitude", "longitude"),
    geography_and_cargo = c("latitude", "longitude", "cargo_t")
  ),
  all = list(geography_and_cargo = c("latitude", "longitude", "cargo_t")),
  geography = list(geography = c("latitude", "longitude"))
)

set.seed(RANDOM_SEED)
model_results <- list()
result_index <- 1L

# Fit and select clusters separately within each island. This is a hard
# operational boundary: no cluster can contain both North and South trucks.
for (island_name in unique(trucks$island)) {
  row_ids <- which(trucks$island == island_name)
  island_trucks <- trucks[row_ids, , drop = FALSE]

  for (feature_name in names(candidate_feature_sets)) {
    features <- candidate_feature_sets[[feature_name]]
    x <- scale(island_trucks[, features, drop = FALSE])
    rownames(x) <- seq_len(nrow(island_trucks))

    if (identical(features, c("latitude", "longitude"))) {
      # Build actual great-circle distances in kilometres. PAM then groups
      # nearby trucks around real observed truck locations (medoids).
      to_rad <- pi / 180
      lat <- island_trucks$latitude * to_rad
      lon <- island_trucks$longitude * to_rad
      distance_matrix <- outer(seq_along(lat), seq_along(lat), Vectorize(function(i, j) {
        dlat <- lat[j] - lat[i]
        dlon <- lon[j] - lon[i]
        a <- sin(dlat / 2)^2 + cos(lat[i]) * cos(lat[j]) * sin(dlon / 2)^2
        6371.0088 * 2 * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
      }))
      distances <- as.dist(distance_matrix)
    } else {
      # Optional sensitivity mode only: standardised location plus cargo.
      distances <- dist(x)
    }

    for (k in island_k_range(nrow(island_trucks))) {
      fit <- cluster::pam(distances, k = k, diss = TRUE)
      cluster_assignment <- fit$clustering
      silhouette_values <- cluster::silhouette(cluster_assignment, distances)
      model_results[[result_index]] <- list(
        island = island_name,
        row_ids = row_ids,
        feature_set = feature_name,
        features = features,
        k = k,
        average_silhouette = mean(silhouette_values[, "sil_width"]),
        fit = fit,
        cluster_assignment = cluster_assignment,
        scaled_data = x,
        silhouette = silhouette_values
      )
      result_index <- result_index + 1L
    }
  }
}

selection_table <- bind_rows(lapply(model_results, function(x) {
  data.frame(
    island = x$island,
    feature_set = x$feature_set,
    features = paste(x$features, collapse = " + "),
    k = x$k,
    average_silhouette = x$average_silhouette
  )
})) %>%
  arrange(island, desc(average_silhouette), k)

# Select the best model independently for each island.
best_rows <- selection_table %>%
  group_by(island) %>%
  arrange(desc(average_silhouette), k, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

trucks$cluster <- NA_character_
trucks$cluster_number <- NA_integer_
trucks$silhouette_width <- NA_real_
best_models <- list()
cluster_prefixes <- island_prefixes(trucks$island)
for (i in seq_len(nrow(best_rows))) {
  chosen <- best_rows[i, ]
  best_index <- which(vapply(model_results, function(x) {
    x$island == chosen$island && x$feature_set == chosen$feature_set &&
      x$k == chosen$k
  }, logical(1)))[1]
  model <- model_results[[best_index]]
  prefix <- unname(cluster_prefixes[model$island])
  trucks$cluster[model$row_ids] <- paste0(prefix, model$cluster_assignment)
  trucks$cluster_number[model$row_ids] <- model$cluster_assignment
  local_order <- as.integer(rownames(model$silhouette))
  trucks$silhouette_width[model$row_ids[local_order]] <-
    model$silhouette[, "sil_width"]
  best_models[[model$island]] <- model
}
trucks$cluster <- factor(trucks$cluster)
trucks$cluster_number <- factor(trucks$cluster_number)

# ---- Variable relevance -----------------------------------------------------
# Relevance is the proportion of a standardised variable's total variation
# occurring between clusters (eta-squared). Higher values mean the variable
# distinguishes the selected clusters more strongly.
feature_relevance <- bind_rows(lapply(names(best_models), function(island_name) {
  rows <- trucks$island == island_name
  scaled <- as.data.frame(scale(trucks[rows, c("latitude", "longitude", "cargo_t")]))
  bind_rows(lapply(names(scaled), function(variable) {
    values <- scaled[[variable]]
    clusters <- droplevels(trucks$cluster[rows])
    overall_mean <- mean(values)
    cluster_means <- tapply(values, clusters, mean)
    cluster_sizes <- table(clusters)
    between_ss <- sum(cluster_sizes * (cluster_means - overall_mean)^2)
    total_ss <- sum((values - overall_mean)^2)
    data.frame(
      island = island_name,
      variable = variable,
      relevance = between_ss / total_ss,
      selected_for_clustering = variable %in% best_models[[island_name]]$features
    )
  }))
})) %>% arrange(island, desc(relevance))

# ---- Tables saved for reproducibility --------------------------------------
brand_cluster_summary <- trucks %>%
  group_by(island, cluster, truck_brand) %>%
  summarise(
    trucks = n(),
    total_cargo_t = sum(cargo_t),
    mean_cargo_t = mean(cargo_t),
    mean_silhouette = mean(silhouette_width),
    .groups = "drop"
  )

write.csv(display_truck_columns(trucks),
          file.path(output_dir, "truck_cluster_assignments.csv"),
          row.names = FALSE)
# Stable data object consumed by the Shiny app. Cluster IDs therefore remain
# fixed until this analysis is deliberately rerun.
saveRDS(trucks, file.path(data_dir, "trucks_clustered.rds"))
write.csv(display_truck_columns(trucks),
          file.path(data_dir, "trucks_clustered.csv"), row.names = FALSE)
write.csv(selection_table, file.path(output_dir, "cluster_model_selection.csv"),
          row.names = FALSE)
write.csv(feature_relevance, file.path(output_dir, "cluster_feature_relevance.csv"),
          row.names = FALSE)
write.csv(brand_cluster_summary,
          file.path(output_dir, "cluster_brand_summary.csv"), row.names = FALSE)

# ---- Plot styling -----------------------------------------------------------
brand_colours <- category_colours(trucks$truck_brand)
plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", colour = "#173b57"),
    plot.subtitle = element_text(colour = "#526575"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

save_svg <- function(plot, filename, width = 10, height = 7) {
  ggsave(
    filename = file.path(output_dir, filename), plot = plot,
    device = grDevices::svg, width = width, height = height, units = "in",
    bg = "white"
  )
}

# 1. New Zealand map: a convex-hull polygon shows the geographic territory of
# each cluster. Points remain visible for individual trucks; point colour shows
# brand and point area shows cargo. No truck labels are drawn.
nz_map <- ggplot2::map_data("nz")

# Convex hulls are a simple, transparent cluster boundary. Clusters with fewer
# than three unique locations cannot form a polygon and remain visible as
# points. Island grouping guarantees a polygon never spans Cook Strait.
cluster_hulls <- trucks %>%
  group_by(island, cluster) %>%
  filter(n_distinct(longitude, latitude) >= 3) %>%
  slice(chull(longitude, latitude)) %>%
  ungroup()

map_plot <- ggplot() +
  geom_polygon(
    data = nz_map,
    aes(long, lat, group = group),
    # Neutral land background only; grey is deliberately not part of the
    # cluster palette, so coloured overlays are unambiguously cluster areas.
    fill = "#d1d5d8", colour = "#687680", linewidth = 0.35
  ) +
  geom_polygon(
    data = cluster_hulls,
    aes(longitude, latitude, group = cluster, fill = cluster),
    colour = "#526575", linewidth = 0.45, alpha = 0.22
  ) +
  geom_point(
    data = trucks,
    aes(longitude, latitude, colour = truck_brand, size = cargo_t),
    alpha = 0.82, stroke = 0.8
  ) +
  scale_colour_manual(values = brand_colours, name = "Truck brand") +
  scale_fill_discrete(name = "Cluster area") +
  scale_size_area(max_size = 7, name = "Cargo (t)") +
  coord_quickmap(xlim = c(166, 179), ylim = c(-48, -34), expand = FALSE) +
  labs(
    title = "New Zealand truck clusters",
    subtitle = paste(vapply(best_models, function(model) sprintf(
      "%s: %d clusters using %s (silhouette %.3f)", model$island,
      model$k, paste(model$features, collapse = ", "),
      model$average_silhouette
    ), character(1)), collapse = " | "),
    x = "Longitude", y = "Latitude"
  ) +
  plot_theme
save_svg(map_plot, "truck_clusters_nz_map.svg", 9, 10)

# 2. Automatic selection diagnostic. The peak marks the chosen combination of
# feature set and cluster count.
selection_plot <- ggplot(selection_table,
                         aes(k, average_silhouette, colour = feature_set)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_point(
    data = best_rows, aes(k, average_silhouette),
    shape = 21, fill = "gold", colour = "#173b57", size = 5, stroke = 1.2,
    inherit.aes = FALSE
  ) +
  facet_wrap(~island) +
  scale_x_continuous(breaks = K_RANGE) +
  labs(
    title = "Automatic cluster-model selection",
    subtitle = "The highlighted model has the highest average silhouette score",
    x = "Number of clusters (k)", y = "Average silhouette score",
    colour = "Candidate features"
  ) +
  plot_theme
save_svg(selection_plot, "cluster_model_selection.svg", 9, 6)

# 3. Parameter relevance for interpreting what defines the selected clusters.
relevance_plot <- ggplot(
  feature_relevance,
  aes(reorder(variable, relevance), relevance, fill = selected_for_clustering)
) +
  geom_col(width = 0.65) +
  coord_flip() +
  facet_wrap(~island) +
  scale_fill_manual(values = c(`TRUE` = "#167d8d", `FALSE` = "#b9c3ca"),
                    labels = c(`TRUE` = "Selected", `FALSE` = "Not selected")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(
    title = "Variables defining the selected clusters",
    subtitle = "Relevance is the share of each variable's variation explained by cluster membership",
    x = NULL, y = "Between-cluster variation", fill = "Automatic model"
  ) +
  plot_theme
save_svg(relevance_plot, "cluster_feature_relevance.svg", 8, 5)

# 4. Cargo composition by cluster and brand.
cargo_plot <- ggplot(
  brand_cluster_summary,
  aes(cluster, total_cargo_t, fill = truck_brand)
) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = brand_colours, name = "Truck brand") +
  labs(
    title = "Cargo represented in each cluster by truck brand",
    subtitle = "Cargo is summed after cluster assignment; it is not an observation weight",
    x = "Cluster", y = "Total cargo (tonnes)"
  ) +
  plot_theme
save_svg(cargo_plot, "cluster_cargo_by_brand.svg", 9, 6)

selection_message <- paste(vapply(best_models, function(model) sprintf(
  "%s: k=%d, features=%s, silhouette=%.3f", model$island, model$k,
  paste(model$features, collapse = "+"), model$average_silhouette
), character(1)), collapse = "; ")
message("Cluster analysis complete. ", selection_message,
        ". Results saved in: ",
        normalizePath(output_dir, winslash = "/", mustWork = FALSE))
