# Dynamic Truck Recruitment Network ------------------------------------------
# A complete, self-contained Shiny application.
#
# Install the required packages once with:
# install.packages(c("shiny", "leaflet", "dplyr", "DT"))
#
# Run from this folder with:
# shiny::runApp()

library(shiny)
library(leaflet)
library(dplyr)
library(DT)

# ---- Dummy data -------------------------------------------------------------
# Keep data generation in one shared file so every analysis uses exactly the
# same population as the application.
source("truck_data.R", local = TRUE)
clustered_data_path <- file.path("data", "trucks_clustered.rds")
if (file.exists(clustered_data_path)) {
  trucks <- readRDS(clustered_data_path)
} else {
  trucks <- make_dummy_trucks()
  trucks$cluster <- NA_character_
  warning("Clustered data not found. Run truck_cluster_analysis.R to enable cluster mode.")
}

# ---- Geographic and simulation functions ----------------------------------
# Vectorised Haversine great-circle distance in kilometres. Using this instead
# of planar latitude/longitude differences keeps the geographic logic correct.
haversine_km <- function(lat1, lon1, lat2, lon2) {
  earth_radius_km <- 6371.0088
  to_rad <- pi / 180
  phi1 <- lat1 * to_rad
  phi2 <- lat2 * to_rad
  dphi <- (lat2 - lat1) * to_rad
  dlambda <- (lon2 - lon1) * to_rad
  a <- sin(dphi / 2)^2 + cos(phi1) * cos(phi2) * sin(dlambda / 2)^2
  2 * earth_radius_km * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

# Run one simulation. Keeping this as a plain R function makes it reusable for
# a future grid search over every start truck, radius, seed, or stopping rule.
run_network_simulation <- function(data, start_ids, radius_km, transfer_prob,
                                   old_speed_kmh, target_cargo_t,
                                   max_time_h, min_wave_cargo_t, seed) {
  set.seed(seed)
  n <- nrow(data)
  start_idx <- match(start_ids, data$truck_id)
  stopifnot(length(start_idx) > 0, !anyNA(start_idx))

  # Precompute all pairwise distances once. Cargo is deliberately absent from
  # this calculation: connections depend only on geographic proximity.
  distance_matrix <- matrix(0, n, n)
  for (i in seq_len(n)) {
    distance_matrix[i, ] <- haversine_km(
      data$latitude[i], data$longitude[i], data$latitude, data$longitude
    )
  }

  speed_multiplier <- unname(AGE_SPEED_MULTIPLIERS[data$truck_age])
  if (anyNA(speed_multiplier)) {
    stop("Missing age/category speed multiplier in project_config.R.")
  }
  speed <- old_speed_kmh * speed_multiplier
  recruited <- rep(FALSE, n)
  recruited[start_idx] <- TRUE
  frontier <- start_idx

  # Truck-level state. The starting truck is available at time zero and its
  # cargo is included immediately; it has no recruitment edge or travel.
  result <- data.frame(
    truck_id = data$truck_id[start_idx],
    island = data$island[start_idx],
    cluster = as.character(data$cluster[start_idx]),
    home_region = data$home_region[start_idx],
    truck_brand = data$truck_brand[start_idx],
    truck_age = data$truck_age[start_idx],
    cargo_t = data$cargo_t[start_idx],
    recruited_by = rep(NA_character_, length(start_idx)),
    wave = rep(0L, length(start_idx)),
    distance_km = rep(0, length(start_idx)),
    speed_kmh = speed[start_idx],
    travel_time = rep(0, length(start_idx)),
    recruitment_time = rep(0, length(start_idx)),
    cumulative_cargo = cumsum(data$cargo_t[start_idx]),
    stringsAsFactors = FALSE
  )
  edges <- data.frame(
    from = character(), to = character(), distance_km = numeric(),
    wave = integer(), stringsAsFactors = FALSE
  )
  wave_summary <- data.frame(
    wave = 0L, new_trucks = length(start_idx), new_cargo = sum(data$cargo_t[start_idx]),
    cumulative_trucks = length(start_idx), cumulative_cargo = sum(data$cargo_t[start_idx]),
    elapsed_time = 0
  )

  reason <- NULL
  wave <- 0L

  # Conditions that can be satisfied before the first expansion.
  if (sum(result$cargo_t) >= target_cargo_t) {
    reason <- "Target cargo reached by the starting trucks"
  } else if (max_time_h <= 0) {
    reason <- "Maximum simulation time reached"
  }

  while (is.null(reason) && length(frontier) > 0) {
    wave <- wave + 1L
    candidates <- which(!recruited)

    # Construct every geographically eligible offer from the current frontier.
    offers <- data.frame(parent = integer(), child = integer(),
                         distance_km = numeric(), stringsAsFactors = FALSE)
    for (parent in frontier) {
      # Water crossings are prohibited: even a large interaction radius cannot
      # connect a North Island truck to a South Island truck.
      eligible <- candidates[
        distance_matrix[parent, candidates] <= radius_km &
          data$island[candidates] == data$island[parent]
      ]
      if (length(eligible) > 0) {
        offers <- rbind(offers, data.frame(
          parent = rep(parent, length(eligible)),
          child = eligible,
          distance_km = distance_matrix[parent, eligible]
        ))
      }
    }

    if (nrow(offers) == 0) {
      reason <- "No additional trucks are within recruitment range"
      break
    }

    # Each eligible parent-child link receives one independent attempt. If more
    # than one parent succeeds, use the offer with the earliest arrival time.
    offers <- offers[runif(nrow(offers)) <= transfer_prob, , drop = FALSE]
    if (nrow(offers) == 0) {
      reason <- "No additional recruitment attempts succeeded"
      break
    }

    parent_times <- result$recruitment_time[
      match(data$truck_id[offers$parent], result$truck_id)
    ]
    offers$travel_time <- offers$distance_km / speed[offers$child]
    offers$arrival_time <- parent_times + offers$travel_time
    offers <- offers[order(offers$arrival_time, offers$distance_km,
                           data$truck_id[offers$parent]), , drop = FALSE]
    offers <- offers[!duplicated(offers$child), , drop = FALSE]

    # Offers arriving after the horizon do not recruit a truck. This keeps the
    # displayed elapsed time at or below the user's maximum.
    within_time <- offers$arrival_time <= max_time_h
    if (!any(within_time)) {
      reason <- "Maximum simulation time reached before the next recruitment"
      break
    }
    offers <- offers[within_time, , drop = FALSE]
    offers <- offers[order(offers$arrival_time, data$truck_id[offers$child]), ]

    new_rows <- data.frame(
      truck_id = data$truck_id[offers$child],
      island = data$island[offers$child],
      cluster = as.character(data$cluster[offers$child]),
      home_region = data$home_region[offers$child],
      truck_brand = data$truck_brand[offers$child],
      truck_age = data$truck_age[offers$child],
      cargo_t = data$cargo_t[offers$child],
      recruited_by = data$truck_id[offers$parent],
      wave = wave,
      distance_km = offers$distance_km,
      speed_kmh = speed[offers$child],
      travel_time = offers$travel_time,
      recruitment_time = offers$arrival_time,
      cumulative_cargo = NA_real_,
      stringsAsFactors = FALSE
    )
    # Cumulative cargo follows chronological recruitment order within a wave.
    new_rows$cumulative_cargo <- sum(result$cargo_t) + cumsum(new_rows$cargo_t)
    result <- rbind(result, new_rows)
    edges <- rbind(edges, data.frame(
      from = data$truck_id[offers$parent],
      to = data$truck_id[offers$child],
      distance_km = offers$distance_km,
      wave = wave,
      stringsAsFactors = FALSE
    ))
    recruited[offers$child] <- TRUE
    frontier <- offers$child

    new_cargo <- sum(new_rows$cargo_t)
    wave_summary <- rbind(wave_summary, data.frame(
      wave = wave,
      new_trucks = nrow(new_rows),
      new_cargo = new_cargo,
      cumulative_trucks = nrow(result),
      cumulative_cargo = sum(result$cargo_t),
      elapsed_time = max(result$recruitment_time)
    ))

    # Stop checks occur after completing a wave. Consequently target cargo may
    # be exceeded by the final wave; this is intentional for wave simultaneity.
    if (sum(result$cargo_t) >= target_cargo_t) {
      reason <- "Target cargo reached"
    } else if (max(result$recruitment_time) >= max_time_h) {
      reason <- "Maximum simulation time reached"
    } else if (new_cargo < min_wave_cargo_t) {
      reason <- sprintf(
        "Marginal cargo fell below threshold (%.1f < %.1f t)",
        new_cargo, min_wave_cargo_t
      )
    }
  }

  if (is.null(reason)) reason <- "No additional trucks can be recruited"

  list(
    trucks = result,
    waves = wave_summary,
    edges = edges,
    recruited_ids = result$truck_id,
    stop_reason = reason,
    total_distance_km = sum(result$distance_km),
    elapsed_time_h = max(result$recruitment_time),
    total_available_cargo = sum(data$cargo_t)
  )
}

# ---- User interface ---------------------------------------------------------
ui <- fluidPage(
  tags$head(
    tags$style(HTML("\n      body { background: #f3f6f9; color: #203040; }\n      .app-title { background: linear-gradient(110deg,#153a5b,#167d8d);\n        color:white; padding:20px 26px; margin:0 -15px 18px;\n        box-shadow:0 2px 9px rgba(0,0,0,.16); }\n      .app-title h2 { margin:0 0 4px; font-weight:700; }\n      .control-card, .panel-card, .kpi { background:white; border-radius:9px;\n        box-shadow:0 2px 8px rgba(27,54,78,.10); padding:16px; margin-bottom:16px; }\n      .control-card { border-top:4px solid #167d8d; }\n      .panel-card { min-height:120px; }\n      .kpi { border-left:5px solid #19a0ae; min-height:105px; }\n      .kpi-label { color:#6b7b88; text-transform:uppercase; font-size:11px;\n        letter-spacing:.7px; font-weight:700; }\n      .kpi-value { color:#173b57; font-size:23px; font-weight:700;\n        margin-top:8px; overflow-wrap:anywhere; }\n      .reason .kpi-value { font-size:15px; line-height:1.35; }\n      .btn-primary { background:#e57b25; border-color:#d66d18; font-weight:700; }\n      .nav-tabs>li.active>a { color:#167d8d; font-weight:700; }\n      .help-block { font-size:12px; }\n    "))
    , tags$style(HTML("\n      body { font-size:13px; }\n      .container-fluid { padding-left:10px; padding-right:10px; }\n      .app-title { padding:9px 18px; margin:0 -10px 7px; }\n      .app-title h2 { margin:0 0 1px; font-size:23px; }\n      .app-title div { font-size:11px; opacity:.9; }\n      .control-card, .panel-card, .kpi { border-radius:7px; padding:9px; margin-bottom:7px; }\n      .control-card { border-top-width:3px; padding:9px; }\n      .control-card h4 { margin:2px 0 7px; font-size:15px; }\n      .control-card .form-group { margin-bottom:7px; }\n      .control-card .control-label { margin-bottom:1px; font-size:11px; }\n      .control-card .well { padding:6px; margin-bottom:7px; }\n      .control-card details { border-top:1px solid #dce3e8; margin-top:5px; padding-top:5px; }\n      .control-card summary { color:#167d8d; cursor:pointer; font-weight:700; margin-bottom:6px; }\n      .panel-card { min-height:70px; }\n      .kpi-grid { display:grid; grid-template-columns:repeat(12,minmax(0,1fr)); gap:7px; }\n      .kpi-cell { grid-column:span 3; min-width:0; }\n      .kpi-cell:nth-child(5), .kpi-cell:nth-child(6), .kpi-cell:nth-child(7) { grid-column:span 2; }\n      .kpi-cell:nth-child(8) { grid-column:span 6; }\n      .kpi { border-left-width:4px; min-height:72px; height:calc(100% - 7px); }\n      .kpi-label { font-size:9px; letter-spacing:.4px; }\n      .kpi-value { font-size:16px; margin-top:2px; line-height:1.15; }\n      .kpi-value small { font-size:9px; font-weight:500; }\n      .reason .kpi-value { font-size:11px; line-height:1.2; }\n      .nav-tabs>li>a { padding:6px 9px; font-size:11px; }\n      .help-block { font-size:9px; margin:3px 0 0; }\n      #map { height:clamp(330px,calc(100vh - 275px),520px) !important; }\n      #progress { height:clamp(320px,calc(100vh - 285px),480px) !important; }\n      @media (max-width:1100px) { .kpi-cell { grid-column:span 6 !important; } }\n      @media (max-width:767px) {\n        .kpi-cell { grid-column:span 12 !important; }\n        #map, #progress { height:390px !important; }\n      }\n    "))
  ),
  div(class = "app-title",
      h2("Dynamic Truck Recruitment Network"),
      div("Stochastic, geographic, wave-by-wave transport network simulation")),
  sidebarLayout(
    sidebarPanel(width = 3, class = "control-card",
      h4("Simulation controls"),
      radioButtons(
        "analysis_level", "Start network from",
        choices = c("Individual trucks" = "truck", "Truck clusters" = "cluster"),
        selected = "truck", inline = TRUE
      ),
      conditionalPanel(
        condition = "input.analysis_level == 'truck'",
        selectInput("start_ids", "Starting truck(s)", choices = trucks$truck_id,
                    selected = trucks$truck_id[1], multiple = TRUE,
                    selectize = TRUE)
      ),
      conditionalPanel(
        condition = "input.analysis_level == 'cluster'",
        selectInput(
          "start_clusters", "Starting cluster(s)",
          choices = sort(unique(na.omit(as.character(trucks$cluster)))),
          selected = sort(unique(na.omit(as.character(trucks$cluster))))[1],
          multiple = TRUE, selectize = TRUE
        )
      ),
      sliderInput("radius", "Interaction radius (km)", min = 5, max = 500,
                  value = 100, step = 5),
      sliderInput("probability", "Information transfer probability",
                  min = 0, max = 1, value = 0.70, step = 0.05),
      tags$details(
        tags$summary("Speed and stopping parameters"),
        numericInput("old_speed", "Base truck speed (km/h)", value = BASE_SPEED_KMH,
                     min = 1, max = 150, step = 1),
        wellPanel(
          tags$small("Calculated speeds by category"),
          h4(textOutput("new_speed", inline = TRUE))
        ),
        numericInput("target_cargo", "Target cargo (tonnes)", value = 300,
                     min = 0, step = 10),
        numericInput("max_time", "Maximum simulation time (hours)", value = 2,
                     min = 0, step = 0.25),
        numericInput("min_marginal", "Minimum cargo per wave (tonnes)", value = 10,
                     min = 0, step = 1),
        numericInput("seed", "Random seed", value = 123, min = 0, step = 1)
      ),
      actionButton("run", "Run Simulation", class = "btn-primary btn-block",
                   icon = icon("play")),
      helpText("Inputs take effect when Run Simulation is clicked.")
    ),
    mainPanel(width = 9,
      div(class = "kpi-grid",
        div(class = "kpi-cell", uiOutput("kpi_start")),
        div(class = "kpi-cell", uiOutput("kpi_trucks")),
        div(class = "kpi-cell", uiOutput("kpi_cargo")),
        div(class = "kpi-cell", uiOutput("kpi_time")),
        div(class = "kpi-cell", uiOutput("kpi_distance")),
        div(class = "kpi-cell", uiOutput("kpi_available")),
        div(class = "kpi-cell", uiOutput("kpi_capture")),
        div(class = "kpi-cell", uiOutput("kpi_reason"))
      ),
      tabsetPanel(
        tabPanel("Network map", div(class = "panel-card", leafletOutput("map", height = 430))),
        tabPanel("Progress chart", div(class = "panel-card", plotOutput("progress", height = 410))),
        tabPanel("Truck results", div(class = "panel-card", DTOutput("truck_table"))),
        tabPanel("Brand summary", div(class = "panel-card", DTOutput("brand_table"))),
        tabPanel("Cluster results", div(class = "panel-card", DTOutput("cluster_table"))),
        tabPanel("Wave summary", div(class = "panel-card", DTOutput("wave_table"))),
        tabPanel("Assumptions",
          div(class = "panel-card",
            h3("Model assumptions"),
            tags$ol(
              tags$li("Geographic Haversine distance alone determines eligible links; cargo never affects connectivity."),
              tags$li("Recruitment is restricted to the same island; trucks never form North-to-South Island links."),
              tags$li("Each frontier truck makes one independent recruitment attempt per eligible unrecruited neighbour in that wave."),
              tags$li("When several trucks recruit the same neighbour, the successful offer with the earliest calculated arrival wins."),
              tags$li("A recruited truck becomes a recruiter in the next wave; there are no repeated attempts on failed links from an old frontier."),
              tags$li("Each truck age/category uses the speed multiplier configured in project_config.R."),
              tags$li("Recruitment travel time equals edge distance divided by the recruited truck's speed. Total distance is the sum of successful recruitment edges."),
              tags$li("All selected starting trucks are active at time zero, and their cargo is counted immediately."),
              tags$li("In cluster mode, all trucks belonging to the selected geographic clusters are starting trucks at time zero."),
              tags$li("Target and marginal rules are checked after a complete wave, so the last wave remains recruited and may overshoot the target."),
              tags$li("The random seed makes a stochastic run reproducible.")),
            h4("Future extension points"),
            p("The simulation engine is separate from Shiny and returns truck-, edge-, and wave-level data. It can be called for every starting truck and parameter combination, then ranked by cargo per elapsed hour, cargo per kilometre, or a multi-objective score."),
            p("A later model can separate communication delay from physical travel, route trucks to a depot, add road-network travel times, permit repeated contact attempts, impose capacities/costs, or optimize decisions dynamically.")
          )
        )
      )
    )
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {
  output$new_speed <- renderText({
    paste(sprintf("%s: %.0f km/h", names(AGE_SPEED_MULTIPLIERS),
                  input$old_speed * AGE_SPEED_MULTIPLIERS), collapse = " | ")
  })

  simulation <- eventReactive(input$run, {
    selected_start_ids <- if (identical(input$analysis_level, "cluster")) {
      trucks$truck_id[as.character(trucks$cluster) %in% input$start_clusters]
    } else {
      input$start_ids
    }
    validate(
      need(length(selected_start_ids) > 0,
           if (identical(input$analysis_level, "cluster"))
             "Select at least one starting cluster."
           else "Select at least one starting truck."),
      need(input$old_speed > 0, "Base truck speed must be positive."),
      need(input$target_cargo >= 0, "Target cargo cannot be negative."),
      need(input$max_time >= 0, "Maximum time cannot be negative."),
      need(input$min_marginal >= 0, "Marginal threshold cannot be negative.")
    )
    sim <- run_network_simulation(
      data = trucks,
      start_ids = selected_start_ids,
      radius_km = input$radius,
      transfer_prob = input$probability,
      old_speed_kmh = input$old_speed,
      target_cargo_t = input$target_cargo,
      max_time_h = input$max_time,
      min_wave_cargo_t = input$min_marginal,
      seed = as.integer(input$seed)
    )
    sim$selection_mode <- input$analysis_level
    sim$selected_clusters <- if (identical(input$analysis_level, "cluster"))
      input$start_clusters else character()
    sim
  }, ignoreNULL = FALSE)

  # A shared brand-level view keeps every output consistent and makes it easy
  # to add more brands later without changing the simulation algorithm.
  brand_summary <- reactive({
    sim <- simulation()
    brands <- sort(unique(trucks$truck_brand))
    data.frame(
      truck_brand = brands,
      available_trucks = sapply(brands, function(b) sum(trucks$truck_brand == b)),
      available_cargo_t = sapply(brands, function(b) sum(trucks$cargo_t[trucks$truck_brand == b])),
      recruited_trucks = sapply(brands, function(b) sum(sim$trucks$truck_brand == b)),
      recruited_cargo_t = sapply(brands, function(b) sum(sim$trucks$cargo_t[sim$trucks$truck_brand == b])),
      distance_km = sapply(brands, function(b) sum(sim$trucks$distance_km[sim$trucks$truck_brand == b])),
      elapsed_time_h = sapply(brands, function(b) {
        values <- sim$trucks$recruitment_time[sim$trucks$truck_brand == b]
        if (length(values)) max(values) else 0
      }),
      stringsAsFactors = FALSE
    ) %>%
      mutate(cargo_capture_pct = 100 * recruited_cargo_t / available_cargo_t)
  })

  brand_lines <- function(values, suffix = "") {
    summary <- brand_summary()
    tagList(lapply(seq_len(nrow(summary)), function(i) {
      tags$small(sprintf("%s: %s%s", summary$truck_brand[i], values[i], suffix),
                 style = "display:block; line-height:1.35;")
    }))
  }

  kpi <- function(label, value, extra_class = "") {
    div(class = paste("kpi", extra_class),
        div(class = "kpi-label", label), div(class = "kpi-value", value))
  }
  output$kpi_start <- renderUI({
    sim <- simulation()
    starts <- sim$trucks %>% filter(wave == 0)
    if (identical(sim$selection_mode, "cluster")) {
      kpi("Starting clusters", tagList(
        paste(sim$selected_clusters, collapse = ", "),
        tags$small(sprintf("%d member trucks active", nrow(starts)),
                   style = "display:block; line-height:1.35;")
      ))
    } else {
      kpi("Starting trucks", tagList(
        length(starts$truck_id),
        tags$small(paste(starts$truck_id, collapse = ", "),
                   style = "display:block; line-height:1.35;")
      ))
    }
  })
  output$kpi_trucks <- renderUI({
    bs <- brand_summary()
    kpi("Trucks recruited", tagList(
      nrow(simulation()$trucks), brand_lines(bs$recruited_trucks)
    ))
  })
  output$kpi_cargo <- renderUI({
    bs <- brand_summary()
    kpi("Cargo recruited", tagList(
      sprintf("%.1f t", sum(simulation()$trucks$cargo_t)),
      brand_lines(sprintf("%.1f", bs$recruited_cargo_t), " t")
    ))
  })
  output$kpi_time <- renderUI(kpi("Elapsed time", sprintf("%.2f h", simulation()$elapsed_time_h)))
  output$kpi_distance <- renderUI(kpi("Total distance", sprintf("%.1f km", simulation()$total_distance_km)))
  output$kpi_available <- renderUI({
    bs <- brand_summary()
    kpi("Total cargo available", tagList(
      sprintf("%.1f t", simulation()$total_available_cargo),
      brand_lines(sprintf("%.1f", bs$available_cargo_t), " t")
    ))
  })
  output$kpi_capture <- renderUI({
    bs <- brand_summary()
    pct <- 100 * sum(simulation()$trucks$cargo_t) / simulation()$total_available_cargo
    kpi("Available cargo reached", tagList(
      sprintf("%.1f%%", pct),
      lapply(seq_len(nrow(bs)), function(i) {
        tags$small(sprintf("%s: %.1f / %.1f t (%.1f%%)",
                           bs$truck_brand[i], bs$recruited_cargo_t[i],
                           bs$available_cargo_t[i], bs$cargo_capture_pct[i]),
                   style = "display:block; line-height:1.35;")
      })
    ))
  })
  output$kpi_reason <- renderUI(kpi("Network stopping reason", simulation()$stop_reason, "reason"))

  output$map <- renderLeaflet({
    sim <- simulation()
    map_data <- trucks %>%
      mutate(status = case_when(
        truck_id %in% sim$trucks$truck_id[sim$trucks$wave == 0] ~ "Starting truck",
        truck_id %in% sim$recruited_ids ~ "Recruited",
        TRUE ~ "Not recruited"
      ))
    palette <- colorFactor(c("#e4572e", "#168aad", "#aeb8c2"),
                           levels = c("Starting truck", "Recruited", "Not recruited"))

    # Leaflet radius controls the visual marker radius. Scaling by square root
    # makes marker AREA (rather than diameter) approximately proportional to
    # cargo, while retaining a practical 5-14 pixel range for map interaction.
    cargo_range <- range(map_data$cargo_t)
    if (diff(cargo_range) == 0) {
      map_data$marker_radius <- 9
    } else {
      scaled_cargo <- (map_data$cargo_t - cargo_range[1]) / diff(cargo_range)
      map_data$marker_radius <- 5 + 9 * sqrt(scaled_cargo)
    }
    # Give the selected starting truck a visible outline without changing its
    # cargo-driven bubble size.
    map_data$marker_weight <- ifelse(map_data$status == "Starting truck", 4, 2)

    m <- leaflet(map_data) %>%
      addProviderTiles(providers$CartoDB.Positron) %>%
      fitBounds(min(map_data$longitude), min(map_data$latitude),
                max(map_data$longitude), max(map_data$latitude))

    # Draw edges first so markers remain easy to click.
    if (nrow(sim$edges) > 0) {
      for (i in seq_len(nrow(sim$edges))) {
        from_row <- map_data[map_data$truck_id == sim$edges$from[i], ]
        to_row <- map_data[map_data$truck_id == sim$edges$to[i], ]
        m <- m %>% addPolylines(
          lng = c(from_row$longitude, to_row$longitude),
          lat = c(from_row$latitude, to_row$latitude),
          color = "#315c78", weight = 2, opacity = 0.65,
          label = sprintf("Wave %d: %s → %s (%.1f km)",
                          sim$edges$wave[i], sim$edges$from[i],
                          sim$edges$to[i], sim$edges$distance_km[i])
        )
      }
    }
    m %>%
      addCircleMarkers(
        lng = ~longitude, lat = ~latitude,
        radius = ~marker_radius,
        color = ~palette(status), fillColor = ~palette(status), fillOpacity = 0.9,
        weight = ~marker_weight,
        label = ~sprintf("%s | Cluster %s | %s | %s | %.1f t | %s",
                         truck_id, cluster, truck_brand, truck_age, cargo_t, status),
        popup = ~sprintf("<b>%s</b><br>Island: %s<br>Cluster: %s<br>Brand: %s<br>Age: %s<br>Cargo: %.1f t<br>Status: %s",
                         truck_id, island, cluster, truck_brand, truck_age,
                         cargo_t, status)
      ) %>%
      addLegend("bottomright", pal = palette, values = ~status,
                title = HTML("Network status<br><small>Bubble size = cargo</small>"),
                opacity = 1)
  })

  output$progress <- renderPlot({
    sim <- simulation()
    waves <- sim$waves$wave
    brands <- sort(unique(trucks$truck_brand))
    brand_colours <- category_colours(brands)
    brand_shapes <- rep(c(16, 17, 15, 18, 8, 4, 3, 7), length.out = length(brands))
    cargo_by_brand <- do.call(cbind, lapply(brands, function(brand) {
      sapply(waves, function(w) {
        sum(sim$trucks$cargo_t[
          sim$trucks$truck_brand == brand & sim$trucks$wave <= w
        ])
      })
    }))
    colnames(cargo_by_brand) <- brands
    total_cargo <- rowSums(cargo_by_brand)

    old_par <- par(mar = c(4.5, 4.5, 2.5, 4.5) + 0.1)
    on.exit(par(old_par))
    matplot(waves, cargo_by_brand, type = "o", pch = brand_shapes,
            lwd = 3, lty = 1, col = brand_colours[brands],
            xlab = "Recruitment wave", ylab = "Cumulative cargo (tonnes)",
            xaxt = "n", main = "Cumulative recruited cargo by truck brand",
            ylim = c(0, max(total_cargo) * 1.08))
    lines(waves, total_cargo, lwd = 2, lty = 2, col = "#384955")
    axis(1, at = waves)
    legend("topleft", c(brands, "All brands"),
           col = c(brand_colours[brands], "#384955"),
           pch = c(brand_shapes, NA),
           lty = c(rep(1, length(brands)), 2),
           lwd = c(rep(3, length(brands)), 2), bty = "n")
    grid(col = "#dce3e8")
  })

  output$truck_table <- renderDT({
    formatted <- simulation()$trucks %>%
      mutate(across(c(cargo_t, distance_km, speed_kmh, travel_time,
                      recruitment_time, cumulative_cargo), ~round(.x, 3)))
    display_truck_columns(formatted)
  }, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)

  output$brand_table <- renderDT({
    brand_summary() %>%
      mutate(across(where(is.numeric), ~round(.x, 2)))
  }, options = list(pageLength = 10, dom = "t", ordering = FALSE), rownames = FALSE)

  output$cluster_table <- renderDT({
    sim <- simulation()
    clusters <- sort(unique(as.character(trucks$cluster)))
    data.frame(
      cluster = clusters,
      island = sapply(clusters, function(x) unique(trucks$island[trucks$cluster == x])[1]),
      available_trucks = sapply(clusters, function(x) sum(trucks$cluster == x)),
      available_cargo_t = sapply(clusters, function(x) sum(trucks$cargo_t[trucks$cluster == x])),
      reached_trucks = sapply(clusters, function(x) sum(sim$trucks$cluster == x)),
      reached_cargo_t = sapply(clusters, function(x) sum(sim$trucks$cargo_t[sim$trucks$cluster == x])),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        cargo_reached_pct = 100 * reached_cargo_t / available_cargo_t,
        selected_start_cluster = cluster %in% sim$selected_clusters,
        across(where(is.numeric), ~round(.x, 2))
      )
  }, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)

  output$wave_table <- renderDT({
    sim <- simulation()
    result <- sim$waves
    for (brand in sort(unique(trucks$truck_brand))) {
      result[[paste0(brand, "_new_trucks")]] <- sapply(result$wave, function(w) {
        sum(sim$trucks$wave == w & sim$trucks$truck_brand == brand)
      })
      result[[paste0(brand, "_new_cargo_t")]] <- sapply(result$wave, function(w) {
        sum(sim$trucks$cargo_t[sim$trucks$wave == w & sim$trucks$truck_brand == brand])
      })
    }
    result %>% mutate(across(where(is.numeric), ~round(.x, 3)))
  }, options = list(pageLength = 12, dom = "tip"), rownames = FALSE)
}

shinyApp(ui, server)
