# Interactive Fixed-Farm Spatial Cascade Explorer ----------------------------
# Run from the project root with:
# shiny::runApp("spatial_cascade_app.R")

library(shiny)
library(leaflet)

config <- jsonlite::fromJSON("config.json", simplifyVector = TRUE)
`%||%` <- function(x, fallback) if (is.null(x) || !length(x) || is.na(x[1])) fallback else x

column_map <- unlist(config$input_data$column_mappings, use.names = TRUE)
required_fields <- c("id", "lon", "lat", "barrier_group", "region",
                     "production_type", "is_sheltered", "production_yield")
raw_farms <- read.csv(config$input_data$file, stringsAsFactors = FALSE,
                      check.names = FALSE)
farms <- raw_farms[, unname(column_map[required_fields]), drop = FALSE]
names(farms) <- required_fields
farms$id <- as.character(farms$id)
farms$lon <- as.numeric(farms$lon)
farms$lat <- as.numeric(farms$lat)
farms$production_yield <- as.numeric(farms$production_yield)
farms$is_sheltered <- as.logical(farms$is_sheltered)

haversine_distance <- function(lon1, lat1, lon2, lat2) {
  radius <- 6371.0088
  radians <- pi / 180
  phi1 <- lat1 * radians
  phi2 <- lat2 * radians
  delta_phi <- (lat2 - lat1) * radians
  delta_lon <- (lon2 - lon1) * radians
  a <- sin(delta_phi / 2)^2 + cos(phi1) * cos(phi2) * sin(delta_lon / 2)^2
  2 * radius * atan2(sqrt(a), sqrt(pmax(0, 1 - a)))
}

distance_matrix <- matrix(haversine_distance(
  rep(farms$lon, times = nrow(farms)), rep(farms$lat, times = nrow(farms)),
  rep(farms$lon, each = nrow(farms)), rep(farms$lat, each = nrow(farms))
), nrow = nrow(farms))

read_boundary_rings <- function(path) {
  if (!file.exists(path)) return(list())
  geo <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  rings <- list()
  add_polygon <- function(polygon, region) {
    for (ring in polygon) {
      xy <- do.call(rbind, ring)
      rings[[length(rings) + 1L]] <<- list(
        lng = as.numeric(xy[, 1]), lat = as.numeric(xy[, 2]), region = region)
    }
  }
  for (feature in geo$features) {
    region <- unlist(feature$properties, use.names = FALSE)[1] %||% "Region"
    if (identical(feature$geometry$type, "Polygon")) {
      add_polygon(feature$geometry$coordinates, region)
    } else if (identical(feature$geometry$type, "MultiPolygon")) {
      for (polygon in feature$geometry$coordinates) add_polygon(polygon, region)
    }
  }
  rings
}

boundary_rings <- read_boundary_rings(
  config$map_options$boundary_file %||% "data/nz_regional_council_2025.geojson")

run_cascade <- function(start_ids, days, lambda, incubation_days,
                        active_days, shelter_efficiency,
                        cross_barrier_multiplier, seed) {
  set.seed(seed)
  n <- nrow(farms)
  starts <- match(start_ids, farms$id)
  if (!length(starts) || anyNA(starts)) stop("Select at least one valid starting farm.")

  state <- rep("A", n)
  timer_b <- timer_c <- integer(n)
  reached_day <- lost_day <- rep(NA_integer_, n)
  parent <- rep(NA_integer_, n)
  state[starts] <- "C"
  reached_day[starts] <- 0L

  same_group <- outer(as.character(farms$barrier_group),
                      as.character(farms$barrier_group), `==`)
  barrier <- ifelse(same_group, 1, cross_barrier_multiplier)

  for (day in seq_len(days)) {
    active <- which(state == "C")
    targets <- which(state == "A")
    existing_b <- which(state == "B")
    existing_c <- active

    if (length(active) && length(targets)) {
      for (target in targets) {
        raw <- exp(-lambda * distance_matrix[target, active])
        protection <- if (farms$is_sheltered[target]) 1 - shelter_efficiency else 1
        source_probability <- pmin(1, pmax(0,
          raw * protection * barrier[target, active]))
        combined_probability <- 1 - prod(1 - source_probability)
        if (runif(1) <= combined_probability) {
          successful_source <- sample(active, 1L, prob = source_probability)
          state[target] <- "B"
          reached_day[target] <- day
          parent[target] <- successful_source
        }
      }
    }

    if (length(existing_b)) {
      timer_b[existing_b] <- timer_b[existing_b] + 1L
      new_c <- existing_b[timer_b[existing_b] >= incubation_days]
      if (length(new_c)) state[new_c] <- "C"
    }
    if (length(existing_c)) {
      timer_c[existing_c] <- timer_c[existing_c] + 1L
      new_d <- existing_c[timer_c[existing_c] >= active_days]
      if (length(new_d)) {
        state[new_d] <- "D"
        lost_day[new_d] <- day
      }
    }
  }

  farm_result <- farms
  farm_result$is_start <- seq_len(n) %in% starts
  farm_result$reached_day <- reached_day
  farm_result$lost_day <- lost_day
  farm_result$parent_index <- parent
  farm_result$parent_id <- ifelse(is.na(parent), NA_character_, farms$id[parent])
  farm_result$final_state <- state
  farm_result
}

farm_choices <- setNames(farms$id, paste(farms$id, farms$region,
                                          farms$production_type, sep = " | "))
default_lambda <- config$model_parameters$spatial_decay_factor

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { background:#f3f6f9; color:#203040; font-size:13px; }
    .container-fluid { padding:0 10px; }
    .app-title { background:linear-gradient(110deg,#153a5b,#167d8d);
      color:white; padding:10px 18px; margin:0 -10px 8px;
      box-shadow:0 2px 9px rgba(0,0,0,.16); }
    .app-title h2 { margin:0 0 2px; font-size:23px; font-weight:700; }
    .app-title div { font-size:11px; opacity:.9; }
    .control-card, .panel-card, .kpi { background:white; border-radius:7px;
      box-shadow:0 2px 8px rgba(27,54,78,.10); padding:10px; margin-bottom:8px; }
    .control-card { border-top:3px solid #167d8d; }
    .control-card h4 { margin:2px 0 8px; font-size:15px; }
    .control-card .form-group { margin-bottom:8px; }
    .control-card .control-label { font-size:11px; margin-bottom:2px; }
    .btn-primary { background:#e57b25; border-color:#d66d18; font-weight:700; }
    .kpi-grid { display:grid; grid-template-columns:repeat(4,1fr); gap:8px; }
    .kpi { border-left:4px solid #19a0ae; min-height:70px; }
    .kpi-label { color:#6b7b88; text-transform:uppercase; font-size:9px;
      letter-spacing:.5px; font-weight:700; }
    .kpi-value { color:#173b57; font-size:18px; font-weight:700; margin-top:3px; }
    #cascade_map { height:clamp(430px,calc(100vh - 210px),680px) !important; }
    .legend-note { font-size:10px; color:#687985; margin-top:5px; }
    @media(max-width:900px) { .kpi-grid { grid-template-columns:repeat(2,1fr); }
      #cascade_map { height:480px !important; } }
  "))),
  div(class = "app-title",
      h2("Fixed-Farm Spatial Cascade Explorer"),
      div("Choose starting farms and inspect how a single loss chain develops")),
  sidebarLayout(
    sidebarPanel(width = 3, class = "control-card",
      h4("Scenario controls"),
      selectInput("start_ids", "Starting farm(s)", choices = farm_choices,
                  selected = farms$id[1], multiple = TRUE, selectize = TRUE),
      sliderInput("display_day", "Display state on day", min = 0,
                  max = config$model_parameters$simulation_days,
                  value = config$model_parameters$simulation_days, step = 1),
      numericInput("duration", "Scenario duration (days)",
                   config$model_parameters$simulation_days, min = 1, max = 365),
      numericInput("lambda", "Spatial decay factor", default_lambda,
                   min = 0, max = 1, step = .005),
      numericInput("incubation", "Incubation period (days)",
                   config$model_parameters$incubation_period_days, min = 1),
      numericInput("active_period", "Active period (days)",
                   config$model_parameters$active_period_days, min = 1),
      sliderInput("protection", "Shelter protection efficiency", min = 0,
                  max = 1, value = config$housing_attributes$barn_protection_efficiency,
                  step = .05),
      sliderInput("island_multiplier", "Cross-island transmission multiplier",
                  min = 0, max = 1,
                  value = config$model_parameters$cross_barrier_transmission_multiplier,
                  step = .01),
      numericInput("seed", "Random seed", config$simulation_controls$random_seed,
                   min = 1, step = 1),
      actionButton("run", "Run scenario", class = "btn-primary btn-block"),
      div(class = "legend-note",
          "Higher spatial decay means transmission falls faster with distance. A cross-island multiplier of 0 blocks island crossings.")
    ),
    mainPanel(width = 9,
      div(class = "kpi-grid",
        div(class = "kpi", div(class = "kpi-label", "Farms reached"),
            div(class = "kpi-value", textOutput("reached_kpi", inline = TRUE))),
        div(class = "kpi", div(class = "kpi-label", "Farms in lost state"),
            div(class = "kpi-value", textOutput("lost_kpi", inline = TRUE))),
        div(class = "kpi", div(class = "kpi-label", "Production capacity lost"),
            div(class = "kpi-value", textOutput("loss_kpi", inline = TRUE))),
        div(class = "kpi", div(class = "kpi-label", "Displayed day"),
            div(class = "kpi-value", textOutput("day_kpi", inline = TRUE)))
      ),
      div(class = "panel-card", leafletOutput("cascade_map"))
    )
  )
)

server <- function(input, output, session) {
  scenario <- eventReactive(input$run, {
    validate(need(length(input$start_ids) > 0, "Select at least one starting farm."))
    validate(need(input$lambda >= 0, "Spatial decay must be non-negative."))
    updateSliderInput(session, "display_day", max = input$duration,
                      value = input$duration)
    run_cascade(
      start_ids = input$start_ids, days = as.integer(input$duration),
      lambda = input$lambda, incubation_days = as.integer(input$incubation),
      active_days = as.integer(input$active_period),
      shelter_efficiency = input$protection,
      cross_barrier_multiplier = input$island_multiplier,
      seed = as.integer(input$seed)
    )
  }, ignoreNULL = FALSE)

  day_view <- reactive({
    result <- scenario()
    day <- min(input$display_day, input$duration)
    result$display_group <- ifelse(
      result$is_start, "Initiating farm",
      ifelse(!is.na(result$reached_day) & result$reached_day <= day,
             "In transmission chain", "Not affected"))
    result$lost_by_day <- !is.na(result$lost_day) & result$lost_day <= day
    result
  })

  output$reached_kpi <- renderText({
    result <- day_view()
    sum(result$display_group != "Not affected")
  })
  output$lost_kpi <- renderText(sum(day_view()$lost_by_day))
  output$loss_kpi <- renderText({
    format(round(sum(day_view()$production_yield[day_view()$lost_by_day])),
           big.mark = ",")
  })
  output$day_kpi <- renderText(paste("Day", min(input$display_day, input$duration)))

  output$cascade_map <- renderLeaflet({
    result <- day_view()
    day <- min(input$display_day, input$duration)
    colours <- c("Initiating farm" = "#E76F51",
                 "In transmission chain" = "#F2B134",
                 "Not affected" = "#7B8790")
    map <- leaflet(options = leafletOptions(zoomControl = TRUE,
                                             minZoom = 4)) |>
      addProviderTiles(providers$CartoDB.PositronNoLabels)

    show_regions <- isTRUE(config$map_options$show_region_boundaries)
    for (ring in boundary_rings) {
      map <- addPolygons(map, lng = ring$lng, lat = ring$lat,
        fillColor = "#E5E7E9", fillOpacity = .55,
        color = if (show_regions) "#87929A" else "#A6AFB5",
        weight = if (show_regions) .7 else .4, opacity = .8,
        group = "New Zealand boundaries", options = pathOptions(interactive = FALSE))
    }

    chain <- result[!is.na(result$parent_index) & result$reached_day <= day, ]
    if (nrow(chain)) {
      for (i in seq_len(nrow(chain))) {
        parent <- result[chain$parent_index[i], ]
        map <- addPolylines(map,
          lng = c(parent$lon, chain$lon[i]), lat = c(parent$lat, chain$lat[i]),
          color = "#C78A1B", weight = 1.5, opacity = .6,
          options = pathOptions(interactive = FALSE))
      }
    }

    popup <- sprintf(
      "<b>%s</b><br>%s<br>%s<br>Production: %s<br>Sheltered: %s<br>Reached: %s<br>Lost: %s",
      result$id, result$region, result$production_type,
      format(round(result$production_yield), big.mark = ","),
      ifelse(result$is_sheltered, "Yes", "No"),
      ifelse(is.na(result$reached_day), "Not reached", paste("Day", result$reached_day)),
      ifelse(is.na(result$lost_day), "Not lost", paste("Day", result$lost_day)))
    radius <- 5 + 5 * sqrt(result$production_yield / max(result$production_yield))
    map |>
      addCircleMarkers(lng = result$lon, lat = result$lat,
        radius = radius, color = "white", weight = 1,
        fillColor = unname(colours[result$display_group]), fillOpacity = .9,
        popup = popup, label = result$id) |>
      addLegend(position = "bottomright", colors = colours,
        labels = names(colours), title = paste("State on day", day), opacity = 1) |>
      fitBounds(lng1 = config$spatial_bounds$minimum_longitude,
                lat1 = config$spatial_bounds$minimum_latitude,
                lng2 = config$spatial_bounds$maximum_longitude,
                lat2 = config$spatial_bounds$maximum_latitude)
  })
}

shinyApp(ui, server)
