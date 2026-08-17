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
                        cross_barrier_multiplier, control_strategy,
                        control_coverage, response_delay_days,
                        seed) {
  set.seed(seed)
  n <- nrow(farms)
  starts <- match(start_ids, farms$id)
  if (!length(starts) || anyNA(starts)) stop("Select at least one valid starting farm.")

  state <- rep("A", n)
  timer_b <- timer_c <- integer(n)
  reached_day <- lost_day <- rep(NA_integer_, n)
  active_day <- control_day <- rep(NA_integer_, n)
  parent <- rep(NA_integer_, n)
  selected_for_control <- controlled <- rep(FALSE, n)
  state[starts] <- "C"
  reached_day[starts] <- 0L
  active_day[starts] <- 0L
  if (control_strategy != "none") {
    selected_for_control[starts] <- runif(length(starts)) <= control_coverage
    control_day[starts[selected_for_control[starts]]] <- response_delay_days
  }

  same_group <- outer(as.character(farms$barrier_group),
                      as.character(farms$barrier_group), `==`)
  barrier <- ifelse(same_group, 1, cross_barrier_multiplier)

  for (day in seq_len(days)) {
    # Apply culling before today's transmission. It stops outgoing
    # transmission and immediately counts the farm's capacity as lost.
    due <- which(state == "C" & selected_for_control & !controlled &
                   !is.na(control_day) & control_day <= day)
    if (length(due)) {
      controlled[due] <- TRUE
      state[due] <- "D"
      lost_day[due] <- day
    }

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
          # sample() treats a single numeric value as 1:value. Sample the
          # position instead so one active farm remains one valid source.
          successful_source <- active[
            sample.int(length(active), 1L, prob = source_probability)
          ]
          state[target] <- "B"
          reached_day[target] <- day
          parent[target] <- successful_source
        }
      }
    }

    if (length(existing_b)) {
      timer_b[existing_b] <- timer_b[existing_b] + 1L
      new_c <- existing_b[timer_b[existing_b] >= incubation_days]
      if (length(new_c)) {
        state[new_c] <- "C"
        active_day[new_c] <- day
        if (control_strategy != "none") {
          selected_for_control[new_c] <- runif(length(new_c)) <= control_coverage
          selected <- new_c[selected_for_control[new_c]]
          if (length(selected)) control_day[selected] <- day + response_delay_days
        }
      }
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
  farm_result$active_day <- active_day
  farm_result$control_day <- control_day
  farm_result$selected_for_control <- selected_for_control
  farm_result$controlled <- controlled
  farm_result$control_strategy <- control_strategy
  farm_result$parent_index <- parent
  farm_result$parent_id <- ifelse(is.na(parent), NA_character_, farms$id[parent])
  farm_result$final_state <- state
  farm_result
}

farm_choices <- setNames(farms$id, paste(farms$id, farms$region,
                                          farms$production_type, sep = " | "))
default_lambda <- config$model_parameters$spatial_decay_factor
help_label <- function(text, note) {
  tags$span(text, class = "parameter-label", title = note)
}

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
    .parameter-label { cursor:help; border-bottom:1px dotted #6f818d; }
    .btn-primary { background:#e57b25; border-color:#d66d18; font-weight:700; }
    .kpi-grid { display:grid; grid-template-columns:repeat(5,1fr); gap:8px; }
    .loss-grid { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
      gap:8px; margin-bottom:8px; }
    .kpi { border-left:4px solid #19a0ae; min-height:70px; }
    .kpi-label { color:#6b7b88; text-transform:uppercase; font-size:9px;
      letter-spacing:.5px; font-weight:700; }
    .kpi-value { color:#173b57; font-size:18px; font-weight:700; margin-top:3px; }
    .kpi-detail { color:#687985; font-size:9px; font-weight:600; margin-top:2px; }
    .loss-grid .kpi { min-height:62px; margin-bottom:0; }
    .loss-grid .kpi:nth-child(2), .loss-grid .kpi:nth-child(3) {
      border-left-color:#E76F51; }
    #cascade_map { height:clamp(390px,calc(100vh - 285px),620px) !important; }
    .legend-note { font-size:10px; color:#687985; margin-top:5px; }
    .nav-tabs { margin-bottom:8px; }
    .nav-tabs>li>a { color:#31536a; font-weight:700; padding:7px 14px; }
    .nav-tabs>li.active>a { color:#167d8d; }
    .about-wrap { max-width:1000px; margin:0 auto 20px; }
    .about-card { background:white; border-radius:8px; padding:16px 20px;
      margin-bottom:10px; box-shadow:0 2px 8px rgba(27,54,78,.10); }
    .about-card h3 { color:#173b57; font-size:18px; margin:0 0 8px; }
    .about-card h4 { color:#167d8d; font-size:14px; margin:12px 0 5px; }
    .about-card p, .about-card li { line-height:1.5; }
    @media(max-width:900px) { .kpi-grid { grid-template-columns:repeat(2,1fr); }
      #cascade_map { height:480px !important; } }
  "))),
  div(class = "app-title",
      h2("Fixed-Farm Spatial Cascade Explorer"),
      div("Choose starting farms and inspect how a single loss chain develops")),
  tabsetPanel(id = "app_tab",
    tabPanel("Explore scenario",
      sidebarLayout(
    sidebarPanel(width = 3, class = "control-card",
      h4("Scenario controls"),
      selectInput("start_ids", help_label("Starting farm(s)",
        "The farm or farms where the scenario begins. Select more than one to test several starting points at the same time."), choices = farm_choices,
                  selected = farms$id[1], multiple = TRUE, selectize = TRUE),
      sliderInput("display_day", help_label("Display state on day",
        "Moves the map through the completed scenario. It changes the day shown, but does not rerun the model."), min = 0,
                  max = config$model_parameters$simulation_days,
                  value = config$model_parameters$simulation_days, step = 1),
      numericInput("duration", help_label("Scenario duration (days)",
        "The number of days modelled. A longer period gives a transmission chain more time to develop."),
                   config$model_parameters$simulation_days, min = 1, max = 365),
      numericInput("lambda", help_label("Spatial decay factor",
        "Controls how quickly transmission probability falls with distance. A larger value makes long-distance spread less likely."), default_lambda,
                   min = 0, max = 1, step = .005),
      numericInput("incubation", help_label("Incubation period (days)",
        "Days between a farm being reached and becoming able to spread to other farms."),
                   config$model_parameters$incubation_period_days, min = 1),
      numericInput("active_period", help_label("Active period (days)",
        "Days an active farm can spread before its production capacity moves to the lost state."),
                   config$model_parameters$active_period_days, min = 1),
      sliderInput("protection", help_label("Shelter protection efficiency",
        "Reduces the chance that a sheltered farm is reached. Zero gives no protection; one gives full protection under this model."), min = 0,
                  max = 1, value = config$housing_attributes$barn_protection_efficiency,
                  step = .05),
      sliderInput("island_multiplier", help_label("Cross-island transmission multiplier",
        "Adds an island barrier after distance is considered. Zero blocks cross-island spread; one applies no extra island penalty."),
                  min = 0, max = 1,
                  value = config$model_parameters$cross_barrier_transmission_multiplier,
                  step = .01),
      radioButtons("control_strategy", help_label("Spread-control policy",
        "None shows the baseline. Cull removes covered active farms from transmission but counts their production as lost."),
        choices = c("None" = "none", "Cull" = "cull"),
        selected = "none", inline = TRUE),
      conditionalPanel("input.control_strategy == 'cull'",
        sliderInput("control_coverage", help_label("Farms covered by policy",
          "The share of active farms the culling programme can reach. Selection is random within each scenario."), min = 0,
                    max = 1, value = .8, step = .05),
        numericInput("response_delay", help_label("Response delay after active (days)",
          "Days between a covered farm becoming active and being culled. A shorter delay stops onward spread sooner."),
                     value = 1, min = 0, max = 60, step = 1)
      ),
      numericInput("seed", help_label("Random seed",
        "Fixes the random choices in a run. Use the same seed when comparing policies so differences are easier to interpret."), config$simulation_controls$random_seed,
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
        div(class = "kpi", div(class = "kpi-label", "Farms not affected"),
            div(class = "kpi-value", textOutput("unaffected_kpi", inline = TRUE))),
        div(class = "kpi", div(class = "kpi-label", "Farms controlled"),
            div(class = "kpi-value", textOutput("controlled_kpi", inline = TRUE))),
        div(class = "kpi", div(class = "kpi-label", "Displayed day"),
            div(class = "kpi-value", textOutput("day_kpi", inline = TRUE)))
      ),
      uiOutput("loss_cards"),
      div(class = "panel-card", leafletOutput("cascade_map"))
    ))),
    tabPanel("About the app",
      div(class = "about-wrap",
        div(class = "about-card",
          h3("What this app is for"),
          p("This app shows how a harmful event could move between a fixed set of New Zealand farms. It helps you compare possible starting points, transmission settings and a culling response."),
          p("It is a scenario explorer, not a forecast. One run shows one possible chain under the settings you choose.")),
        div(class = "about-card",
          h3("How the model works"),
          p("Each dot is one farm from the project data. Its location, production type, shelter status and production capacity stay fixed."),
          tags$ul(
            tags$li(tags$b("Unaffected: "), "the farm has not been reached."),
            tags$li(tags$b("Incubating: "), "the farm has been reached but cannot yet spread the event."),
            tags$li(tags$b("Active: "), "the farm can spread the event to other farms."),
            tags$li(tags$b("Lost: "), "the farm's production capacity is counted as lost.")),
          p("Nearby farms are more likely to be reached than distant farms. Shelter lowers the chance that a farm is reached. North-to-South Island spread is controlled separately by the cross-island setting.")),
        div(class = "about-card",
          h3("How culling works"),
          p("Choose Cull to apply a response after a farm becomes active. A culled farm stops spreading immediately, but all of its production capacity is counted as lost."),
          tags$ul(
            tags$li(tags$b("Farms covered by policy "), "is the share of active farms the response programme can reach."),
            tags$li(tags$b("Response delay "), "is the number of days between a farm becoming active and being culled.")),
          p("A faster response can reduce onward spread. Wider coverage can reach more active farms. Both may also increase the production directly lost through culling, so compare the result with the None setting.")),
        div(class = "about-card",
          h3("How to use the app"),
          tags$ol(
            tags$li("Select one or more starting farms."),
            tags$li("Set the scenario duration and transmission assumptions."),
            tags$li("Choose None for a baseline, or Cull and set its coverage and delay."),
            tags$li("Choose a random seed. Reuse the same seed when comparing settings."),
            tags$li("Select Run scenario."),
            tags$li("Move the day slider to see how the chain develops.")),
          h4("Reading the map"),
          p("Red farms start the scenario, yellow farms join the transmission chain, and grey farms are not affected by the displayed day. Lines show the modelled link from one farm to the next. Larger dots have more production capacity."),
          h4("Reading the figures"),
          p("The cards show farms reached, farms lost, farms culled, total national capacity lost, and losses by production type. Loss includes farms that naturally reach the lost state and farms removed through culling.")),
        div(class = "about-card",
          h3("Important limits"),
          p("The result depends on the data and settings supplied. Straight-line distance does not describe every real pathway, such as animal movements, shared workers, vehicles or equipment. Use the app to compare scenarios and support discussion, not as a stand-alone operational decision."))
      )
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
      control_strategy = input$control_strategy,
      control_coverage = input$control_coverage,
      response_delay_days = as.integer(input$response_delay),
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
  output$unaffected_kpi <- renderText(
    sum(day_view()$display_group == "Not affected"))
  output$controlled_kpi <- renderText({
    result <- day_view()
    day <- min(input$display_day, input$duration)
    sum(result$controlled & !is.na(result$control_day) & result$control_day <= day)
  })
  output$day_kpi <- renderText(paste("Day", min(input$display_day, input$duration)))

  output$loss_cards <- renderUI({
    result <- day_view()
    total_capacity <- sum(result$production_yield)
    total_lost <- sum(result$production_yield[result$lost_by_day])
    total_lost_pct <- if (total_capacity > 0) 100 * total_lost / total_capacity else 0

    card <- function(label, value, detail = NULL) {
      div(class = "kpi",
          div(class = "kpi-label", label),
          div(class = "kpi-value", value),
          if (!is.null(detail)) div(class = "kpi-detail", detail))
    }

    type_cards <- lapply(sort(unique(result$production_type)), function(type) {
      rows <- result$production_type == type
      type_capacity <- sum(result$production_yield[rows])
      type_lost <- sum(result$production_yield[rows & result$lost_by_day])
      type_pct <- if (type_capacity > 0) 100 * type_lost / type_capacity else 0
      card(paste(type, "capacity lost"),
           format(round(type_lost), big.mark = ","),
           sprintf("%.1f%% of %s capacity", type_pct, type))
    })

    div(class = "loss-grid",
        card("National production capacity",
             format(round(total_capacity), big.mark = ","),
             paste(format(nrow(result), big.mark = ","), "farms")),
        card("National capacity lost",
             format(round(total_lost), big.mark = ","),
             "Absolute production capacity"),
        card("Lost / national capacity",
             sprintf("%.1f%%", total_lost_pct),
             "Share of national production capacity"),
        type_cards)
  })

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
      "<b>%s</b><br>%s<br>%s<br>Production: %s<br>Sheltered: %s<br>Reached: %s<br>Control: %s<br>Lost: %s",
      result$id, result$region, result$production_type,
      format(round(result$production_yield), big.mark = ","),
      ifelse(result$is_sheltered, "Yes", "No"),
      ifelse(is.na(result$reached_day), "Not reached", paste("Day", result$reached_day)),
      ifelse(is.na(result$control_day), "Not selected",
             paste(tools::toTitleCase(result$control_strategy), "on day", result$control_day)),
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
