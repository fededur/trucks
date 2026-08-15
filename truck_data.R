# Shared reproducible truck population used by both the Shiny application and
# the standalone cluster analysis. One row always represents one truck.
source("project_config.R", local = TRUE)

make_dummy_trucks <- function(n = 120, seed = 2026) {
  set.seed(seed)

  centres <- data.frame(
    city = c("Whangarei", "Auckland", "Hamilton", "Tauranga", "Rotorua",
             "Gisborne", "Napier", "New Plymouth", "Palmerston North",
             "Wellington", "Nelson", "Christchurch", "Dunedin",
             "Queenstown", "Invercargill"),
    island = c(rep("North Island", 10), rep("South Island", 5)),
    lat = c(-35.725, -36.849, -37.787, -37.687, -38.137,
            -38.662, -39.492, -39.055, -40.356, -41.286,
            -41.270, -43.532, -45.878, -45.031, -46.413),
    lon = c(174.323, 174.763, 175.279, 176.165, 176.249,
            178.018, 176.912, 174.076, 175.611, 174.776,
            173.284, 172.637, 170.503, 168.662, 168.354),
    weight = c(0.03, 0.25, 0.08, 0.06, 0.04,
               0.03, 0.04, 0.03, 0.04, 0.10,
               0.04, 0.12, 0.07, 0.03, 0.04)
  )
  group <- sample(seq_len(nrow(centres)), n, replace = TRUE,
                  prob = centres$weight)

  data.frame(
    truck_id = sprintf("TRK-%03d", seq_len(n)),
    island = centres$island[group],
    home_region = centres$city[group],
    latitude = centres$lat[group] + rnorm(n, 0, 0.10),
    longitude = centres$lon[group] + rnorm(n, 0, 0.13),
    cargo_t = round(runif(n, 4, 24), 1),
    truck_brand = sample(c("Fiat", "Scania"), n, replace = TRUE),
    truck_age = sample(c("Old", "New"), n, replace = TRUE,
                       prob = c(0.45, 0.55)),
    stringsAsFactors = FALSE
  )
}

# Use a user-supplied CSV when present; otherwise use the reproducible dummy
# population. Column aliases are controlled entirely by project_config.R.
load_truck_source <- function(path = file.path("data", "trucks_input.csv")) {
  if (file.exists(path)) {
    raw_data <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    normalise_truck_columns(raw_data)
  } else {
    make_dummy_trucks()
  }
}
