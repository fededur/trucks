# Create the reproducible example farm population used by the fixed-geography
# Monte Carlo model. Run this only when intentionally rebuilding data/farms.csv.
set.seed(2026)

centres <- data.frame(
  region = c("Northland", "Auckland", "Waikato", "Bay of Plenty",
             "Taranaki", "Manawatu", "Hawke's Bay", "Wellington",
             "Nelson", "Canterbury", "Otago", "Southland"),
  island = c(rep("North Island", 8), rep("South Island", 4)),
  lat = c(-35.72, -36.85, -37.79, -37.69, -39.06, -40.36,
          -39.49, -41.29, -41.27, -43.53, -45.88, -46.41),
  lon = c(174.32, 174.76, 175.28, 176.17, 174.08, 175.61,
          176.91, 174.78, 173.28, 172.64, 170.50, 168.35),
  weight = c(0.05, 0.12, 0.22, 0.08, 0.07, 0.10,
             0.06, 0.05, 0.04, 0.11, 0.06, 0.04)
)

population_size <- 150L
centre_id <- sample(seq_len(nrow(centres)), population_size, replace = TRUE,
                    prob = centres$weight)

farms <- data.frame(
  farm_id = sprintf("FARM-%03d", seq_len(population_size)),
  longitude = centres$lon[centre_id] + rnorm(population_size, 0, 0.12),
  latitude = centres$lat[centre_id] + rnorm(population_size, 0, 0.10),
  island = centres$island[centre_id],
  region = centres$region[centre_id],
  production_type = sample(c("Eggs", "Meat"), population_size, replace = TRUE,
                           prob = c(0.55, 0.45)),
  is_sheltered = sample(c(TRUE, FALSE), population_size, replace = TRUE,
                        prob = c(0.45, 0.55)),
  production_yield = pmax(0, rnorm(population_size, 1000, 150)),
  stringsAsFactors = FALSE
)

if (!dir.exists("data")) dir.create("data", recursive = TRUE)
write.csv(farms, file.path("data", "farms.csv"), row.names = FALSE)
message("Created data/farms.csv with ", nrow(farms), " fixed farms.")
