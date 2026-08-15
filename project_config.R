# Project data configuration --------------------------------------------------
# Edit this file when source column names or category values change. The app
# and analysis use stable internal names, so their simulation code does not
# need to be rewritten.

# Left side: internal field required by the project.
# Right side: column name in data/trucks_input.csv.
TRUCK_COLUMNS <- c(
  truck_id = "truck_id",
  island = "island",
  home_region = "home_region",
  latitude = "latitude",
  longitude = "longitude",
  cargo_t = "cargo_t",
  truck_brand = "truck_brand",
  truck_age = "truck_age"
)

# The numeric value is the speed multiplier applied to BASE_SPEED_KMH.
# Add or rename categories here, for example c(Vintage = 0.8, Standard = 1,
# Electric = 1.3). Every age/category present in the input data should appear.
AGE_SPEED_MULTIPLIERS <- c(Old = 1, New = 2)
BASE_SPEED_KMH <- 40

# Optional preferred brand colours. New brands not listed here automatically
# receive a distinct colour, so adding a brand never requires code changes.
BRAND_COLOURS <- c(Fiat = "#167d8d", Scania = "#e57b25")

# Optional stable island prefixes for cluster IDs. Unlisted island categories
# receive I1, I2, ... automatically.
ISLAND_PREFIXES <- c(`North Island` = "N", `South Island` = "S")

# Convert configured source columns to the stable names used internally.
normalise_truck_columns <- function(data, mapping = TRUCK_COLUMNS) {
  missing_columns <- setdiff(unname(mapping), names(data))
  if (length(missing_columns)) {
    stop("Input data is missing configured column(s): ",
         paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  selected <- data[, unname(mapping), drop = FALSE]
  names(selected) <- names(mapping)
  selected$truck_id <- as.character(selected$truck_id)
  selected$island <- as.character(selected$island)
  selected$home_region <- as.character(selected$home_region)
  selected$truck_brand <- as.character(selected$truck_brand)
  selected$truck_age <- as.character(selected$truck_age)
  selected$latitude <- as.numeric(selected$latitude)
  selected$longitude <- as.numeric(selected$longitude)
  selected$cargo_t <- as.numeric(selected$cargo_t)

  if (anyDuplicated(selected$truck_id)) stop("truck_id values must be unique.")
  if (anyNA(selected)) stop("Required truck fields cannot contain missing values.")
  unknown_age <- setdiff(unique(selected$truck_age), names(AGE_SPEED_MULTIPLIERS))
  if (length(unknown_age)) {
    stop("Add speed multipliers in project_config.R for: ",
         paste(unknown_age, collapse = ", "), call. = FALSE)
  }
  selected
}

category_colours <- function(categories, preferred = BRAND_COLOURS) {
  categories <- sort(unique(as.character(categories)))
  colours <- grDevices::hcl.colors(length(categories), palette = "Dark 3")
  names(colours) <- categories
  matching <- intersect(categories, names(preferred))
  colours[matching] <- preferred[matching]
  colours
}

island_prefixes <- function(categories) {
  categories <- sort(unique(as.character(categories)))
  prefixes <- paste0("I", seq_along(categories))
  names(prefixes) <- categories
  matching <- intersect(categories, names(ISLAND_PREFIXES))
  prefixes[matching] <- ISLAND_PREFIXES[matching]
  prefixes
}

# Restore configured source labels for human-facing tables and CSV exports.
# Internal RDS objects retain stable internal names for safe model operation.
display_truck_columns <- function(data) {
  matched <- intersect(names(TRUCK_COLUMNS), names(data))
  names(data)[match(matched, names(data))] <- unname(TRUCK_COLUMNS[matched])
  data
}
