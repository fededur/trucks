# New Zealand Truck Recruitment and Fuel-Risk Zones

An R project for exploring a truck transport network across New Zealand.

It contains two related analyses:

- a geographic cluster analysis that groups nearby trucks into local
  fuel-shortage risk zones; and
- a Shiny simulation that starts from selected trucks or clusters and models
  wave-by-wave recruitment of nearby trucks.

North and South Island trucks are always kept separate.

## Main files

- `app.R` — Shiny application
- `truck_cluster_analysis.R` — geographic cluster analysis and SVG outputs
- `truck_data.R` — reproducible dummy truck population or custom-data loader
- `project_config.R` — source column mappings, brand colours, and speed settings
- `truck_analysis_report.qmd` — plain-language Quarto report
- `truck_analysis_report.html` — rendered self-contained report
- `data/trucks_clustered.rds` — stable clustered data loaded by the app

## Required R packages

```r
install.packages(c(
  "shiny", "leaflet", "dplyr", "DT", "ggplot2", "cluster", "maps"
))
```

Quarto is required only to rebuild the HTML report.

## Run the project

Open `trucks.Rproj` in RStudio.

Rebuild geographic clusters and analysis outputs:

```r
source("truck_cluster_analysis.R")
```

Run the Shiny application:

```r
shiny::runApp()
```

Render the report:

```bash
quarto render truck_analysis_report.qmd
```

## Use another dataset

Save it as `data/trucks_input.csv`, then update `TRUCK_COLUMNS` in
`project_config.R` if its column names differ. The input file is ignored by Git
by default because operational fleet data may be private.

New truck brands are handled automatically. Age or truck-type speed categories
are configured through `AGE_SPEED_MULTIPLIERS` in `project_config.R`.

After changing the source data or configuration, rerun
`truck_cluster_analysis.R` before starting the Shiny app.

## Important limitation

The included data are generated examples. Distances are direct geographic
distances, not road-network travel times. Results are for scenario exploration,
not live fleet operations.
