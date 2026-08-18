# New Zealand Truck Recruitment and Fuel-Risk Zones

An R project for exploring a truck transport network across New Zealand.

It contains two related analyses:

- a geographic cluster analysis that groups nearby trucks into local
  fuel-shortage risk zones; and
- a Shiny simulation that starts from selected trucks or clusters and models
  wave-by-wave recruitment of nearby trucks.

It also includes an independent, domain-neutral Monte Carlo engine for spatial
cascading state transitions. Its labels and parameters come from JSON, allowing
the same A → B → C → D model to represent biosecurity, urban systems, financial
risk networks, utility grids, or other spatial scenarios.

The farm-risk map includes New Zealand boundaries. Set
`map_options.show_region_boundaries` in `config.json` to `true` or `false` to
show or hide regional council outlines without changing the simulation.

## Interactive fixed-farm cascade app

Run `shiny::runApp("spatial_cascade_app.R")` to explore one fixed-farm
scenario. Select one or more initiating farms, adjust transmission and timing
assumptions, run the scenario, and move the day slider to inspect the chain.
The app colours initiating farms, farms reached through the chain, and farms
not affected. This app supports scenario exploration; use
`spatial_cascade_monte_carlo.R` for probability estimates across many runs.

The app also supports a culling response policy. `Cull` stops outgoing
transmission from a covered farm after the response delay but counts that
farm's production capacity as lost. Policy coverage represents the proportion
of active farms that can be reached by the response programme. `None` provides
the no-policy comparison.

The Monte Carlo analysis also runs a matched no-control versus culling
comparison. Configure response coverage and delay under `cull_response` in
`config.json`. Covered farms are culled as soon as the
response delay expires. Its run-level and summary results, plus
the comparison chart, are written to `assets/` and interpreted in the Quarto
report.

North and South Island trucks are always kept separate.

## Main files

- `app.R` — Shiny application
- `truck_cluster_analysis.R` — geographic cluster analysis and SVG outputs
- `truck_data.R` — reproducible dummy truck population or custom-data loader
- `project_config.R` — source column mappings, brand colours, and speed settings
- `truck_analysis_report.qmd` — plain-language Quarto report
- `truck_analysis_report.html` — rendered self-contained report
- `data/trucks_clustered.rds` — stable clustered data loaded by the app
- `spatial_cascade_monte_carlo.R` — fixed-population spatial cascade engine
- `spatial_cascade_random_geography_legacy.R` — preserved random-layout version
- `create_fixed_farm_population.R` — rebuilds the example fixed farm table
- `data/farms.csv` — fixed farm population used by the cascade model
- `config.json` — external parameters and labels for the cascade engine
- `spatial_cascade_report.qmd` — short method and results report
- `assets/` — saved ggplot charts and reusable simulation results

## Required R packages

```r
install.packages(c(
  "shiny", "leaflet", "dplyr", "DT", "ggplot2", "cluster", "maps",
  "jsonlite"
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

Run the generic spatial cascade simulation:

```r
source("spatial_cascade_monte_carlo.R")
```

Edit `config.json` to change the scenario, population, spatial bounds, phase
durations, protection settings, production value, production-type categories,
and Monte Carlo run count. Production type is a reporting split only and does
not alter transmission.

By default, population size comes from `model_parameters.population_size`. To
use the number of rows in a farm dataset, set
`input_data.use_row_count_for_population_size` to `true` and save the CSV at the
configured `input_data.file` path. Only its row count is used; the Monte Carlo
engine still generates new random coordinates and attributes for every run.
Every run uses the same farm geography and attributes. Randomness comes from
the event origin and transmission outcomes. The script saves population,
production-type and farm-level risk results under `assets/`.

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
