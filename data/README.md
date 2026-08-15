# Clustered truck data

Run `source("truck_cluster_analysis.R")` from the project root to create:

- `trucks_clustered.rds` — the stable data object loaded by the Shiny app;
- `trucks_clustered.csv` — the same data in a human-readable format.

Rerunning the cluster analysis deliberately refreshes the cluster assignments.

To use another dataset, save it as `data/trucks_input.csv`. If its column names
differ, edit `TRUCK_COLUMNS` in `project_config.R`, then rerun the analysis.

Brand values are open-ended. Add brands directly to the input data. Optional
preferred colours and age/category speed multipliers are also configured in
`project_config.R`.
