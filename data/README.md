# Data and provenance

The numerical inputs are kept separate from the analysis scripts and generated
results. The files use original length units unless explicitly indicated.

## TSPLIB

- `tsplib_optima.csv`: the 77 instance names, exact tour lengths, number of
  nodes, bounding-box area, Clark–Evans statistic and random-tour SD extracted
  from the authors' supplied `default_tsp_model.rds` (training data).
- `tsplib_coords.csv`: 138,085 node coordinates taken from the TSPLIB files
  in [mastqe/tsplib](https://github.com/mastqe/tsplib), commit
  `775541365dc3f7f3cd6faa8d0e206036d2c45e0e`, in their original node order.
  TSPLIB: Reinelt (1991), ORSA Journal on Computing, 3(4), 376–384.

The 77 lengths agree with the published solutions in that mirror. The number
of nodes and bounding-box area also agree with the authors' archived data for
every instance. Five explicit-distance instances and the duplicate linhp318
are not included, as specified in the manuscript.

The TSPLIB analysis uses the archived random-tour SD by default and calculates
the other descriptors from the coordinates. Set `recompute_sd = TRUE` in
`analysis/05_tsplib.R` to repeat the 1,000 random tours; the script compares the
recomputed values with the archived data.

## Main simulation

- `simulated_instances.csv`: 150 simulated instances, extracted from the
  three archived R data files in `source/`. The original heuristic lengths
  are recorded as `opt_heuristic` and are not used as response values.
- `simulation_optima.csv`: optimal lengths matched by instance name from
  `source/results_concorde.xlsx` and divided by 1,000,000.
- `source/`: the three original R data files from
  [albrizre/Spatial_TSP_Regression](https://github.com/albrizre/Spatial_TSP_Regression),
  commit `f524203cda3ca32bc8cf52150eb7b613a361cba2`, retained for provenance,
  and the updated `results_concorde.xlsx` supplied by the authors.

The updated spreadsheet supplies exact lengths for all 148 retained
instances. Only `clust_ins_3` and `clust_ins_7` have missing lengths; both
are excluded by the manuscript's SD filter (SD < 1). No manual import is
needed: the converted CSV is included. The analysis stops if a retained
instance has no exact length.

## Sensitivity experiment

The authors' original `data_sensitivity_clust.rda` is included as
`cluster_sensitivity.rda` (object `sens`). It contains all 360 instances,
with 10 replicates per combination of parent intensity, cluster scale and
sample size. All lengths are positive and marked as Concorde solutions.
The analysis uses this archive by default.

The optional generation scripts produce a new experiment in `generated/`.
Do not combine reference lengths from an earlier experiment with regenerated
instances merely because their names are the same.
