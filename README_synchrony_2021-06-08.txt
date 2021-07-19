Generating the synchrony vs. universality raw materials

1) Run 'analysis_interpolate_series.R'
     This generates the 'host_mean_predictions_adjusted.rds' file.
2) Run batch_between_within_jobs.pl
     This runs 'analysis_synchrony_distribution_generation.R', which a) uses 
     'host_mean_predictions_adjusted.rds' and b) generates the 
     'within_between_distros-XXX' output files.
3) Run 'analysis_stitch_within_between.R'
     This stitches those files into one 'within_between_distros.rds' file.
4) Run 'figures_synchrony_all_pairs.R'
     This plots using 'within_between_distros.rds'.

