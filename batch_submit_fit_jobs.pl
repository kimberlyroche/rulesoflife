use strict;
use warnings;

my @baboons = qw(ACA ALE CAI COB COO DAG DAS DUI DUN DUX EAG ECH ELV FAX GAB 
                 HOL HON LAD LAZ LEB LIW LOB LOG LOL LUI LYE MBE NAP NIK NOB 
                 ODE OFR ONY OPH ORI OTI OXY PEB SEB TAL THR VAI VAP VEI VEL 
                 VET VEX VIG VIN VOG VOT WAD WHE WIP WRI WYN);

my $filename = "job.slurm";
my $tax_level = "family";
my $MAP = "FALSE";
my $days_min_cor = 90;
my $diet_weight = 0.25;
my $var_scale_taxa = 1;
my $var_scale_samples = 1;
my $output_dir = "fam_days90_diet25_scale1";
my $use_adam = "FALSE";

for my $baboon (@baboons) {
  open(my $fh, '>', $filename);
  print $fh '#!/bin/bash'."\n";
  print $fh '#SBATCH -J '.$baboon."\n";
  print $fh '#SBATCH --mem=16GB'."\n";
  print $fh '#SBATCH --get-user-env'."\n";
  print $fh '#SBATCH --time=00:30:00'."\n";
  print $fh '#'."\n\n";

  print $fh 'cd /data/mukherjeelab/roche/rulesoflife'."\n\n";

  print $fh 'srun Rscript analysis_fit_model.R --sname='.$baboon.
                                             ' --tax_level='.$tax_level.
                                             ' --MAP='.$MAP.
                                             ' --output_dir='.$output_dir.
                                             ' --days_min_cor='.$days_min_cor.
                                             ' --diet_weight='.$diet_weight.
                                             ' --var_scale_taxa='.$var_scale_taxa.
                                             ' --var_scale_samples='.$var_scale_samples.
                                             ' --scramble_sample=FALSE'.
                                             ' --use_adam='.$use_adam."\n\n";

  close $fh;

  my $call_str = "sbatch $filename";
  print("Calling: ".$call_str."\n");
  `$call_str`;

  # the world's laziest delay
  my $lazy = 0;
  while($lazy < 1000000) { $lazy++; }
}
