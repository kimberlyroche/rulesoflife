use strict;
use warnings;
use POSIX;
use List::Util qw(min);

my $no_pairs = 8911;
my $segment_size = 100;

my $offset = 0;

my $job_no = ceil($no_pairs / $segment_size);

my $filename = "job.slurm";
my $start = -1;
my $end = -1;

for my $jn (1..$job_no) {
  $start = (($jn - 1) * $segment_size) + 1;
  $end = $start + $segment_size - 1;
  $end = min($end, $no_pairs);

  $start = $start + $offset;
  $end = $end + $offset;

  print("START: ".$start.", END: ".$end."\n");

  open(my $fh, '>', $filename);
  print $fh '#!/bin/bash'."\n";
  print $fh '#SBATCH -J pairs'.$start.'-'.$end."\n";
  print $fh '#SBATCH --mem=16GB'."\n";
  print $fh '#SBATCH --get-user-env'."\n";
  print $fh '#SBATCH --time=1:30:00'."\n";
  print $fh '#'."\n\n";

  print $fh 'cd /data/mukherjeelab/roche/rulesoflife'."\n\n";

  print $fh 'srun Rscript analysis_synchrony_distribution_generation.R --start='.$start.' --end='.$end."\n\n";

  close $fh;

  my $call_str = "sbatch $filename";
  print("Calling: ".$call_str."\n");
  `$call_str`;

  # the world's laziest delay
  my $lazy = 0;
  while($lazy < 1000000) { $lazy++; }
}
