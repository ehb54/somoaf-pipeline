#!/usr/bin/perl

my $targetbase = "/scratch/00451/tg457210/af/w";

$slurm = <<__EOF;
#!/bin/bash

#SBATCH -J af_wxxx    # Job name
#SBATCH -o slurm_%j.out
#SBATCH -e slurm_%j.err
#SBATCH -p normal          # Queue (partition) name
#SBATCH -N 1               # Total # of nodes (must be 1 for serial)
#SBATCH -n 1               # Total # of mpi tasks (should be 1 for serial)
#SBATCH -t 04:00:00        # Run time (hh:mm:ss)
#SBATCH --mail-type=all    # Send email at begin and end of job
#SBATCH --mail-user=brookes@uthscsa.edu

# Any other commands must follow all #SBATCH directives...
cd $targetbase/w_xxx/afpp
module load tacc-singularity
pwd
date
singularity/start_mongo.sh
echo "mongo started, now run"    
singularity run ~/somoafpipe.sif bash -c "export HOME=$targetbase/w_xxx/afpp && singularity/xargs.pl 136 ./run.pl $targetbase/w_xxx/ids >xxx.out 2>xxx.err"

__EOF
;

$notes = "
usage: $0 ids

id is a number of the split job

builds slurm scripts for ids
targetbase for jobs is $targetbase

";


die $notes if !@ARGV;

while ( @ARGV ) {
    my $id = shift;
    my $uid = '0'x(3 - length($id)) . $id;
    my $uslurm = $slurm;

    $uslurm =~ s/xxx/$uid/mg;

    my $f = "somoaf_${uid}.slurm";
    open my $fh, ">$f";
    print $fh $uslurm;
    close $fh;
    print "$f\n";
}
  
