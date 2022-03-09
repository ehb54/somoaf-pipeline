#!/usr/bin/perl

## 

$notes = "
usage: $0 basedir targetdir count file

split list into multiple sets for smaller file system

basedir is location of original work directory
targetdir is location of new work directories


";

use POSIX;
use File::Temp qw(tempfile);

$bd  = shift || die $notes;
$td  = shift || die $notes;
$cnt = shift || die $notes;
$f   = shift || die $notes;


die "count must be > 1\n"       if $cnt <= 1;
die "$bd is not a directory\n"  if !-d $bd;
die "$td is not a directory\n"  if !-d $td;
die "$f not found\n"            if !-e $f;
die "$f not readable\n"         if !-r $f;

open IN, $f;
@l = <IN>;
close IN;
grep chomp, @l;

## split remaining args into 

my $jpc = ceil( @l / $cnt );

$tmpdir = "/scratch/00451/tg457210/af/work/af_tmp";

@cds = (
    "af_blast"
    ,"af_cd"
    ,"af_chains"
    ,"af_cif"
    ,"af_cperrors"
    ,"af_cpnotes"
    ,"af_csv"
    ,"af_features"
    ,"af_ftp"
    ,"af_mmcif"
    ,"af_mongocmds"
    ,"af_pdb"
    ,"af_pdbnotes"
    ,"af_pdbstage1"
    ,"af_pdbstage2"
    ,"af_pdb_tfrev"
    ,"af_pr"
    ,"af_tmp"
    ,"af_txz"
    ,"af_zip"
    );

@cpd = (
    "af_blast"
    ,"af_pdb"
    ,"af_features"
    );

my $i = 0;

while ( @l ) {
    my ( $fh, $ft ) = tempfile( "$tmpdir/xargs.XXXXXX", UNLINK => 0 );
    my @args = splice @l, 0, $jpc;
    print $fh join "\n", @args;
    close $fh;

    my $iu = '0'x(3-length($i)) . $i;
    
    my $prefix = "$td/w_$iu";

    my $cmd =
        "mkdir $td/w_$iu && cd $td/w_$iu && \\\n"
        . "mkdir " . ( join ' ',@cds ) . " && \\\n";

    for my $cpd ( @cpd ) {
        $cmd .= "cd $bd/$cpd && ls | grep -f $ft | xargs -P4 cp -t $td/w_$iu/$cpd/ && \\\n";
    }
    $cmd .= "echo $w_$iu done\n";

    print $cmd;
    print `$cmd`;
    ++$i;
}

       
        
    
