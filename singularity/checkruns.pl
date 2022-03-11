#!/usr/bin/perl

my $targetbase = "/scratch/00451/tg457210/af/w";

$notes = "
usage: $0 ids

id is a number of the split job

for each id various sanity checks:
reports count of ids vs non pp txzs found



";


die $notes if !@ARGV;

$fmt = "%3s | %7s | %7s | %7s | %7s | %7s | %7s | %7s | %7s | %7s | %7s | %7s | %s\n";
    
print sprintf( $fmt
               ,"uid"
               ,"ids"
               ,"pps"
               ,"pdbs"
               ,"pdbtf"
               ,"mmcif"
               ,"pr"
               ,"csv"
               ,"cd"
               ,"txz"
               ,"zip"
               ,"mongo"
               ,"warnings"
    );

while ( @ARGV ) {
    my $id = shift;
    
    my $uid   = '0'x(3 - length($id)) . $id;

    my $base  = "$targetbase/w_$uid";
    my $ids   = "$base/ids";
    my $pdb   = "$base/af_pdbstage2";
    my $pdbtf = "$base/af_pdb_tfrev";
    my $cif   = "$base/af_mmcif";
    my $pr    = "$base/af_pr";
    my $csv   = "$base/af_csv";
    my $cd    = "$base/af_cd";
    my $txz   = "$base/af_txz";
    my $zip   = "$base/af_zip";
    my $mongo = "$base/af_mongocmds";

    die "$ids does not exist\n" if !-e $ids;

    my @ids   = `cat $ids`;
    my @pdb   = `cd $pdb && ls -1 | grep pdb`;
    my @pdbtf = `cd $pdbtf && ls -1`;
    my @cif   = `cd $cif && ls -1`;
    my @pr    = `cd $pr && ls -1`;
    my @cd    = `cd $cd && ls -1`;
    my @csv   = `cd $csv && ls -1`;
    my @txz   = `cd $txz && ls -1`;
    my @zip   = `cd $zip && ls -1`;
    my @mongo = `cd $mongo && ls -1`;

    @pps = grep /pp/, @pdb;
    
    my $warn = "";

    if ( scalar @ids + scalar @pps != scalar @pdb ) {
        $warn .= "pdb ";
    }

    if ( scalar @ids + scalar @pps != scalar @pdbtf ) {
        $warn .= "pdbtf ";
    }

    if ( scalar @ids + scalar @pps != scalar @cif ) {
        $warn .= "mmcif ";
    }

    if ( scalar @ids + scalar @pps != scalar @pr ) {
        $warn .= "pr ";
    }

    if ( scalar @ids + scalar @pps != scalar @cd ) {
        $warn .= "cd ";
    }

    if ( scalar @ids + scalar @pps != scalar @csv ) {
        $warn .= "csv ";
    }

    if ( scalar @ids + scalar @pps != scalar @txz ) {
        $warn .= "txz ";
    }

    if ( scalar @ids + scalar @pps != scalar @zip ) {
        $warn .= "zip ";
    }

    if ( scalar @ids + scalar @pps != scalar @mongo ) {
        $warn .= "mongo ";
    }

    $warn = "ok" if !$warn;

    print sprintf( $fmt
                   ,$uid
                   ,scalar @ids
                   ,scalar @pps
                   ,scalar @pdb
                   ,scalar @pdbtf
                   ,scalar @cif
                   ,scalar @pr
                   ,scalar @csv
                   ,scalar @cd
                   ,scalar @txz
                   ,scalar @zip
                   ,scalar @mongo
                   ,$warn
        );
    
    $sum{ids}   += scalar @ids;
    $sum{pps}   += scalar @pps;
    $sum{pdb}   += scalar @pdb;
    $sum{pdbtf} += scalar @pdbtf;
    $sum{cif}   += scalar @cif;
    $sum{pr}    += scalar @pr;
    $sum{csv}   += scalar @csv;
    $sum{cd}    += scalar @cd;
    $sum{txz}   += scalar @txz;
    $sum{mongo} += scalar @mongo;
    $sum{zip}   += scalar @zip;
}


my $warn = "";

if ( $sum{ids} + $sum{pps} != $sum{pdb} ) {
    $warn .= "pdb ";
}

if ( $sum{ids} + $sum{pps} != $sum{pdbtf} ) {
    $warn .= "pdbtf ";
}

if ( $sum{ids} + $sum{pps} != $sum{cif} ) {
    $warn .= "mmcif ";
}

if ( $sum{ids} + $sum{pps} != $sum{pr} ) {
    $warn .= "pr ";
}

if ( $sum{ids} + $sum{pps} != $sum{cd} ) {
    $warn .= "cd ";
}

if ( $sum{ids} + $sum{pps} != $sum{csv} ) {
    $warn .= "csv ";
}

if ( $sum{ids} + $sum{pps} != $sum{txz} ) {
    $warn .= "txz ";
}

if ( $sum{ids} + $sum{pps} != $sum{zip} ) {
    $warn .= "zip ";
}

if ( $sum{ids} + $sum{pps} != $sum{mongo} ) {
    $warn .= "mongo ";
}

$warn = "ok" if !$warn;

print sprintf( $fmt
               ,"tot"
               ,$sum{ids}
               ,$sum{pps}
               ,$sum{pdb}
               ,$sum{pdbtf}
               ,$sum{cif}
               ,$sum{pr}
               ,$sum{csv}
               ,$sum{cd}
               ,$sum{txz}
               ,$sum{zip}
               ,$sum{mongo}
               ,$warn
    );
