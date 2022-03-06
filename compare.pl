#!/usr/bin/perl

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

compares v1 & v2 files

";

die $notes if !@ARGV;

p_config_init();

$errors = "";
use File::Temp qw(tempfile);

sub diffrpt {
    my ( $fh0, $ft0 ) = tempfile( "$$p_config{tempdir}/$id-v1.XXXXXX", UNLINK => 1 );
    my ( $fh1, $ft1 ) = tempfile( "$$p_config{tempdir}/$id-v2.XXXXXX", UNLINK => 1 );

    print $fh0 join "\n", @pdb0;
    print $fh1 join "\n", @pdb1;

    print $fh0 "\n";
    print $fh1 "\n";
    
    close $fh0;
    close $fh1;
    
    my $diffs = run_cmd( "diff $ft0 $ft1", true );

    my @diffs = split /\n/, $diffs;
    @diffs             = grep /^. ATOM/, @diffs;
    $diffs_count       = ( scalar @diffs ) / 2;
    @diffs             = grep !/NH. ARG /, @diffs;
    $diffs_not_ng_only = scalar @diffs;

    "$id\n$diffs\n";
}


while ( $fid = shift ) {
    $id = uniprot_id( $fid );
    # log_append( "$0: processing ascension $id\n" );

    my @pdbs   = af_pdbs   ( $id );
    @pdbs = grep /-F1-/, @pdbs;

    if ( @pdbs != 2 ) {
        my $error = "$0: ERROR - $id only one version present\n";
        $errors .= $error;
        next;
    }

    @pdb0 = get_pdb_lines( $pdbs[0] );
    @pdb0    = grep /^ATOM/, @pdb0;

    @pdb1 = get_pdb_lines( $pdbs[1] );
    @pdb1    = grep /^ATOM/, @pdb1;

    if ( @pdb0 != @pdb1 ) {
        print STDERR diffrpt();
        log_append( sprintf( "$0: $id differ ($diffs_count) %s\n", $diffs_not_ng_only ? "many types" : "NG only" ) );
        log_flush();
        next;
    }
        
    my $match = 1;
    for ( my $i = 0; $i < @pdb0; ++$i ) {
        if ( $pdb0[$i] cmp $pdb1[$i] ) {
            $match = 0;
            last;
        }            
    }
    if ( !$match ) {
        print STDERR diffrpt();
        log_append( sprintf( "$0: $id differ ($diffs_count) %s\n", $diffs_not_ng_only ? "many types" : "NG only" ) );
        log_flush();
        next;
    }

    log_append( "$id match\n" );
    log_flush();
}

log_flush();

if ( $errors ) {
    print STDERR $errors;
}

exit length( $errors );    

        
        
