#!/usr/bin/perl

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

runs chain of processes described in $configfile

";

die $notes if !@ARGV;

p_config_init();

sub stage_process {
    my $s  = shift;
    my $id = shift || error_exit( "$0: stage_process() requires 2 arguments" );
    my @reqs = ( "name", "cmd" );
    for my $k ( @reqs ) {
        error_exit( "$0: stage_process(): $configfile stages entry missing $k" ) if !$s->{$k};
    }
        
    log_append( "stage $s->{name}", "-" );
    if ( exists $s->{active} && !$s->{active} ) {
        log_append( "stage not active, skipped\n" );
        log_flush();
        return;
    }

    my $cmd = $s->{cmd};
    $cmd =~ s/__base__/$$p_config{base}/g;
    $cmd =~ s/__id__/$id/g;

    log_append( $cmd, "+" );
    $res = run_cmd( $cmd, true );
    log_append( "$res\n" );

    log_flush();

    run_cmd_last_error();
}

my $stages = $$p_config{stages};

while ( $fid = shift ) {
    my $id = uniprot_id( $fid );
    log_append( "processing ascension $id", "=" );
    for ( my $i = 0; $i < scalar @$stages; ++$i ) {
        my $s = $$stages[$i];
        if ( stage_process( $s, $id ) ) {
            log_append( "errors returned for $id, skipping further processing", "*" );
            last;
        }
    }
}
        
log_flush();

# print log_final();

