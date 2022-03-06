#!/usr/bin/perl

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

helper program to return unretrieved ids for 'collect.pl'

returns unprocessed ids

";

die $notes if !@ARGV;

p_config_init();

while ( $fid = shift ) {
    my $id = uniprot_id( $fid );
    
    my $fasta = "$$p_config{fasta}/$id.fasta";
    my $ffeat = "$$p_config{unifeat}/${id}.uf.txt";
    next if -e $fasta && -e $ffeat;

    print "$id\n";
}
