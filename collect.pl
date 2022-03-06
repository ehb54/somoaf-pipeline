#!/usr/bin/perl

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

verify all downloaded data for processing & download if needed

";

die $notes if !@ARGV;

p_config_init();

$errors = "";

my $http = HTTP::Tiny->new();

while ( $fid = shift ) {
    my $id = uniprot_id( $fid );
    log_append( "$0: processing ascension $id\n" );

    ## verify or get fasta

    {
        my $f = "$$p_config{fasta}/$id.fasta";
        if ( !-e $f || -z $f ) {
            ## get fasta
            my $url = "$$p_config{fasta_url}/$id.fasta";
            my $json = $http->get( $url );
            if ( !$json->{success} ) {
                my $error = "could not retrive $f from $url\n";
                log_append( "$0: ERROR $error" );
                $errors .= $error;
                next;
            }
            my $fasta = $json->{content};
            open OUT, ">$f";
            print OUT $fasta;
            close OUT;
            log_append( "$0: ok : $f retrieved\n" );
        } else {
            log_append( "$0: ok : $f exists\n" );
        }
    }

    ## do we have the pdb ?

    {
        my @f = af_pdbs( $id );
        grep chomp, @f;
        if ( !@f ) {
            my $error = "no pdb found in $$p_config{afpdb}\n";
            log_append( "$0: ERROR $error" );
            $errors .= $error;
            next;
        }
        log_append( "$0: ok : " . scalar @f . " pdb(s) found for $id\n" );
    }

    ## verify or get uniprot features
    {
        my $f = "$$p_config{unifeat}/$id.uf.txt";
        if ( !-e $f ) {
            my $url = "$$p_config{unifeat_url}$id";
            my $json = $http->get($url, {
                headers => { 'Accept' => 'application/json' }
                                  });

            if ( !$json->{success} ) {
                my $error = "could not retrive $f from $url\n";
                log_append( "$0: ERROR $error" );
                $errors .= $error;
                next;
            }

            my $res = decode_json( $json->{content} );

            my $out = "[]";
            if ( exists $$res[0]->{features} ) {
                my $feat = $$res[0]->{features};
                $out = encode_json( $feat );
            }

            open OUT, ">$f";
            print OUT $out;
            close OUT;
            log_append( "$0: ok : $f retrieved\n" );
        } else {
            log_append( "$0: ok : $f exists\n" );
        }
    }
}

log_flush();
if ( $errors ) {
    print STDERR $errors;
}

exit length( $errors );    
    
