#!/usr/bin/perl

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

check chains files, determine variants and create variants for processing

";

die $notes if !@ARGV;

p_config_init();
require "$ppdir/pdbutil.pm";

$errors = "";

use Math::Combinatorics;

sub pdb_rechain_summary {
    my $pdb = shift || error_exit( "$0: pdb_rechain_summary requires an argument" );

    ## update chain id at resseq breaks and return string of ranges present in pdb

    my $nextresseq = 0;
    my $result = "";
    my $chain  = "A";
    my $resseq;
    
    for my $l ( @$pdb ) {
        next if $l !~ /^(TER|ATOM|ENDMDL)/;
        if ( $l =~ /^TER/ ) {
            $l = '';
            next;
        }
        if ( $l =~ /^ENDMDL/ ) {
            $l = "TER\n$l";
            next;
        }
        $resseq = mytrim( substr( $l, 22, 4 ) );
        if ( !$nextresseq ) {
            $result     .= "$chain:$resseq";
            $nextresseq  = $resseq + 1;
            next;
        }
        my $addter;
        if ( $resseq > $nextresseq ) {
            $chain   = ++$chain;
            $result .= "-" . ( $nextresseq - 1 ) . "; $chain:$resseq";
            $addter  = 1;
        }
        $nextresseq = $resseq + 1;
        substr( $l, 21, 1 ) = $chain;
        $l = "TER\n$l" if $addter;
    }
    "$result-$resseq";
}

while ( $fid = shift ) {
    my $id = uniprot_id( $fid );
    log_append( "$0: processing ascension $id\n" );

    my @chains = get_chains( $id );
    my @pdbs   = af_pdbs   ( $id );

    my @types;
    my @begins;
    my @ends;
    
    my $fn       = "$$p_config{pdbnotes}/$id.notes";
    unlink $fn if -e $fn;

    my @pdbcomments;

    for my $l ( @chains ) {
        if ( $l =~ /^\s*#/ ) {
            push @pdbcomments, $l;
            next;
        }
        my @l = split /\s+/, $l;
        if ( @l != 3 ) {
            my $error = "$0: ERROR - $id $f chains line $l does not contain 3 tokens\n";
            log_append( $error );
            $errors .= $error;
            next;
        }
        push @types,  $l[0];
        push @begins, $l[1];
        push @ends,   $l[2];
    }

    log_append( "$0: pdbs found", "-" );
    log_append( join "\n", @pdbs );
    log_append( "\n" );

    log_append( "$0: chain summary", "-" );
    log_append( join "\n", @pdbcomments );
    log_append( "\n" ) if @pdbcomments;
    for ( my $i = 0; $i < @types; ++$i ) {
        log_append( "$types[$i] $begins[$i] $ends[$i]\n" );
    }

    ## setup additonal flags

    my $peptide_only = grep /only peptides/, @pdbcomments;

    ## create cut sets
    my $cuts = \{};

    $$cuts->{always} = {};
    $$cuts->{always}->{type}  = [];
    $$cuts->{always}->{begin} = [];
    $$cuts->{always}->{end}   = [];

    ### ignore end_assumed for now
    grep s/~end_assumed//, @types;

    ### always cuts

    my %alwaysinitcuts = (
        "init_met" => "initiator methionine"
        ,"signal"  => "signal peptide"
        ,"transit" => "transit"
        );

    while( @types && $alwaysinitcuts{ $types[0] } ) {
        push @{$$cuts->{always}->{type}} , shift @types;
        push @{$$cuts->{always}->{begin}}, shift @begins;
        push @{$$cuts->{always}->{end}}  , shift @ends;
    }

    if ( $debug > 1 ) {
        log_append( debug_json( "cuts after always", $cuts ) );
        log_append( "$0: chain summary after cuts always", "-" );
        log_append( join "\n", @pdbcomments );
        log_append( "\n" ) if @pdbcomments;

        for ( my $i = 0; $i < @types; ++$i ) {
            log_append( "$0: $types[$i] $begins[$i] $ends[$i]\n" );
        }
    }

    if ( !$peptide_only ) {
        for ( my $i = 0; $i < @types; ++$i ) {
            if ( $types[$i] eq 'propep' ) {
                if ( !exists $$cuts->{permute} ) {
                    $$cuts->{permute} = {};
                    $$cuts->{permute}->{type}  = [];
                    $$cuts->{permute}->{begin} = [];
                    $$cuts->{permute}->{end}   = [];
                }
                push @{$$cuts->{permute}->{type}} , $types [$i];
                push @{$$cuts->{permute}->{begin}}, $begins[$i];
                push @{$$cuts->{permute}->{end}}  , $ends  [$i];
            }
        }
    }

    if ( $debug > 1 ) {
        log_append( debug_json( "cuts after permute", $cuts ) );

        log_append( sprintf( "size of permutes list %d\n", scalar @{$$cuts->{permute}->{type}} ) );
    }

    ## find variants, make pdbs

    ### always cuts for base pdb
    my %always_cuts;
    my @always_remarks;

    for ( my $i = 0; $i < @{$$cuts->{always}->{type}}; ++$i ) {
        for ( my $j = $$cuts->{always}->{begin}[$i]; $j <= $$cuts->{always}->{end}[$i]; ++$j ) {
            $always_cuts{$j}++;
        }
        if ( $$cuts->{always}->{type}[$i] eq 'init_met' ) {
            push @always_remarks, "REMARK   2 UNIPROT $alwaysinitcuts{$$cuts->{always}->{type}[$i]} removed\n";
        } else {
            push @always_remarks, "REMARK   2 UNIPROT $alwaysinitcuts{$$cuts->{always}->{type}[$i]} seq. $$cuts->{always}->{begin}[$i]-$$cuts->{always}->{end}[$i] removed\n";
        }
    }

    my %permutetags = (
        "propep" => "propeptide"
    );

    log_append( "build pdbs", "-" );

    for my $fpdb ( @pdbs ) {
        my $pdb_ver       = pdb_ver  ( $fpdb );
        my $pdb_frame     = pdb_frame( $fpdb );

        my @lpdb          = get_pdb_lines( $fpdb );
        my $pdb_res_count = get_pdb_residue_count( $fpdb );

        if ( $debug > 1 ) {
            log_append( "$fpdb - pdb_res count $pdb_res_count\n" );
        }
        
        # file with always cuts

        {
            my $fo = "$id-${pdb_frame}-$pdb_ver.pdb";
            if ( $debug > 1 ) {
                log_append( "$fpdb ver '$pdb_ver' frame '$pdb_frame'\nto create as '$fo'", "-" );
                log_append( ( join '', @always_remarks ) . line() );
            }

            # create pdb

            my @newpdb = @always_remarks;
            for my $l ( @lpdb ) {
                if ( $l !~ /^ATOM/ ) {
                    push @newpdb, "$l\n";
                    next;
                }
                my $resseq = mytrim( substr( $l, 22, 4 ) );
                next if $always_cuts{ $resseq };
                push @newpdb, "$l\n";
                next;
            }
            my $range = pdb_rechain_summary( \@newpdb );
            unshift @newpdb, "REMARK   2 RESIDUE seq. $range\n";
            write_file( "$$p_config{pdbstage1}/$fo", join '', @newpdb );
            log_append( "$fo\n" );
        }

        # create all permutations as needed

        if ( exists $$cuts->{permute} ) {
            my @words = (0 .. (@{$$cuts->{permute}->{type}} - 1 ));
            for my $count (1 .. @words) {
                my $comb = Math::Combinatorics->new(count => $count,
                                                    data  => [@words]);
                while (my @combo = $comb->next_combination) {
                    @combo = sort { $a <=> $b } @combo;

                    my $pps = "pp";
                    for my $v ( @combo ) {
                        $pps .= ($v+1) . '_';
                    }
                    $pps =~ s/_$/-/;

                    my $fo = "$id-${pdb_frame}-$pps$pdb_ver.pdb";
                    if ( $debug > 1 ) {
                        log_append( "$fpdb ver '$pdb_ver' frame '$pdb_frame'\nto create as '$fo'", "-" );
                    }

                    ## build up cut range
                    my %permute_cuts;
                    my @permute_remarks;
                    
                    for my $i ( @combo ) {
                        for ( my $j = $$cuts->{permute}->{begin}[$i]; $j <= $$cuts->{permute}->{end}[$i]; ++$j ) {
                            $permute_cuts{$j}++;
                        }
                        push @permute_remarks, "REMARK   2 UNIPROT $permutetags{$$cuts->{permute}->{type}[$i]} seq. $$cuts->{permute}->{begin}[$i]-$$cuts->{permute}->{end}[$i] removed\n";
                    }
                    if ( $debug > 1 ) {
                        log_append( ( join '', @always_remarks ) . ( join '', @permute_remarks ) . line() );
                    }

                    ## create pdb
                    my @newpdb = @always_remarks;
                    push @newpdb, @permute_remarks;
                    for my $l ( @lpdb ) {
                        if ( $l !~ /^ATOM/ ) {
                            push @newpdb, "$l\n";
                            next;
                        }
                        my $resseq = mytrim( substr( $l, 22, 4 ) );
                        next if $always_cuts{ $resseq } || $permute_cuts{ $resseq };
                        push @newpdb, "$l\n";
                        next;
                    }
                    my $range = pdb_rechain_summary( \@newpdb );
                    unshift @newpdb, "REMARK   2 RESIDUE seq. $range\n";
                    write_file( "$$p_config{pdbstage1}/$fo", join '', @newpdb );
                    log_append( "$fo\n" );
                }
            }
        }
    }

    write_file( $fn, log_final() );
}

log_flush();
if ( $errors ) {
    print STDERR $errors;
}

exit length( $errors );    
    
