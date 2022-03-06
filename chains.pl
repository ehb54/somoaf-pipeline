#!/usr/bin/perl

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

builds chains files if empty or non-existant

";

die $notes if !@ARGV;

p_config_init();

$errors = "";

while ( $fid = shift ) {
    my $id = uniprot_id( $fid );
    log_append( "$0: processing ascension $id\n" );

    my $f        = "$$p_config{chains}/$id.chains";
    my $fn       = "$$p_config{cpnotes}/$id.notes";
    my $fe       = "$$p_config{cperrs}/$id.error";

    my $oktoskip = 1;

    if ( -e $f && !-z $f ) {
        log_append( "$0: ok : $f exists\n" );
    } else {
        $oktoskip = 0;
    }
    if ( -e $fn ) {
        log_append( "$0: ok : $fn exists\n" );
    } else {
        $oktoskip = 0;
    }
    next if $oktoskip;

    unlink $fe if -e $fe;

    my $cpnotes;
    my @outnotes;

    ## fasta, pdb & features should exist

    my $fasta_res_count   = get_fasta_residue_count( $id );
    my $pdb_res_count     = get_residue_count      ( $id );

    log_append( "$0: fasta residue count $fasta_res_count\n" );
    log_append( "$0: pdb residue count $pdb_res_count\n" );

    if ( $pdb_res_count eq 'multiframe' ) {
        my $error = "$0: notice - multiframe pdb ignored\n";
        log_append( $error );
        $errors .= $error;
        next;
    }
    
    if ( $fasta_res_count != $pdb_res_count ) {
        my $error = "$0: warning - $id fasta ($fasta_res_count) & pdb residue ($pdb_res_count) count mismatch\n";
        write_file( $fe, $error );
        log_append( $error );
        # write_file( $fn, log_final() );
        my $warnornotice = abs( $fasta_res_count - $pdb_res_count ) > 2 ? "warning" : "notice";
        push @outnotes, "$warnornotice - alphafold pdb residue count ($pdb_res_count) mismatches fasta residue count ($fasta_res_count)";
        # $errors .= $error;
        # next;
    }        

    log_append( "$0: ok - fasta & pdb residues match $id\n" );

    ## process features to make chain info

    my $featstr = read_features( $id );

    if ( !$featstr ) {
        my $error = "$0: ERROR - empty features for $id\n";
        write_file( $fe, $error );
        log_append( $error );
        write_file( $fn, log_final() );
        $errors .= $error;
        next;
    }        
        
    eval {
        $features = decode_json( $featstr );
        1;
    } or do {
        my $error = "$0: ERROR - json decoding features for $id $@\n";
        write_file( $fe, $error );
        log_append( $error );
        write_file( $fn, log_final() );
        $errors .= $error;
        next;
    };
        
    log_append( "$0: ok - found features for $id\n" );

    ## ok now we can features
    my @types;
    my @begins;
    my @ends;
    
    my %counts;

    for my $feature ( @$features ) {
        log_append( "$feature->{category}\n" );
        next if $feature->{category} ne 'MOLECULE_PROCESSING';
        print Dumper( $feature ) if $debug > 1;

        my $type  = exists $feature->{type}  ? lc( $feature->{type} )  : 'unknown';
        my $begin = exists $feature->{begin} ? $feature->{begin}       : 'unknown';
        my $end   = exists $feature->{end}   ? $feature->{end}         : 'unknown';
        
        $counts{ $type }++;
        
        $type = "chain" if $type eq "peptide";
        $begin =~ s/<//;
        $begin =~ s/~(\d+)/$1/;
        $end   =~ s/>//;
        $end   =~ s/~(\d+)/$1/;

        if ( $type eq 'init_met' &&
             $begin > 1 ) {
            my $warning = "$0: notice $id - internal initial methionine ignored\n";
            log_append( $warning );
            push @outnotes, "notice - internal initial methionine at residue $begin ignored";
            next;
        }

        if ( $type eq 'init_met' &&
             $begin == 1 &&
             @types ) {
            @types  = ();
            @begins = ();
            @ends   = ();
            log_append( "$0: notice $id - reinitialized molecule processing at init_met\n" );
        }

        push @types , $type;
        push @begins, $begin;
        push @ends,   $end;
        
        last if $end == $fasta_res_count;

        if ( $end > $fasta_res_count ) {
            my $error = "$0: ERROR $id - ERROR end past fasta_res_count\n";
            log_append( $error );
            $errors .= $error;
            last;
        }
    }
    
    if ( !@types ) {
        log_append( "$0: $id no molecule processing - as complete chain\n" );
        my $outnotes;
        if ( @outnotes ) {
            grep s/^/# /, @outnotes;
            $outnotes = join "\n", @outnotes;
            $outnotes .= "\n";
        }
        write_file( $f, "${outnotes}chain 1 $fasta_res_count\n" );
        log_append( "$0: wrote $f\n" );
        log_flush();
        write_file( $fn, log_final() );
        next;
    } else {
        log_append( line() . "initial molecule processing chains\n" . line() );
        for ( my $i = 0; $i < @types; ++$i ) {
            log_append( "$types[$i] $begins[$i] $ends[$i]\n" );
        }
    }

    log_append( line() . "checking for sequential chains\n" . line() );

    ## special handling

    ### remove all propeptides if all peptides
    my $peptides_only = 0;

    if ( $counts{'peptide'} &&
         $counts{'propep'} &&
         !$counts{'chain'} ) {
        
        push @outnotes, "notice - only peptides present";
        $peptides_only = 1;
        
        ## logic for checking prior or next chain commented

        my $removed_propeptide = 0;

        my @ntypes;
        my @nbegins;
        my @nends;

        ## needed this when our 'for' started at index 1
        # my @ntypes  = @types [0];
        # my @nbegins = @begins[0];
        # my @nends   = @ends  [0];
        for ( my $i = 0; $i < @types; ++$i ) {
            ## previous types on stack

            # my $ptype  = $ntypes [-1];
            # my $pbegin = $nbegins[-1];
            # my $pend   = $nends  [-1];

            ## this type

            my $type  = $types [$i];
            my $begin = $begins[$i];
            my $end   = $ends  [$i];

            ## next type

            # my $ntype  = $types [$i + 1];
            # my $nbegin = $begins[$i + 1];
            # my $nend   = $ends  [$i + 1];
            
            if (
                $type eq 'propep' &&
                $begin ne '~' &&
                $end ne '~'
                ## alternate logic if before or after chain
                # &&
                # ( $ptype eq 'chain' ||
                # $ntype eq 'chain')
                ) {
                $removed_propeptide = 1;
                next;
            }

            push @ntypes,  $type;
            push @nbegins, $begin;
            push @nends,   $end;
        }
        if ( $removed_propeptide ) {
            ## needed this when our for stopped at @types - 1
            # push @ntypes,  $types [-1];
            # push @nbegins, $begins[-1];
            # push @nends,   $ends  [-1];
            
            @types  = @ntypes;
            @begins = @nbegins;
            @ends   = @nends;

            log_append( "special handling - $id removed propeptide from egg basket\n" );
        }
    }

    ## try to join sequential chain types
    {
        my @ntypes  = @types [0];
        my @nbegins = @begins[0];
        my @nends   = @ends  [0];
        for ( my $i = 1; $i < @types; ++$i ) {
            ## previous types on stack

            my $ptype  = $ntypes [-1];
            my $pbegin = $nbegins[-1];
            my $pend   = $nends  [-1];

            ## this type

            my $type  = $types [$i];
            my $begin = $begins[$i];
            my $end   = $ends  [$i];
            if (
                ( $ptype eq 'chain' &&
                  $type  eq 'chain' ) ||
                ( $ptype eq 'transit' &&
                  $type  eq 'transit' )

                ) {
                # log_append( "sequential chain processing\n" );
                if ( $end <= $pend &&
                     $begin >= $pbegin &&
                     $begin ne '~' &&
                     $end ne '~' &&
                     $pbegin ne '~' &&
                     $pend ne '~'
                    ) {
                    log_append( "sequential chain processing - skipped subseq sequential chain\n" );
                    next;
                }
                if ( $pend eq '~' && $begin eq '~' ) {
                    log_append( "sequential chain processing - joining tilde chains\n" );
                    $nends[ -1 ] = $end;
                    next;
                }
                if ( $pend + 1 == $begin &&
                     $pend ne '~' ) {
                    log_append( "sequential chain processing - joining seq chains\n" );
                    $nends[ -1 ] = $end;
                    next;
                }

                if ( $end ne '~' ) {
                    log_append( "sequential chain processing - joining seq $type with gap or overlap\n" );
                    $nends[ -1 ] = $end;
                    next;
                }

                if ( $pbegin <= $begin &&
                     $end eq '~' ) {
                    log_append( "sequential chain processing - skipped subseq sequential chain with tilde end\n" );
                    next;
                }
            }

            push @ntypes,  $type;
            push @nbegins, $begin;
            push @nends,   $end;
        }
        @types  = @ntypes;
        @begins = @nbegins;
        @ends   = @nends;
    }

    ## special handling

    ### no chains but signal, propep or other
    ### fill with dummy chain at beginning or end
    if ( !grep /chain/, @types ) {
        if ( $begins[0] == 1 && $ends[-1] == $fasta_res_count ) {
            ## do nothing, caught in warn
        } else {
            if ( $begins[0] == 1 ) {
                ## add chain at end
                push @types,  "chain";
                push @begins, $ends[-1] + 1;
                push @ends,   $fasta_res_count;
                log_append( "special handling - $id filled with chain at end\n" );
            } else {
                ## add chain at beginning
                unshift( @types, "chain" );
                unshift( @begins, 1 );
                unshift( @ends,   $begins[1] -1 );
                log_append( "special handling - $id filled with chain at beginning\n" );
            }
        }
    }

    ### gaps between chain and non-chain -> extend chain
    
    for ( my $i = 0; $i < @types; ++$i ) {
        my $type  = $types [$i];
        my $begin = $begins[$i];
        my $end   = $ends  [$i];

        my $extended_chain = 0;

        if ( $type eq 'chain' ) {
            if ( $i &&
                 $begin > $ends[ $i - 1 ] + 1 &&
                 $begin ne '~' &&
                 $ends[ $i - 1 ] ne '~' ) {
                $begins[ $i ] = $ends[ $i - 1 ] + 1;
                $extended_chain = 1;
            }
            if ( $i < @types - 1 &&
                 $end < $begins[ $i + 1 ] - 1 &&
                 $end ne '~' &&
                 $begins[ $i + 1 ] ne '~' ) {
                $ends[ $i ] = $begins[ $i + 1 ] - 1;
                $extended_chain = 1;
            }
        }            

        if ( $extended_chain ) {
            log_append( "special handling - $id extended chains to fill gaps\n" );
        }
    }

    ### overlap of non-chain to chain -> cut non-chain

    for ( my $i = 0; $i < @types; ++$i ) {
        my $type  = $types [$i];
        my $begin = $begins[$i];
        my $end   = $ends  [$i];

        my $cut_chain = 0;

        if ( $type eq 'chain' ) {
            if ( $i &&
                 $begin < $ends[ $i - 1 ] + 1 &&
                 $begin ne '~' &&
                 $ends[ $i - 1 ] ne '~'
                ) {
                $ends[ $i - 1 ] = $begin - 1;
                $cut_chain = 1;
            }
            if ( $i < @types - 1 &&
                 $end > $begins[ $i + 1 ] - 1 &&
                 $end ne '~' &&
                 $begins[ $i + 1 ] ne '~'
                ) {
                $begins[ $i + 1 ] = $end + 1;
                $cut_chain = 1;
            }
        }            

        if ( $cut_chain ) {
            log_append( "special handling - $id chain due to overlap\n" );
        }
    }
    
    ### extend last chain to end if short

    if ( $types[-1] eq 'chain' &&
         $ends[-1] < $fasta_res_count ) {
        $ends[-1] = $fasta_res_count;
        log_append( "special handling - $id last chain extended to match fasta count\n" );
    }

    ### special - remove chain 1 0 if present at top
    if ( $types[0] eq 'chain' &&
         $begins[0] == 1 &&
         $ends[0] == 0 &&
         $ends[0] ne '~' ) {
        ## remove
        shift @types;
        shift @begins;
        shift @ends;
        log_append( "special handling - $id first chain removed\n" );
    }

    ### special - cut 2nd propeptide or transit if covered by prior
    {
        my @ntypes  = @types [0];
        my @nbegins = @begins[0];
        my @nends   = @ends  [0];
        for ( my $i = 1; $i < @types; ++$i ) {
            ## previous types on stack

            my $ptype  = $ntypes [-1];
            my $pbegin = $nbegins[-1];
            my $pend   = $nends  [-1];

            ## this type

            my $type  = $types [$i];
            my $begin = $begins[$i];
            my $end   = $ends  [$i];

            if (
                ( $ptype eq 'propep' &&
                  $type  eq 'propep'  ) ||
                ( $ptype eq 'transit' &&
                  $type  eq 'transit' )
                ) {
                if ( 
                    $begin ne '~' &&
                    $end ne '~' &&
                    $pbegin ne '~' &&
                    $pend ne '~' &&
                    $end <= $pend &&
                    $begin >= $pbegin ) {
                    log_append( "special handling - $id sequential covered $type removed (a)\n" );
                    next;
                }
                if ( 
                    $begin ne '~' &&
                    $end ne '~' &&
                    $pbegin ne '~' &&
                    $pend ne '~' &&
                    $end > $pend &&
                    $pend >= $begin &&
                    $begin >= $pbegin ) {
                    log_append( "special handling - $id sequential covered $type removed (b)\n" );
                    $nends[-1] = $end;
                    next;
                }
            }

            push @ntypes,  $type;
            push @nbegins, $begin;
            push @nends,   $end;
        }
        @types  = @ntypes;
        @begins = @nbegins;
        @ends   = @nends;
    }

    ### overlap of transit to chain -> cut transit

    for ( my $i = 0; $i < @types; ++$i ) {
        my $type  = $types [$i];
        my $begin = $begins[$i];
        my $end   = $ends  [$i];

        my $cut_transit = 0;

        if ( $type eq 'transit' ) {
            if ( $i &&
                 $begin < $ends[ $i - 1 ] + 1 &&
                 $begin ne '~' &&
                 $ends[ $i - 1 ] ne '~' ) {
                $begins[ $i ] = $ends[ $i - 1 ] + 1;
                $cut_transit = 1;
            }
            if ( $i < @types - 1 &&
                 $end > $begins[ $i + 1 ] - 1 &&
                 $end ne '~' &&
                 $begins[ $i + 1 ] ne '~' ) {
                $ends[ $i ] = $begins[ $i + 1 ] - 1;
                $cut_transit = 1;
            }
        }            

        if ( $cut_transit ) {
            log_append( "special handling - $id cut transit due to overlap\n" );
        }
    }

    ### special handling - tilde merge

    my @t_begin = grep /^~$/, @begins;
    my @t_end   = grep /^~$/, @ends;

    if ( @t_begin || @_tend ) {
        my @ntypes  = @types [0];
        my @nbegins = @begins[0];
        my @nends   = @ends  [0];

        for ( my $i = 1; $i < @types; ++$i ) {
            ## previous types on stack

            my $ptype  = $ntypes [-1];
            my $pbegin = $nbegins[-1];
            my $pend   = $nends  [-1];

            ## this type

            my $type  = $types [$i];
            my $begin = $begins[$i];
            my $end   = $ends  [$i];

            if ( $pend eq '~' &&
                 $begin eq '~' ) {
                if ( $ptype eq $type ) {
                    log_append( "special handling - $id joined sequential $type using tildes\n" );
                    $nends [ -1 ] = $end;
                    next;
                }
                if ( $ptype ne 'chain' &&
                     $type  eq 'chain' ) {
                    log_append( "special handling - $id joined tidle-ending $ptype into chain\n" );
                    $ntypes[ -1 ] = "$ptype~chain";
                    $nends [ -1 ] = $end;
                    next;
                }
                if ( $ptype eq 'chain' &&
                     $type  ne 'chain' ) {
                    log_append( "special handling - $id joined tidle-ending $type into chain\n" );
                    $nends [ -1 ] = $end;
                    $ntypes[ -1 ] = "chain~$type";
                    next;
                }                    
                if ( $ptype eq 'signal' &&
                     $type  eq 'propep' ) {
                    log_append( "special handling - $id joined tidle-ending $ptype into propeptide\n" );
                    $ntypes[ -1 ] = 'signal~propep';
                    $nends [ -1 ] = $end;
                    next;
                }                    
                if ( $ptype eq 'transit' &&
                     $type  eq 'propep' ) {
                    log_append( "special handling - $id joined tidle-ending $ptype into propeptide\n" );
                    $ntypes[ -1 ] = 'transit~propep';
                    $nends [ -1 ] = $end;
                    next;
                }                    
            }       
            push @ntypes,  $type;
            push @nbegins, $begin;
            push @nends,   $end;
        }
        @types  = @ntypes;
        @begins = @nbegins;
        @ends   = @nends;
    }

    ## special handling - chain tilde start

    @t_begin = grep /^~$/, @begins;
    @t_end   = grep /^~$/, @ends;
    
    if ( @t_begin || @t_end ) {
        for ( my $i = 1; $i < @types; ++$i ) {
            ## previous types on stack

            my $ptype  = $types [ $i - 1 ];
            my $pbegin = $begins[ $i - 1 ];
            my $pend   = $ends  [ $i - 1 ];

            ## this type

            my $type  = $types [$i];
            my $begin = $begins[$i];
            my $end   = $ends  [$i];

            if ( $type eq 'chain' &&
                 $begin eq '~' &&
                 $pend ne '~' ) {
                $begins[ $i ] = $pend + 1;
                log_append( "special handling - $id tilde start chain set to prior end plus one\n" );
            }
        }
    }

    ## special handling - first end ~, next begin numbered
    if ( $ends[0] eq '~' &&
         @types > 1 &&
         $begins[1] ne '~' ) {
        $ends[0]   = $begins[1] - 1;
        $types[0] .= "~end_assumed";
        log_append( "special handling - first chain ended with ~, next numbered so end of first chain assumed\n" );
    }

    ## special handling - start at 1
    if ( $begins[ 0 ] != 1 ) {
        $begins[ 0 ] = 1;
        log_append( "special handling - extended 1st chain of type $types[0] to start at residue 1\n" );
    }
        
    ## various checks
    ## 1. overlaps or gaps
    ## 2. remaining tildes
    ## 3. fasta end matched
    ## 4. starts with 1 

    my $warn;
    warn_sig_init();
    
    my $tilde_msg = 0;
    
    for ( my $i = 0; $i < @types; ++$i ) {
        my $type  = $types [$i];
        my $begin = $begins[$i];
        my $end   = $ends  [$i];
        if ( $begin eq '~' ||
             $end eq '~' ) {
            if ( !$tilde_msg ) {
                $warn .= "warning - $id tildes present\n";
                warn_sig_add( "tilde" );
                $tilde_msg = 1;
            }
            next;
        }
        if ( $i &&
             $ends[ $i - 1 ] + 1 != $begin
            ) {
            $warn .= sprintf( "warning - $id gaps present end %d begin %d\n", $ends[ $i - 1 ] + 1 , $begin );
            warn_sig_add( "gap" );
            next;
        }
    }
    if ( $ends[-1] != $fasta_res_count ) {
        $warn .= sprintf( "warning - $id fasta sequence count %d does not match end %d\n", $fasta_res_count, $ends[-1] );
        warn_sig_add( "end" );
    }

    if ( $begins[ 0 ] != 1 ) {
        $warn .= "warning - $id does not start at residue 1 but starts at $begins[0]\n";
        warn_sig_add( "end" );
    }
        
    if ( $warn ) {
        log_append( line() . "*** warnings ***\n" . line() . $warn );
        log_flush();
    }

    ## assemble output
    my $out;
    for ( my $i = 0; $i < @types; ++$i ) {
        $out .= "$types[$i] $begins[$i] $ends[$i]\n";
    }

    log_append( line() . sprintf( "signature: %s" . ( join '+', @types ) . "\n", $peptides_only ? "peptides_only:" : "" ) );
    my $wsigs = warn_sig();
    if ( $wsigs ) {
        log_append( "warning_sig: $wsigs\n" );
    }

    ## write chain info
    my $outnotes;
    if ( @outnotes ) {
        grep s/^/# /, @outnotes;
        $outnotes = join "\n", @outnotes;
        $outnotes .= "\n";
    }

    log_append( line() . "final results\n" . line() . $outnotes . $out );

    write_file( $f, $outnotes . $out );

    ## we can change $cpnotes later if needed as a db field
    ## possibly split full log & notes ?

    $cpnotes = log_final();
    write_file( $fn, $cpnotes );
    log_append( "$0: wrote $fn\n" );
}

log_flush();
if ( $errors ) {
    print STDERR $errors;
}

exit length( $errors );    
    
