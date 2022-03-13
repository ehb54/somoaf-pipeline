#!/usr/bin/perl

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

takes pdbs and does further processing
ssbonds via somo
helix, sheet via chimera

";

die $notes if !@ARGV;

use File::Temp qw(tempfile);

require "$ppdir/pdbutil.pm";

p_config_init();

$errors = "";

sub map_tf {
    my $tf = shift; # || die "map_tf requires an argument\n";

    if ( $tf > 90 ) {
        return ( ( $tf - 90 ) * (20/10) ) + 80;
    } elsif ( $tf > 70 ) {
        return ( ( $tf - 70 ) * (30/20) ) + 50;
    } elsif ( $tf > 50 ) {
        return ( ( $tf - 50 ) * (25/20) ) + 25;
    } else {
        return ( ( $tf ) * (25/50) ) + 0;
    }
}

sub seqres {
    my $pdb = shift || error_exit( "$0: seqres() requires an argument" );

    my $seqs = \{};
    my %lastseq;
    
    my $lastresseq;
    my $lastchain;


    for my $l ( @$pdb ) {
        next if $l !~ /^ATOM/;
        my $resseq  = mytrim( substr( $l, 22, 4 ) );
        my $resname = substr( $l, 17, 3 );
        my $chain   = substr( $l, 21, 1 );

        if ( $resseq ne $lastresseq ) {
            $$seqs->{$chain} = [] if !exists $$seqs->{$chain};
            push @{$$seqs->{$chain}}, $resname;
            $lastseq{ $chain } = $resseq;
            $lastresseq        = $resseq;
        }
    }        

    my $line = 1;
    my @out;

    for my $c ( sort keys %{$$seqs} ) {
        ## build seq

        push @out, sprintf( "SEQRES%4d $c%5d ", $line++, $lastseq{$c} );
    
        my $added = 0;
        for my $res ( @{$$seqs->{$c}} ) {
            $out[-1] .= " $res";
            if ( ++$added == 13 ) {
                push @out, sprintf( "SEQRES%4d $c%5d ", $line++, $lastseq{$c} );
                $added = 0;
            }                
        }
        pop @out if !$added;
    }

    grep s/$/\n/, @out;
    @out;
}

while ( $fid = shift ) {
    my $id = uniprot_id( $fid );
    log_append( "$0: processing ascension $id", "-" );

    my @pdbs   = get_pdbs1( $id );

    if ( $debug > 1 ) {
        log_append( "$0: pdbs found", "-" );
        log_append( join "\n", @pdbs );
        log_append( "\n" );
    }

    for my $fpdb ( @pdbs ) {
        my $f           = "$$p_config{pdbstage1}/$fpdb";
        my $pdb_ver     = pdb_ver    ( $fpdb );
        my $pdb_frame   = pdb_frame  ( $fpdb );
        my $pdb_variant = pdb_variant( $fpdb );

        my $fo        = "$$p_config{pdbstage2}/AF-${id}-${pdb_frame}${pdb_variant}-model_${pdb_ver}-somo.pdb";

        log_append(
              "pdb_ver     '$pdb_ver'\n"
            . "pdb_frame   '$pdb_frame'\n"
            . "pdb_variant '$pdb_variant'\n"
            . "output      '$fo'\n"
            ) if $debug > 1;

        ## get ssbonds

        log_append( "$fpdb getting ssbonds\n" );

        my $cmd = qq[$$p_config{somocli} json '{"ssbond":1,"pdbfile":"$f"}' 2>/dev/null];
        my $res = run_cmd( $cmd, true );
        if ( run_cmd_last_error() ) {
            my $error = sprintf( "$0: ERROR [%d] - $fpdb running somo $cmd\n", run_cmd_last_error() );
            $errors .= $error;
        }
        $res =~ s/\n/\\n/g;
        my $dj = decode_json( $res );

        my @lpdb          = get_pdb_lines( $f, true );
        my @remarks       = grep /^REMARK   2 /, @lpdb;
        @lpdb             = grep !/^REMARK   2 /, @lpdb;
        @lpdb             = grep !/^SEQRES/, @lpdb;
        my @seqres        = seqres( \@lpdb );
        my $orgoxts       = scalar grep / OXT /, @lpdb;

        my @ssbonds;
        
        if ( length( $$dj{"ssbonds"} ) ) {
            my $ssbc = scalar split /\n/, $$dj{"ssbonds"};
            unshift @remarks, "REMARK   2 SSBOND ($ssbc) record(s) added\n";
            push @ssbonds, $$dj{"ssbonds"};
        } else {
            unshift @remarks, "REMARK   2 no SSBOND records added\n";
        }
        
        log_append( "$fpdb summary:\n"
                    . "remarks:\n" . ( join '', @remarks )
                    . "ssbonds:\n" . ( join '', @ssbonds )
                    , "+" );

        unshift @lpdb, @ssbonds;
        unshift @lpdb, @seqres;

        write_file( $fo, join '', @lpdb );

        ## run chimera
        log_append( "$fpdb running chimera\n" );

        my $mkchimera =
            "open $fo; addh; write format pdb 0 $fo; close all";

        my ( $fh, $ft ) = tempfile( "$$p_config{tempdir}/mkchimera.XXXXXX", UNLINK => 1 );
        print $fh $mkchimera;
        close $fh;
        run_cmd( "chimera --nogui < $ft", true );
        if ( run_cmd_last_error() ) {
            my $error = sprintf( "$0: ERROR [%d] - $fpdb running chimera $cmd\n", run_cmd_last_error() );
            $errors .= $error;
        }
        
        ## reread chimera file again, strip hydrogens
        @lpdb          = get_pdb_lines( $fo, true );
        for my $l ( @lpdb ) {
            next if $l !~/^ATOM/;
            $l = '' if mytrim( substr( $l, 76, 2 ) ) eq 'H';
        }
        write_file( $fo, join '', @lpdb );
        $mkchimera =
            "open $fo; write format pdb 0 $fo; close all";
        ( $fh, $ft ) = tempfile( "$$p_config{tempdir}/mkchimera.XXXXXX", UNLINK => 1 );
        print $fh $mkchimera;
        close $fh;
        run_cmd( "chimera --nogui < $ft", true );
        if ( run_cmd_last_error() ) {
            my $error = sprintf( "$0: ERROR [%d] - $fpdb running chimera $cmd\n", run_cmd_last_error() );
            $errors .= $error;
        }
        

        ## reread chimera file again
        @lpdb          = get_pdb_lines( $fo, true );
        my $helixl     = scalar grep /^HELIX/, @lpdb;
        my $sheetl     = scalar grep /^SHEET/, @lpdb;
        my $conectl    = scalar grep /^CONECT/, @lpdb;
        my $oxtl       = ( scalar grep / OXT /, @lpdb ) - $orgoxts;
        push @remarks, "REMARK   2 HELIX ($helixl), SHEET ($sheetl), CONECT ($conectl) records added using UCSF-CHIMERA\n";
        push @remarks, "REMARK   2 OXT ($oxtl) records added using UCSF-CHIMERA\n" if $oxtl;
        
        if ( grep /removed/, @remarks ) {
            unshift @remarks,
                "REMARK   2\n"
                . "REMARK   2 The following entries were added/removed by the US-SOMO team\n"
                ;
        } else {
            unshift @remarks,
                "REMARK   2\n"
                . "REMARK   2 The following entries were added by the US-SOMO team\n"
                ;
        }            
        push @remarks, "REMARK   2\n";

        log_append( "$fpdb remarks after chimera:\n"
                    . "remarks:\n" . ( join '', @remarks )
                    , "+" );

        for my $l (@lpdb ) {
            if ( $l =~ /^DBREF/ ) {
                $l = ( join '', @remarks ) . $l;
                last;
            }
        }

        write_file( $fo, join '', @lpdb );

        log_append( "$fpdb creating cif, mmcif\n" );
        ## pdb -> cif -> mmcif
        {
            my $fpdbnoext = $fpdb;
            $fpdbnoext    =~ s/\.pdb$//;
            my $cif       = "$$p_config{cifdir}/AF-${id}-${pdb_frame}${pdb_variant}-model_${pdb_ver}-somo.cif";
            my $mmcif     = "$$p_config{mmcifdir}/AF-${id}-${pdb_frame}${pdb_variant}-model_${pdb_ver}-somo.cif";
            my $logf      = "$$p_config{tempdir}/AF-${id}-${pdb_frame}${pdb_variant}-model_${pdb_ver}-somo.log";

            ## make cif
            {
                my $cmd = "$$p_config{maxit} -input $f -output $cif -o 1 -log $logf";
                run_cmd( $cmd, true );
                if ( run_cmd_last_error() ) {
                    my $error = sprintf( "$0: ERROR [%d] - $fpdb running maxit pdb->cif $cmd\n", run_cmd_last_error() );
                    $errors .= $error;
                }
                
            }

            ## make mmcif
            {
                my $cmd = "$$p_config{maxit} -input $cif -output $mmcif -o 8 -log $logf";
                run_cmd( $cmd, true );
                if ( run_cmd_last_error() ) {
                    my $error = sprintf( "$0: ERROR [%d] - $fpdb running maxit mmcif->cif $cmd\n", run_cmd_last_error() );
                    $errors .= $error;
                }
            }

            ## cleanup
            unlink $logf;
        }

        ## pdb -> pdb_tf 
        {
            log_append( "$fpdb creating pdb for jsmol visualization (tf encoded)\n" );
            my @ol;
            push @ol,
                "REMARK   0 **** WARNING: TF CONFIDENCE FACTORS ARE MODIFIED! ****\n"
                ."REMARK   0 **** THIS VERSION IS STRICTLY FOR JSMOL DISPLAY ****\n"
                ;

            my $count = 0;

            for my $l ( @lpdb ) {
                my $r = pdb_fields( $l );
                if ( $r->{"recname"}  =~ /^ATOM$/ ) {
                    my $tf;
                    if ( $count == 2 ) {
                        $tf = "0.00";
                    } elsif ( $count == 3 ) {
                        $tf = "100.00";
                    } else {
                        $tf = sprintf( "%.2f", 100 - map_tf( $r->{"tf"} ) );
                    }
                    $tf = ' 'x(6 - length($tf) ) . $tf;
                    $l = substr( $l, 0, 60 ) . $tf . substr($l, 66 );
                    ++$count;
                }
                push @ol, $l;
            }
            my $ftfo = "$$p_config{pdb_tfrev}/AF-${id}-${pdb_frame}${pdb_variant}-model_${pdb_ver}-somo.pdb";
            write_file( $ftfo, join '', @ol );
        }
    }
}

log_flush();
if ( $errors ) {
    print STDERR $errors;
}

exit length( $errors );    
    
