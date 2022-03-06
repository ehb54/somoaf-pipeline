#!/usr/bin/perl


## set to one for testing (skip CD & somo calcs)

$skipproc = 0;

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

takes finalized pdbs and perform computations

";

die $notes if !@ARGV;

use File::Temp qw(tempdir);
use File::Temp qw(tempfile);
use DateTime qw();

p_config_init();

require "$ppdir/mapping.pm";
require "$ppdir/pdbutil.pm";

$errors = "";

$processing_date = uc `date +'%d-%b-%y'`;
chomp $processing_date;

while ( $fid = shift ) {
    my $id = uniprot_id( $fid );
    log_append( "$0: processing ascension $id", "-" );

    my @pdbs   = get_pdbs2( $id );

    if ( $debug > 1 ) {
        log_append( "$0: pdbs found", "-" );
        log_append( join "\n", @pdbs );
        log_append( "\n" );
    }

    my $cdcmd = "python2 $$p_config{sesca} \@pdb";

    for my $fpdb ( @pdbs ) {
        my $f           = "$$p_config{pdbstage2}/$fpdb";
        my $pdb_ver     = pdb_ver    ( $fpdb );
        my $pdb_frame   = pdb_frame  ( $fpdb );
        my $pdb_variant = pdb_variant( $fpdb );
        my $fpdbnoext   = $fpdb;
        $fpdbnoext =~ s/\.pdb$//;

        log_append(
              "source file '$fpdb'\n"
            . "pdb_ver     '$pdb_ver'\n"
            . "pdb_frame   '$pdb_frame'\n"
            . "pdb_variant '$pdb_variant'\n"
            ) if $debug > 1;

        ## compute CD spectra

        if ( !$skipproc ) {
            my $template = "$$p_config{tempdir}/cdrun.XXXXXXXXX";
            my $dir = tempdir( $template, CLEANUP => 1 );

            my $fb  =  $f;
            $fb     =~ s/-somo\.pdb$//;

            log_append( "$fpdb compute CD spectra\n" );
            log_flush();
            my $cmd = "ln $f $dir/ && cd $dir && $cdcmd $f && grep -v Workdir: CD_comp.out | perl -pe 's/ \\/srv.*SESCA\\// SESCA\\//' > $$p_config{sescadir}/${fpdbnoext}-sesca-cd.dat";
            run_cmd( $cmd, true );
            if ( run_cmd_last_error() ) {
                my $error = sprintf( "$0: ERROR [%d] - $fpdb running SESCA computation $cmd\n", run_cmd_last_error() );
                $errors .= $error;
            }
        } else {
            warn "SESCA CD not computed!\n";
        }

        ## compute hydrodynamics && p(r)

        {
            my $mongoid = "AF-${id}-${pdb_frame}${pdb_variant}";

            ## somo cmds

            warn "WARNING - setup proper somo.config!\n";

            my ( $fh, $ft ) = tempfile( "$$p_config{tempdir}/somocmds.XXXXXX", UNLINK => 1 );
            print $fh
                "threads 1\n"
                . "batch selectall\n"
                . "batch somo_o\n"
                . "batch prr\n"
                . "batch zeno\n"
                . "batch combineh\n"
                . "batch combinehname $fpdb\n"
                . "batch saveparams\n"
                . "somo overwrite\n"
                . "batch overwrite\n"
                . "batch start\n"
                . "exit\n"
                ;
            close $fh;

            ## run somo
            
            my $prfile    = "$$p_config{somordir}/saxs/${fpdbnoext}_1b1.sprr_x";
            my $hydrofile = "$$p_config{somordir}/$fpdb.csv";

            my @expected_outputs =
                (
                 $hydrofile
                 ,$prfile
                );
            
            if ( !$skipproc ) {
                ## clean up before running

                unlink glob "$$p_config{somordir}/$fpdbnoext*";
                unlink glob "$$p_config{somordir}/saxs/$fpdbnoext*";

                my $cmd = "$$p_config{somoenv} && cd $$p_config{pdbstage2} && $$p_config{somorun} -g $ft $fpdb";
                run_cmd( $cmd, true, 2 ); # try 2x since very rarely zeno call crashes and/or hangs?

                ## cleanup extra files
                unlink glob "$$p_config{somordir}/$fpdbnoext*{asa_res,bead_model,hydro_res,bod}";

                ## check run was ok
                if ( run_cmd_last_error() ) {
                    my $error = sprintf( "$0: ERROR [%d] - $fpdb running SOMO computation $cmd\n", run_cmd_last_error() );
                    $errors .= $error;
                } else {
                    for my $eo ( @expected_outputs ) {
                        print "checking for: $eo\n";
                        if ( !-e $eo ) {
                            my $error = "$0: ERROR [%d] - $fpdb SOMO expected result $eo was not created";
                            $errors .= $error;
                            next;
                        }
                    }
                }

                ## rename and move p(r)

                {
                    my $cmd = "mv $prfile $$p_config{prdir}/${fpdbnoext}-pr.dat";
                    run_cmd( $cmd, true );
                    if ( run_cmd_last_error() ) {
                        my $error = sprintf( "$0: ERROR [%d] - $fpdb mv error $cmd\n", run_cmd_last_error() );
                        $errors .= $error;
                    }
                }
            }

            ## build up data for mongo insert

            my %data;

            ## extract csv info for creation of mongo insert

            die "$0: unexpected: $hydrofile does not exist\n" if !-e $hydrofile;
            
            my @hdata = `cat $hydrofile`;

            if ( @hdata != 2 ) {
                my $error = "$0: ERROR - $fpdb SOMO expected result $hydrofile does not contain 2 lines\n";
                $errors .= $error;
                next;
            }

            grep chomp, @hdata;
            
            ## split up csv and validate parameters
            {
                my @headers = split /,/, $hdata[0];
                my @params  = split /,/, $hdata[1];

                grep s/"//g, @headers;

                my %hmap = map { $_ => 1 } @headers;
                
                ## are all headers present?

                for my $k ( keys %csvh2mongo ) {
                    if ( !exists $hmap{$k} ) {
                        my $error = "$0: ERROR - $fpdb SOMO expected result $hydrofile does not contain header '$k'\n";
                        $errors .= $error;
                        next;
                    }
                }

                ## create data
                for ( my $i = 0; $i < @headers; ++$i ) {
                    my $h = $headers[$i];

                    ## skip any extra fields
                    next if !exists $csvh2mongo{$h};

                    $data{ $csvh2mongo{$h} } = $params[$i];
                }

            }

            ## additional fields
            $data{_id}      = "${id}-${pdb_frame}${pdb_variant}";
            $data{name}     = "AF-${id}-${pdb_frame}-model_${pdb_ver}";
            $data{somodate} = $processing_date;

            ### additional fields from the pdb
            {
                my @lpdb     = get_pdb_lines( "$$p_config{pdbstage2}/$fpdb" );

                {
                    my @lheaders = grep /^HEADER/, @lpdb;
                    if ( @lheaders != 1 ) {
                        my $error = "$0: ERROR - $fpdb pdb does not contain exactly one header line\n";
                        $errors .= $error;
                        next;
                    } else {
                        if ( $lheaders[0] =~ /HEADER\s*(\S+)\s*$/ ) {
                            $data{afdate} = $1;
                        } else {
                            $data{afdate} = "unknown";
                        }
                    }
                }
                {
                    my @lsource  = grep /^SOURCE/, @lpdb;
                    grep s/^SOURCE   .//, @lsource;
                    grep s/\s*$//, @lsource;
                    my $src = join '', @lsource;
                    if ( $src ) {
                        $data{source} = $src;
                    } else {
                        $data{source} = "unknown";
                    }
                }
                {
                    my @ltitle  = grep /^TITLE/, @lpdb;
                    grep s/^TITLE   ..//, @ltitle;
                    grep s/\s*$//, @ltitle;
                    my $title = join '', @ltitle;
                    $title =~ s/^\s*//;
                    if ( $title ) {
                        $data{title} = $title;
                    } else {
                        $data{title} = "unknown";
                    }
                }
                {
                    my @lremarks = grep /^REMARK   2 /, @lpdb;
                    {
                        my @res = grep /^REMARK   2 RESIDUE seq\. /, @lremarks;
                        if ( @res != 1 ) {
                            my $error = "$0: ERROR - $fpdb pdb does not contain exactly one RESIDUE seq. line\n";
                            $errors .= $error;
                            next;
                        } else {
                            my $res = $res[0];
                            $res =~ s/^.*seq\. //;
                            $data{res} = $res;
                        }
                    }
                    {
                        my @removed = grep /removed/, @lremarks;
                        @removed = grep !/US-SOMO/, @removed;
                        grep s/^REMARK   2 UNIPROT //, @removed;
                        grep s/ removed.*$//, @removed;
                        if ( @removed ) {
                            my $removed = join ', ', @removed;
                            $removed .= " removed";
                            $data{proc} = $removed;
                        } else {
                            $data{proc} = "none";
                        }                            
                    }
                }

                #### helix/sheet

                {
                    my $lastresseq = 0;
                    my $helixcount = 0;
                    my $sheetcount = 0;

                    for my $l ( @lpdb ) {
                        my $r = pdb_fields( $l );
                        my $recname = $r->{recname};
                        if ( $recname =~ /^HELIX/ ) {
                            my $initseqnum = $r->{initseqnum};
                            my $endseqnum  = $r->{endseqnum};
                            $helixcount += $endseqnum - $initseqnum;
                            next;
                        } elsif ( $recname =~ /^SHEET/ ) {
                            my $initseqnum = $r->{initseqnum};
                            my $endseqnum  = $r->{endseqnum};
                            $sheetcount += $endseqnum - $initseqnum;
                            next;
                        } elsif ( $recname =~ /^ATOM/ ) {
                            my $resseq = $r->{resseq};
                            if ( $lastresseq != $resseq ) {
                                $lastresseq = $resseq;
                                ++$rescount;
                            }
                        }
                    }

                    $data{helix} = sprintf( "%.2f", $helixcount * 100.0 / ( $rescount - 1.0 ) );
                    $data{sheet} = sprintf( "%.2f", $sheetcount * 100.0 / ( $rescount - 1.0 ) );
                }

                #### confidence
                {
                    my $count = 0;
                    my $total = 0;

                    my $lastresseq = 0;

                    for my $l ( @lpdb ) {
                        my $r = pdb_fields( $l );
                        if ( $r->{recname}  =~ /^ATOM$/ ) {
                            my $resseq = $r->{resseq};
                            if ( $lastresseq != $resseq ) {
                                $lastresseq = $resseq;
                                $total += $r->{tf};
                                ++$count;
                            }
                        }
                    }

                    $data{afmeanconf} = sprintf( "%.2f", $total / $count ) + 0;
                }
                    
            }

            ## build mongo command

            {
                my $mongoc = qq[db.$$p_config{mongocoll}.update({_id:"$data{_id}"},{\$set:{];
                my @sets;
                for my $k ( keys %data ) {
                    if ( exists $mongostring{$k} ) {
                        push @sets, qq[$k:"$data{$k}"];
                    } else {
                        push @sets, qq[$k:$data{$k}];
                    }
                }
                $mongoc .= join ',', @sets;
                $mongoc .= qq[}},{upsert:true})\n];
                # print $mongoc;
                write_file( "$$p_config{mongocmds}/$data{_id}.mongo", $mongoc );
            }

            ## add data to mongo

            {
                my $cmd = "mongo $$p_config{mongodb} < $$p_config{mongocmds}/$data{_id}.mongo";
                run_cmd( $cmd, true );
                if ( run_cmd_last_error() ) {
                    my $error = sprintf( "$0: ERROR [%d] - $fpdb inserting into mongo $cmd\n", run_cmd_last_error() );
                    $errors .= $error;
                }
            }
            
            ## get csv and save to file

            {
                my $csv  = dbm_csv_header() . dbm_csv( $data{_id} );
                my $csvf = "$$p_config{csvdir}/${fpdbnoext}.csv";
                write_file( $csvf, $csv );
            }

            if ( $debug > 2 ) {
                print Dumper( \%data );
                print "mongoid: $mongoid\n";
                die "testing\n";
            }

            ## make tar & zip files
            {
                my $template = "$$p_config{tempdir}/zip.XXXXXXXXX";
                my $dir = tempdir( $template, CLEANUP => 1 );
                ### link contents into tar directory
                my $cmd =
                    "cd $dir"
                    . " && ln $$p_config{pdbstage2}/$fpdb $$p_config{mmcifdir}/${fpdbnoext}.cif $$p_config{prdir}/${fpdbnoext}-pr.dat $$p_config{sescadir}/${fpdbnoext}-sesca-cd.dat $$p_config{csvdir}/${fpdbnoext}.csv ."
                    . " && tar Jcf $$p_config{txzdir}/${fpdbnoext}.txz *"
                    . " && zip $$p_config{zipdir}/${fpdbnoext}.zip *"
                    ;
                run_cmd( $cmd, true );
                if ( run_cmd_last_error() ) {
                    my $error = sprintf( "$0: ERROR [%d] - $fpdb error creating txz & zips $cmd\n", run_cmd_last_error() );
                    $errors .= $error;
                }
            }
        }
    }
}

log_flush();
if ( $errors ) {
    print STDERR $errors;
}

exit length( $errors );    
