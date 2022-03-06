# perl defines and utility module for alphafold

use JSON;
use PHP::Functions::File qw(file_get_contents file_put_contents);
use MongoDB;
use Data::Dumper;
use HTTP::Tiny;

die "$0: \$ppdir is not defined\n" if !defined $ppdir;

# configfile
$configfile = "$ppdir/config.json";

# local storage

sub dbm_connect {
    print STDERR "dbm_connect\n" if $debug;
    $mongo_client = MongoDB->connect();
    $mongo_af      = $mongo_client->ns( $__dbm_coll );
}

sub dbm_read {
    my $name   = shift || die "dbm_read requires argument\n";

#    print "dbm_read : '$name'\n";

    dbm_connect() if !defined $mongo_client;
    my $result = $mongo_af->find_one( {"_id" => $name } );
    $result;
}

sub dbm_write {
    my $name    = shift || die "dbm_write requires argument\n";
    my $payload = shift || die "dbm_write requires 2 arguments\n";

    dbm_connect() if !defined $mongo_client;
    $$payload{ "_id" } = $name;
#    print "ref is " . ref( $payload ) . "\n";
#    print "dbm_write : '$name'\n" . Dumper( $payload );
    
    eval {
        $result =
            $mongo_af->insert_one( $payload );
        print "inserted id " . $result->inserted_id . "\n" if $debug;
        1;
    } or do {
        print " error after eval ... " .  Dumper( $@ ) . "\n";
        return false;
    };
    if ( $@ ) {
        print "error\n";
        return false;
    }
    return true;
}

sub dbm_count {
    dbm_connect() if !defined $mongo_client;

    return $mongo_af->count_documents({});
}

sub dbm_print {
    my $name = shift;
    my $res = dbm_read( $name );
    if ( !$res ) {
        print "$name does not exist in the db\n";
        return;
    }
    print
        sprintf(
            "_id                                                      : %s\n"
            . "Model name                                               : %s\n"
            . "Title                                                    : %s\n"
            . "Source                                                   : %s\n"
            . "Alphafold date                                           : %s\n"
            . "Hydrodynamic calculations date                           : %s\n"
            . "Post translational processing                            : %s\n"
            . "UniProt residues present                                 : %s\n"
            . "Molecular mass [Da]                                      : %.2f\n"
            . "Partial specific volume [cm^3/g]                         : %.3f\n"
            . "Sedimentation coefficient s [S]                          : %.5g\n"
            . "Sedimentation coefficient s.d.                           : %.5g\n"
            . "Translational diffusion coefficient D [cm/sec^2]         : %.5g\n"
            . "Translational diffusion coefficient D s.d.               : %.5g\n"
            . "Stokes radius [nm]                                       : %.5g\n"
            . "Stokes radius s.d.                                       : %.5g\n"
            . "Intrinsic viscosity [cm^3/g]                             : %.5g\n"
            . "Intrisic viscosity s.d.                                  : %.5g\n"
            . "Radius of gyration (+r) [A] (from PDB atomic structure)  : %.2f\n"
            . "Maximum extensions X [nm]                                : %.2f\n"
            . "Maximum extensions Y [nm]                                : %.2f\n"
            . "Maximum extensions Z [nm]                                : %.2f\n"
            . "Helix %                                                  : %.2f\n"
            . "Sheet %                                                  : %.2f\n"

            , $$res{ '_id' }
            , $$res{ 'name' }
            , $$res{ 'title' }
            , $$res{ 'source' }
            , $$res{ 'afdate' }
            , $$res{ 'somodate' }
            , $$res{ 'proc' }
            , $$res{ 'res' }
            , $$res{ 'mw' }
            , $$res{ 'psv' }
            , $$res{ 'S' }
            , $$res{ 'S_sd' }
            , $$res{ 'Dtr' }
            , $$res{ 'Dtr_sd' }
            , $$res{ 'Rs' }
            , $$res{ 'Rs_sd' }
            , $$res{ 'Eta' }
            , $$res{ 'Eta_sd' }
            , $$res{ 'Rg' }
            , $$res{ 'ExtX' }
            , $$res{ 'ExtY' }
            , $$res{ 'ExtZ' }
            , $$res{ 'helix' }
            , $$res{ 'sheet' }
        );
    "";
}

sub dbm_csv_header {
    '"UniProt accession"'
        . ',"AlphaFold Model name"'
        . ',"Title"'
        . ',"Source"'
        . ',"Alphafold date"'
        . ',"Mean confidence"'
        . ',"Hydrodynamic calculations date"'
        . ',"Post translational processing"'
        . ',"UniProt residues present"'
        . ',"Molecular mass [Da]"'
        . ',"Partial specific volume [cm^3/g]"'
        . ',"Translational diffusion coefficient D [F]"'
#        . ',"Translational diffusion coefficient D s.d."'
        . ',"Sedimentation coefficient s [S]"'
#        . ',"Sedimentation coefficient s.d."'
        . ',"Stokes radius [nm]"'
#        . ',"Stokes radius s.d."'
        . ',"Intrinsic viscosity [cm^3/g]"'
        . ',"Intrisic viscosity s.d."'
        . ',"Radius of gyration (+r) [A] (from PDB atomic structure)"'
        . ',"Maximum extensions X [nm]"'
        . ',"Maximum extensions Y [nm]"'
        . ',"Maximum extensions Z [nm]"'
        . ',"Helix %"'
        . ',"Sheet %"'
        . "\n"
        ;
}

sub dbm_csv {
    my $name = shift;
    my $res = dbm_read( $name );
    if ( !$res ) {
        print "$name does not exist in the db\n";
        return;
    }
     sprintf(
         qq["%s","%s","%s","%s","%s","%s",%s,"%s","%s",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n]
        , $$res{ '_id' }
        , $$res{ 'name' }
        , $$res{ 'title' }
        , $$res{ 'source' }
        , $$res{ 'afdate' }
        , sprintf( "%.2f", $$res{ 'afmeanconf' } )
        , $$res{ 'somodate' }
        , $$res{ 'proc' }
        , $$res{ 'res' }
        , sprintf( "%.1f", $$res{ 'mw' } )
        , $$res{ 'psv' }
        , sprintf( "%.3g", $$res{ 'Dtr' } * 1e7 )
#        , sprintf( "%.2g", $$res{ 'Dtr_sd' } * 1e7 )
        , sprintf( "%.3g", $$res{ 'S' } )
#        , sprintf( "%.2g", $$res{ 'S_sd' } )
        , sprintf( "%.3g", $$res{ 'Rs' } )
#        , sprintf( "%.2g", $$res{ 'Rs_sd' } )
        , sprintf( "%.3g", $$res{ 'Eta' } )
        , sprintf( "%.2f", $$res{ 'Eta_sd' } )
        , sprintf( "%.3g", $$res{ 'Rg' } )
        , sprintf( "%.2f", $$res{ 'ExtX' } )
        , sprintf( "%.2f", $$res{ 'ExtY' } )
        , sprintf( "%.2f", $$res{ 'ExtZ' } )
        , sprintf( "%.1f", $$res{ 'helix' } )
        , sprintf( "%.1f", $$res{ 'sheet' } )
        );
}

sub dbm_csv_no_multiframe {
    my $name = shift;
    my $res = dbm_read( $name );
    if ( !$res ) {
        print "$name does not exist in the db\n";
        return;
    }

    return "" if $$res{ 'multiframe' };

    sprintf(
         qq["%s","%s","%s","%s","%s",%s,"%s","%s",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n]
        , $$res{ '_id' }
        , $$res{ 'name' }
        , $$res{ 'title' }
        , $$res{ 'source' }
        , $$res{ 'afdate' }
        , sprintf( "%.2f", $$res{ 'afmeanconf' } )
        , $$res{ 'somodate' }
        , $$res{ 'sp' } ? $$res{ 'sp' } : "n/a"
        , sprintf( "%.1f", $$res{ 'mw' } )
        , $$res{ 'psv' }
        , sprintf( "%.3g", $$res{ 'Dtr' } * 1e7 )
#        , sprintf( "%.2g", $$res{ 'Dtr_sd' } * 1e7 )
        , sprintf( "%.3g", $$res{ 'S' } )
#        , sprintf( "%.2g", $$res{ 'S_sd' } )
        , sprintf( "%.3g", $$res{ 'Rs' } )
#        , sprintf( "%.2g", $$res{ 'Rs_sd' } )
        , sprintf( "%.3g", $$res{ 'Eta' } )
        , sprintf( "%.2f", $$res{ 'Eta_sd' } )
        , sprintf( "%.3g", $$res{ 'Rg' } )
        , sprintf( "%.2f", $$res{ 'ExtX' } )
        , sprintf( "%.2f", $$res{ 'ExtY' } )
        , sprintf( "%.2f", $$res{ 'ExtZ' } )
        , sprintf( "%.1f", $$res{ 'helix' } )
        , sprintf( "%.1f", $$res{ 'sheet' } )
        );
}

sub dbm_ids {
    dbm_connect() if !defined $mongo_client;
    my @ids;
    my $print = shift;

    my $cursor = $mongo_af->find();

    # Cursor iteration
    while ( my $doc = $cursor->next ) {
        push @ids, $$doc{'_id'};
        dbm_print( $$doc{'_id'} ) if $print;
    }
    return \@ids
}

$run_cmd_last_error;

sub run_cmd {
    my $cmd       = shift || die "run_cmd() requires an argument\n";
    my $no_die    = shift;
    my $repeattry = shift;
    print "$cmd\n" if $debug;
    $run_cmd_last_error = 0;
    my $res = `$cmd`;
    if ( $? ) {
        $run_cmd_last_error = $?;
        if ( $no_die ) {
            warn "run_cmd(\"$cmd\") returned $?\n";
            if ( $repeattry > 0 ) {
                warn "run_cmd(\"$cmd\") repeating failed command tries left = $repeattry )\n";
                return run_cmd( $cmd, $no_die, --$repeattry );
            }
        } else {
            error_exit( "run_cmd(\"$cmd\") returned $?" );
        }
    }
                
    chomp $res;
    return $res;
}

sub run_cmd_last_error {
    return $run_cmd_last_error;
}

sub dbm_csv_header_spec {
    '"Model Name"'
        . ',"Molecular mass [Da]"'
        . ',"Partial specific volume [cm^3/g]"'
        . ',"Translational diffusion coefficient D [F]"'
        . ',"Sedimentation coefficient s [S]"'
        . ',"Stokes radius [nm]"'
        . ',"Intrinsic viscosity [cm^3/g]"'
        . ',"Intrisic viscosity s.d."'
        . ',"Radius of gyration (+r) [A] (from PDB atomic structure)"'
        . ',"Rg/Rs"'
        . ',"Maximum extensions X [nm]"'
        . ',"Maximum extensions Y [nm]"'
        . ',"Maximum extensions Z [nm]"'
        . "\n"
        ;
}

sub stats_init {
    undef %stats;
}

sub stats_add_val {
    my $key = shift;
    my $val = shift;

    $stats{ "$key:sum" }  += $val;
    $stats{ "$key:sum2" } += $val * $val;

    if ( $stats{ 'count' } ) {
        $stats{ "$key:min" } = $stats{ "$key:min" } > $val ? $val : $stats{ "$key:min" };
        $stats{ "$key:max" } = $stats{ "$key:max" } < $val ? $val : $stats{ "$key:max" };
    } else {
        $stats_keys{ $key }++;
        $stats{ "$key:min" } = $val;
        $stats{ "$key:max" } = $val;
    }
}

sub stats_commit {
    $stats{ 'count' }++;
}


sub dbm_csv_spec {
    my $name = shift;
    my $res = dbm_read( $name );
    if ( !$res ) {
        print "$name does not exist in the db\n";
        return;
    }

    stats_add_val( 'mw', sprintf( "%.1f", $$res{ 'mw' } ) );
    stats_add_val( 'psv', $$res{ 'psv' } );
    stats_add_val( 'Dtr', sprintf( "%.3g", $$res{ 'Dtr' } * 1e7 ) );
    stats_add_val( 'S', sprintf( "%.3g", $$res{ 'S' } ) );
    stats_add_val( 'Rs', sprintf( "%.3g", $$res{ 'Rs' } ) );
    stats_add_val( 'Eta', sprintf( "%.3g", $$res{ 'Eta' } ) );
    stats_add_val( 'Eta_sd', sprintf( "%.2f", $$res{ 'Eta_sd' } ) );
    stats_add_val( 'Rg', sprintf( "%.3g", $$res{ 'Rg' } ) );
    stats_add_val( 'Rg_Rs', sprintf( "%.3g", $$res{ 'Rg_Rs' } ) );
    stats_add_val( 'ExtX', sprintf( "%.2f", $$res{ 'ExtX' } ) );
    stats_add_val( 'ExtY', sprintf( "%.2f", $$res{ 'ExtY' } ) );
    stats_add_val( 'ExtZ', sprintf( "%.2f", $$res{ 'ExtZ' } ) );
                           
    stats_commit();

    sprintf(
         qq["%s",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n]
        , $$res{ '_id' }
        , sprintf( "%.1f", $$res{ 'mw' } )
        , $$res{ 'psv' }
        , sprintf( "%.3g", $$res{ 'Dtr' } * 1e7 )
        , sprintf( "%.3g", $$res{ 'S' } )
        , sprintf( "%.3g", $$res{ 'Rs' } )
        , sprintf( "%.3g", $$res{ 'Eta' } )
        , sprintf( "%.2f", $$res{ 'Eta_sd' } )
        , sprintf( "%.3g", $$res{ 'Rg' } )
        , sprintf( "%.3g", $$res{ 'Rg_Rs' } )
        , sprintf( "%.2f", $$res{ 'ExtX' } )
        , sprintf( "%.2f", $$res{ 'ExtY' } )
        , sprintf( "%.2f", $$res{ 'ExtZ' } )
        );
}

sub stats_compute_final {
    for my $k ( keys %stats_keys ) {
        $stats{ "$k:sd"  } = sqrt( abs( $stats{ "$k:sum2" } - ( ( $stats{ "$k:sum" } * $stats{ "$k:sum" } ) / $stats{ 'count' } ) ) / ( $stats{ 'count' } - 1 ) );
        $stats{ "$k:avg" } = $stats{ "$k:sum" } / $stats{ 'count' };
    }
}

sub dbm_csv_spec_summary {
    my $out;
    for my $kt ( 'avg', 'sd', 'min', 'max' ) {
        $out .=
            sprintf(
                qq["%s",%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n]
                , $kt
                , sprintf( "%.2f", $stats{ "mw:$kt" } )
                , $stats{ "psv:$kt" }
                , sprintf( "%.4g", $stats{ "Dtr:$kt" } )
                , sprintf( "%.4g", $stats{ "S:$kt" } )
                , sprintf( "%.4g", $stats{ "Rs:$kt" } )
                , sprintf( "%.4g", $stats{ "Eta:$kt" } )
                , sprintf( "%.3f", $stats{ "Eta_sd:$kt" } )
                , sprintf( "%.4g", $stats{ "Rg:$kt" } )
                , sprintf( "%.4g", $stats{ "Rg_Rs:$kt" } )
                , sprintf( "%.3f", $stats{ "ExtX:$kt" } )
                , sprintf( "%.3f", $stats{ "ExtY:$kt" } )
                , sprintf( "%.3f", $stats{ "ExtZ:$kt" } )
            );
    }
    $out;
}

sub p_config_init {
    error_exit("$0: p_config_init() : $configfile does not exist or is not readable" ) if !-e $configfile || !-r $configfile;

    my @configdata = `cat $configfile`;

    error_exit( "$0: p_config_init(): $configfile empty" ) if @$configdata;

    grep s/#.*$//, @configdata;

    $p_config = decode_json( join '', @configdata ) || die "$0: $configfile error decoding\n";

    for my $req (
        "base"
        ,"afversion"
        ,"afftp"
        ,"afpdb"
        ,"fasta"
        ,"fasta_url"
        ,"unifeat"
        ,"unifeat_url"
        ,"chains"
        ,"cpnotes"
        ,"cperrs"
        ,"pdbstage1"
        ,"pdbstage2"
        ,"pdb_tfrev"
        ,"somocli"
        ,"sesca"
        ,"sescadir"
        ,"somoenv"
        ,"somorun"
        ,"somordir"
        ,"tempdir"
        ,"cifdir"
        ,"mmcifdir"
        ,"csvdir"
        ,"prdir"
        ,"zipdir"
        ,"txzdir"
        ,"maxit"
        ,"mongodb"
        ,"mongocoll"
        ,"mongocmds"
        ,"stages"
        ) {
        error_exit( "$0: p_config_init() : $configfile missing required definitions for $req" ) if !exists $$p_config{$req};
    }
    for my $d (
        "base"
        ,"afftp"
        ,"afpdb"
        ,"fasta"
        ,"unifeat"
        ,"chains"
        ,"cperrs"
        ,"pdbstage1"
        ,"pdbstage2"
        ,"pdb_tfrev"
        ,"sescadir"
        ,"somordir"
        ,"tempdir"
        ,"cifdir"
        ,"mmcifdir"
        ,"csvdir"
        ,"prdir"
        ,"zipdir"
        ,"txzdir"
        ,"mongocmds"
        ) {
        error_exit( "$0: p_config_init() : $configfile $d $$p_config{$dir} is not a directory" ) if !-d $$p_config{$d};
    }

    $__dbm_coll = "$$p_config{mongodb}.$$p_config{mongocoll}";

    log_init();
    if ( exists $$p_config{debug} ) {
        $debug = $$p_config{debug};
    }
}

sub error_exit {
    my $msg = shift;
    die "$msg\n";
}

sub line {
    my $char = shift;
    $char = '-' if !$char;
    ${char}x80 . "\n";
}


## logging utilities

sub log_init {
    $log         = "";
    $log_flushed = "";
}

sub log_append {
    my $msg = shift;
    my $lchar = shift;
    if ( $lchar ) {
        $log .= line( $lchar ) . $msg . "\n" . line( $lchar );
    } else {
        $log .= $msg;
    }
}

sub log_flush {
    if ( $debug ) {
        $log_flushed .= $log;
        print $log;
        $log = "";
    }
}

sub log_final {
    if ( $debug ) {
        log_flush();
        return $log_flushed;
    } else {
        $log;
    }
}

sub log_write {
    error_exit( "log_write() - not yet" );
    my $fo = shift;
    open OUT, ">$fo" || error_exit( "$0: log_write() error trying write log $fo $!" );
    print OUT log_final();
    close OUT;
}

sub uniprot_id {
    my $id = shift || error_exit( "$0: uniprot_id() requires an argument" );
    $id =~ s/^AF-//;
    $id =~ s/\..*$//;
    $id =~ s/-.*$//;
    $id;
}

sub af_pdbs {
    my $id = shift || error_exit( "$0: af_pdbs() requires an argument" );
    my @f = `cd $$p_config{afpdb} && ls AF-$id-F*-model_v*.pdb.gz`;
    grep chomp, @f;
    @f;
}

sub get_chains {
    my $id = shift || error_exit( "$0: get_chains() requires an argument" );
    my $f  = "$$p_config{chains}/$id.chains";
    error_exit( "$0: get_chains($id) chains file not found" ) if !-e $f;
    my @chains = `cat $f`;
    grep chomp, @chains;
    @chains;
}

sub get_fasta_residue_count {
    my $id = shift || error_exit( "$0: get_fasta_residue_count() requires an argument" );
    my $f = "$$p_config{fasta}/$id.fasta";
    error_exit( "$0: get_fasta_residue_count( $id ) $f does not exist" ) if !-e $f;

    my @seq = `sed '1d' $f`;
    grep chomp, @seq;
    length( join '', @seq );
}

sub get_pdb_residue_count {
    my $f = shift ||  error_exit( "$0: get_pdb_residue_count() requires an argument" );
    $f = "$$p_config{afpdb}/$f";
    error_exit( "$0: get_pdb_residue_count( $f ) pdb not found" ) if !-e $f;
    my $res = `zcat $f | grep ^ATOM | tail -1  | cut -c 23-27 | xargs`;
    chomp $res;
    return $res;
}

sub get_pdb_lines {
    my $f = shift ||  error_exit( "$0: get_pdb_lines() requires an argument" );
    my $nochomp = shift;
    $f = "$$p_config{afpdb}/$f" unless $f =~ /^\//;
    error_exit( "$0: get_pdb_lines( $f ) pdb not found" ) if !-e $f;
    my $cmd = $f =~ /\.gz$/ ? 'zcat' : 'cat';
    my @res = `$cmd $f`;
    grep chomp, @res if !$nochomp;
    @res;
}

sub get_pdbs1 {
    my $id = shift || error_exit( "$0: get_pdbs2() requires an argument" );
    my @f = `cd $$p_config{pdbstage1} && ls $id-F*-v*.pdb`;
    grep chomp, @f;
    if ( $$p_config{afversion} ne "all" ) {
        @f = grep /v$$p_config{afversion}\./,@f;
    }
    @f;
}

sub get_pdbs2 {
    my $id = shift || error_exit( "$0: get_pdbs2() requires an argument" );
    my @f = `cd $$p_config{pdbstage2} && ls AF-${id}-*.pdb`;
    grep chomp, @f;
    if ( $$p_config{afversion} ne "all" ) {
        @f = grep /v$$p_config{afversion}-somo\./,@f;
    }
    @f;
}

sub get_residue_count {
    my $id = shift || error_exit( "$0: get_residue_count() requires an argument" );
    my @fs = af_pdbs( $id );
    error_exit( "$0: get_residue_count() no pdbs found" ) if !@fs;
    if ( @fs == 1 ) {
        my $res = `zcat $$p_config{afpdb}/$fs[0] | grep ^ATOM | tail -1  | cut -c 23-27 | xargs`;
        chomp $res;
        return $res;
    }
    ## get max frame #
    my $maxframe = 1;
    my @fsf = @fs;
    grep s/(^.*F|-model.*$)//g, @fsf;
    @fsf = sort { $a <=> $b } @fsf;
    my $lframe = $fsf[-1];
    my @fs = grep /-F${lframe}-model/, @fs;
    if ( $lframe > 1 ) {
        return "multiframe";
    }

    if ( @fs == 1 ) {
        my $res = `zcat $$p_config{afpdb}/$fs[0] | grep ^ATOM | tail -1  | cut -c 23-27 | xargs`;
        chomp $res;
        return $res;
    }

    ## must be multiple versions (likely most common case)
    ## check them all

    my %rcnts;
    
    for my $f ( @fs ) {
        my $res = `zcat $$p_config{afpdb}/$f | grep ^ATOM | tail -1  | cut -c 23-27 | xargs`;
        chomp $res;
        $rcnts{$res}++;
    }
    
    if ( scalar keys %rcnts == 1 ) {
        return (keys %rcnts)[0];
    }
    error_exit( "$0 : get_residue_count( $id ) differing version residue counts" );
}
   
sub pdb_ver {
    my $pdb = shift || error_exit( "$0: pdb_ver() requires an argument" );
    $pdb =~ /^.*(v\d+)(?:|-somo)\.pdb/;
    $1;
}

sub pdb_frame {
    my $pdb = shift || error_exit( "$0: pdb_frame() requires an argument" );
    $pdb =~ /^.*-(F\d+)-.*\.pdb/;
    $1;
}

sub pdb_variant {
    my $pdb = shift || error_exit( "$0: pdb_variant() requires an argument" );
    $pdb =~ /^.*-F\d+(-[^-]*).*-(?:|model_)v.*\.pdb/;
    $1;
}

sub read_features {
    my $id = shift || error_exit( "$0: read_features() requires an argument" );
    my $f = "$$p_config{unifeat}/${id}.uf.txt";
    error_exit( "$0: read_features( $id ) $f does not exist" ) if !-e $f;
    return `cat $f`;
}

sub warn_sig_init {
    undef %warn_sigs;
}

sub warn_sig_add {
    my $type = shift;
    $warn_sigs{ $type }++;
}

sub warn_sig {
    return "" if !%warn_sigs;
    join '+', sort keys %warn_sigs;
}

sub write_file {
    my $f   = shift || error_exit( "$0: write_file() : missing argument" );
    my $msg = shift || error_exit( "$0: write_file( $f ) : missing 2nd argument" );
    open my $fh, ">$f" || error_exit( "$0: write_file( $f, _ ) : file open error $!" );
    print $fh $msg;
    close $fh;
    error_exit( "$0: error writing file $f, does not exist after writing!\n" ) if !-e $f;
}

sub debug_json {
    my $tag = shift;
    my $msg = shift;
    my $json = JSON->new; # ->allow_nonref;
    
    line()
        . "$tag\n"
        . line()
        . $json->pretty->encode( $$msg )
        . "\n"
        . line()
        ;
}

return true;
