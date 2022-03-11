#!/usr/bin/perl

my $targetbase  = "/scratch/00451/tg457210/af/w";
my $resultsbase = "/work/00451/tg457210/ls6";

$notes = "
usage: $0 ids

id is a number of the split job

copy results unified to work


";

die $notes if !@ARGV;

use JSON;

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
        ,"cpnotes"
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
        error_exit( "$0: p_config_init() : $configfile $d $$p_config{$d} is not a directory" ) if !-d $$p_config{$d};
    }

    $__dbm_coll = "$$p_config{mongodb}.$$p_config{mongocoll}";

    log_init();

    if ( exists $$p_config{debug} ) {
        $debug = $$p_config{debug};
    }

    if ( exists $$p_config{"ln"} ) {
        $ln = $$p_config{"ln"};
    } else {
        $ln = "ln";
    }
}

sub error_exit {
    my $msg = shift;
    die "$msg\n";
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

sub line {
    my $char = shift;
    $char = '-' if !$char;
    ${char}x80 . "\n";
}

while ( @ARGV ) {
    my $id = shift;
    
    my $uid   = '0'x(3 - length($id)) . $id;

    ## copy results for $uid to work base

    my $base    = "$targetbase/w_$uid";
    $configfile = "$base/afpp/config.json";
    p_config_init();

    my $ids     = "$base/ids";

    error_exit( "no $ids found" ) if !-e $ids;

    my @ids = `cat $ids`;

    my $f           = "$$p_config{pdbstage2}";
    my $fpdbtf      = "$$p_config{pdb_tfrev}";
    my $fcd         = "$$p_config{sescadir}";
    my $fmmcif      = "$$p_config{mmcifdir}";
    my $fcsv        = "$$p_config{csvdir}";
    my $fpr         = "$$p_config{prdir}";
    my $fzip        = "$$p_config{zipdir}";
    my $ftxz        = "$$p_config{txzdir}";
    my $fmongo      = "$$p_config{mongocmds}";

    my %exists = (
        $f              => "pdb"
        ,$fpdbtf        => "pdb_tfrev"
        ,$fcd           => "cd"
        ,$fmmcif        => "mmcif"
        ,$fcsv          => "csv"
        ,$fpr           => "pr"
        ,$fzip          => "zip"
        ,$ftxz          => "txz"
        ,$fmongo        => "mongo"
        );

    
    my @pkgdirs;

    my $any_errors = 0;

    for my $check ( keys %exists ) {
        if ( !-d $check ) {
            my $error = "$0: ERROR - $check is missing\n";
            $errors .= $error;
            $any_errors = 1;
        }
        push @pkgdirs, $exists{$check} if $exists{$check} !~ /^>/;
    }
    next if $any_errors;


    {
        #        my $cmd = "cd $$p_config{pkgdir} && rm -fr mongo.comands " . ( join ' ', @pkgdirs ) . " 2> /dev/null";
        #        run_cmd( $cmd );
        my $cmd = "cd $$p_config{pkgdir} && mkdir " . ( join ' ', @pkgdirs ) . " 2> /dev/null";
        run_cmd( $cmd, true );
    }

    ## copy files

    {
        my $cmd;
        
        for my $from ( keys %exists ) {
            my $d = $exists{$from};
            $cmd .=  " && \\\n" if $cmd;
            $cmd .= "cd $from && ls | xargs -P4 cp -t $$p_config{pkgdir}/$d";
        }
        log_append( "$cmd", '-' );
        log_flush();
        run_cmd( $cmd );
    }
}
