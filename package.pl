#!/usr/bin/perl

## set to one for testing (skip CD & somo calcs)

$skipproc = 0;

use File::Basename;
$ppdir = dirname(__FILE__);
require "$ppdir/utility.pm";

$notes = "usage: $0 ids

ids can be uniprot ids and or file names with uniprot ids

builds script to upload these to the server

";

die $notes if !@ARGV;

use File::Temp qw(tempfile);

p_config_init();

warn "clean $$p_config{pkgdir} before starting!!!
--------
rm -fr $$p_config{pkgdir}/*
--------
";

while ( $fid = shift ) {
    my $id = uniprot_id( $fid );
    log_append( "$0: processing ascension $id", "-" );

    my @pdbs   = get_pdbs2( $id );

    for my $fpdb ( @pdbs ) {
        my $pdb_ver     = pdb_ver    ( $fpdb );
        my $pdb_frame   = pdb_frame  ( $fpdb );
        my $pdb_variant = pdb_variant( $fpdb );
        my $fpdbnoext   = $fpdb;
        $fpdbnoext      =~ s/\.pdb//;
        
        my $f           = "$$p_config{pdbstage2}/$fpdb";
        my $fpdbtf      = "$$p_config{pdb_tfrev}/$fpdb";
        my $fcd         = "$$p_config{sescadir}/${fpdbnoext}-sesca-cd.dat";
        my $fmmcif      = "$$p_config{mmcifdir}/${fpdbnoext}.cif";
        my $fcsv        = "$$p_config{csvdir}/${fpdbnoext}.csv";
        my $fpr         = "$$p_config{prdir}/${fpdbnoext}-pr.dat";
        my $fzip        = "$$p_config{zipdir}/${fpdbnoext}.zip";
        my $ftxz        = "$$p_config{txzdir}/${fpdbnoext}.txz";
        my $fmongo      = "$$p_config{mongocmds}/${id}-${pdb_frame}${pdb_variant}.mongo";
        
        my %exists = (
            $f              => "pdb"
            ,$fpdbtf        => "pdb_tfrev"
            ,$fcd           => "cd"
            ,$fmmcif        => "mmcif"
            ,$fcsv          => "csv"
            ,$fpr           => "pr"
            ,$fzip          => "zip"
            ,$ftxz          => "txz"
            ,$fmongo        => ">> $$p_config{pkgdir}/mongo.cmds"
            );

        my @pkgdirs;

        my $any_errors = 0;

        for my $check ( keys %exists ) {
            if ( !-e $check ) {
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

        ## link files

        {
            my $cmd;
            
            for my $f2link ( keys %exists ) {
                my $d = $exists{$f2link};
                $cmd .=  " && \\\n" if $cmd;
                if ( $d !~ /^>/ ) {
                    $cmd .= "ln -f $f2link $$p_config{pkgdir}/$d/";
                } else {
                    $cmd .= "cat $f2link $d";
                }
            }
            # log_append( "$cmd", '-' );
            run_cmd( $cmd );
        }
    }
}

    
log_flush();
if ( $errors ) {
    print STDERR $errors;
}

exit length( $errors );    
