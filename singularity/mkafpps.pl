#!/usr/bin/perl

my $afppbase   = "/scratch/00451/tg457210/af/afpp2";
my $targetbase = "/scratch/00451/tg457210/af/w";

my $json = <<"__EOF";
{
    "debug"           : true

    ,"ln"              : "cp"
    ,"afversion"       : 2                                                                           # version of alpha fold to process, number or "all"
    ,"base"            : "$targetbase/w_000/afpp"                                                         # base directory for called scripts
    ,"afftp"           : "$targetbase/w_000/af_ftp"                                                       # alphafold download
    ,"afpdb"           : "$targetbase/w_000/af_pdb"                                                       # directory of all downloaded alphafold pdbs
    ,"fasta"           : "$targetbase/w_000/af_blast"                                                     # directory of all downloaded fasta
    ,"fasta_url"       : "https://www.uniprot.org/uniprot"                                           # url for fasta 
    ,"unifeat"         : "$targetbase/w_000/af_features"                                                  # directory of all downloaded uniprot features
    ,"unifeat_url"     : "https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&accession="  # REST api url prefix to download features
    ,"chains"          : "$targetbase/w_000/af_chains"                                                    # chain summary files
    ,"cpnotes"         : "$targetbase/w_000/af_cpnotes"                                                   # chain processing notes
    ,"cperrs"          : "$targetbase/w_000/af_cperrors"                                                  # chain processing errors
    ,"pdbnotes"        : "$targetbase/w_000/af_pdbnotes"                                                  # pdb processing notes
    ,"pdbstage1"       : "$targetbase/w_000/af_pdbstage1"                                                 # cut pdbs for processing
    ,"pdbstage2"       : "$targetbase/w_000/af_pdbstage2"                                                 # processed pdbs ready for computations, posting
    ,"pdb_tfrev"       : "$targetbase/w_000/af_pdb_tfrev"                                                 # processed pdbs for jsmol TF
    ,"somocli"         : "/ultrascan3/us_somo/bin64/us_saxs_cmds_t"                                  # directory of US SOMO cli 
    ,"sesca"           : "/somoaf/lib/SESCA/scripts/SESCA_main.py"                                   # directory of SESCA_main.py
    ,"sescadir"        : "$targetbase/w_000/af_cd"                                                        # directory of SESCA computed output
    ,"somoenv"         : "export LD_LIBRARY_PATH=/ultrascan3/us_somo/lib:/qwt-6.1.5/lib"    # us3_somo environment setup
    ,"somorun"         : "timeout 3m xvfb-run -a /ultrascan3/us_somo/bin64/us3_somo"               # us3_somo executable run
    ,"somordir"        : "$targetbase/w_000/afpp/ultrascan/somo"                                     # directory of somo results
    ,"tempdir"         : "$targetbase/w_000/af_tmp"                                                       # directory of temporary/scratch dir
    ,"cifdir"          : "$targetbase/w_000/af_cif"                                                       # directory of generated cif files
    ,"mmcifdir"        : "$targetbase/w_000/af_mmcif"                                                     # directory of generated mmcif files
    ,"csvdir"          : "$targetbase/w_000/af_csv"                                                       # directory of generated csv files
    ,"prdir"           : "$targetbase/w_000/af_pr"                                                        # directory of generated pr files
    ,"zipdir"          : "$targetbase/w_000/af_zip"                                                       # directory of generated zip files
    ,"txzdir"          : "$targetbase/w_000/af_txz"                                                       # directory of generated txz files
    ,"maxit"           : "env RCSBROOT=/somoaf/lib/maxit-v11.100-prod-src /somoaf/bin/maxit"         # path to run maxit
    ,"mongodb"         : "somo"                                                                      # mongo database name
    ,"mongocoll"       : "afd"                                                                       # mongo collection name
    ,"mongocmds"       : "$targetbase/w_000/af_mongocmds"                                                 # directory of mongo update commands
    ,"pkgdir"          : "/work/00451/tg457210/ls6/af_pkg"

    ,"stages" : [
        {
            "name"    : "collect requried files"
            ,"cmd"    : "__base__/collect.pl __id__"
            ,"active" : false
        }
        ,{
            "name"    : "build chain summary"
            ,"cmd"    : "__base__/chains.pl __id__"
            ,"active" : true
        }
        ,{
            "name"    : "produce pdb variants"
            ,"cmd"    : "__base__/pdbs.pl __id__"
            ,"active" : true
        }
        ,{
            "name"    : "finalize pdb variants"
            ,"cmd"    : "__base__/pdb2.pl __id__"
            ,"active" : true
        }
        ,{
            "name"    : "computations"
            ,"cmd"    : "__base__/compute.pl __id__"
            ,"active" : true
        }
        ,{
            "name"    : "package"
            ,"cmd"    : "__base__/package.pl __id__"
            ,"active" : false
        }
    ]
}
__EOF

;

$usrc = <<__EOF;
4.0
/usr/bin/firefox
tar_dummy
zip_dummy
/ultrascan3/us_somo/somo/doc
$targetbase/w_000/afpp/ultrascan/data
$targetbase/w_000/afpp/ultrascan/
$targetbase/w_000/afpp/ultrascan/archive
$targetbase/w_000/afpp/ultrascan/results
0
0.5
$targetbase/w_000/afpp/ultrascan/reports
/ultrascan3/us_somo
Helvetica
10
10
128
$targetbase/w_000/afpp/ultrascan/tmp
__EOF
;


$notes = "
usage: $0 ids

id is a number of the split job
makes an afpp directory with proper config file 
under: $targetbase

afpp will be copied from $afppbase
also makes somo usrc.conf

";



die $notes if !@ARGV;

while ( @ARGV ) {
    my $id = shift;
    my $uid = '0'x(3 - length($id)) . $id;
    
    my $ujson = $json;
    $ujson =~ s/w_000/w_$uid/mg;

    my $td = "$targetbase/w_$uid";
    die "$td not a directory\n" if !-d $td;

    ## mk afppdir

    die "$td/afpp directory exists, remove with:\nrm -fr $td/afpp\n" if -d "$td/afpp";

    my $cmd = "cp -r $afppbase $td/afpp";
    
    print "$cmd\n";
    print `$cmd`;
    die "$cmd failed $?\n" if $?;

    {
        my $f = "$td/afpp/config.json";
        open my $fh, ">$f";
        print $fh $ujson;
        close $fh;
        print "$f\n";
    }

    ## make somo
    {
        my $uusrc = $usrc;
        $uusrc =~ s/w_000/w_$uid/mg;

        my $f = "$td/afpp/ultrascan/etc/usrc.conf";
        open my $fh, ">$f";
        print $fh $uusrc;
        close $fh;

        print "$f\n";
    }
}
  
