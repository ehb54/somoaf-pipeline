#!/usr/bin/perl

use POSIX;

## sort of xargs to handle parallel jobs which doesn't seem to work with -Pn on TACC's compute notes

$notes = "
usage: $0 processes command-line file

quote the command-line if is the more than one token

";

use File::Temp qw(tempfile);

$p   = shift || die $notes;
$cmd = shift || die $notes;
$f   = shift || die $notes;

die "processes must be > 0\n"   if $p <= 0;
die "$f not found\n"            if !-e $f;
die "$f not readable\n"         if !-r $f;

open IN, $f;
@l = <IN>;
close IN;
grep chomp, @l;

die "$f empty\n" if !@l;

## split remaining args into 

my $jpp = ceil( @l / $p );

$tmpdir = "/scratch/00451/tg457210/af/workv2/af_tmp";

$cmds = "";

while ( @l ) {
    my ( $fh, $ft ) = tempfile( "$tmpdir/xargs.XXXXXX", UNLINK => 1 );
    my @args = splice @l, 0, $jpp;
    print $fh join "\n", @args;
    close $fh;
    my $rcmd = "cat $ft | xargs $cmd &";
    $cmds .= "$rcmd\n";
}

$cmds .= "wait\n";

print $cmds;
print `$cmds`;
    
