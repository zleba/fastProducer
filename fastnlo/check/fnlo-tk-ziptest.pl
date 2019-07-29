#!/usr/bin/env perl
#
# Unit test to read a fastNLO table
# Version:
#
# created by K. Rabbertz: 20.10.2016
#
#-----------------------------------------------------------------------
#
use Cwd;
use English;
use Getopt::Std;
use strict;
use warnings;

# Define source dir as required for VPATH builds like in make distcheck
my $src = ".";

# Remove potentially left-over temporary files
my $tab   = "fnl1014_I902309_2-2";
my $tabfl = ${tab}.".tab";
my $tabgz = ${tabfl}.".gz";
foreach my $file ( "ziptest.log", "zipdiff.log", "${tabfl}", "${tabgz}" ) {
    if ( -e $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy gzipped default tables from data/check storage
my $cmd = "cp -f ${src}/../data/check/${tabgz} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-ziptest: Copying test table ${tabgz} failed: $ret, aborted!\n";}

# Evaluate gzipped table directly, keep only result lines
$cmd = "../src/fnlo-tk-cppread ${tabgz} | grep \"Calculate my cross sections\" -B1 -A11 > ziptest.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-ziptest: Evaluating gzipped test table ${tabgz} failed: $ret, aborted!\n";}

# Determine difference to default evaluation output
$cmd = "diff ${src}/../data/check/${tab}_cpp.log ziptest.log > zipdiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-ziptest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "zipdiff.log" ) {
    print "fnlo-tk-ziptest: Evaluation of test table differs from default:\n";
    $cmd = "cat zipdiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    print "fnlo-tk-ziptest: Do you use LHAPDF version 6? Do you have the CT10nlo PDF set installed? Is it found?\n";
    die "fnlo-tk-ziptest: Table evaluation unit test failed, please fix!\n";
}

print "fnlo-tk-ziptest: Table evaluation unit test passed.\n";

exit 0;
