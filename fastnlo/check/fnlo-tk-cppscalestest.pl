#!/usr/bin/env perl
#
# Unit test to evaluate a fastNLO table for 7 scale factor combinations
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
foreach my $file ( "cppscalestest.log", "cppscalesdiff.log", "${tabfl}", "${tabgz}" ) {
    if ( -e $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy gzipped default tables from data/check storage
my $cmd = "cp -f ${src}/../data/check/${tabgz} .";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppscalestest: Copying test table ${tabgz} failed: $ret, aborted!\n";}
$cmd = "gunzip ${tabgz}";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppscalestest: Ungzipping test table ${tabgz} failed: $ret, aborted!\n";}

# Evaluate, keep only result lines
$cmd = "../src/fnlo-tk-cppread ${tabfl} _ 7 | grep \"Calculate my cross sections\" -B1 -A95 > cppscalestest.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppscalestest: Evaluating test table ${tabfl} failed: $ret, aborted!\n";}

# Determine difference to default evaluation output
$cmd = "diff ${src}/../data/check/${tab}_cppscales.log cppscalestest.log > cppscalesdiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-cppscalestest: Result comparison with 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "cppscalesdiff.log" ) {
    print "fnlo-tk-cppscalestest: Evaluation of test table differs from default:\n";
    $cmd = "cat cppscalesdiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    print "fnlo-tk-cppscalestest: Do you use LHAPDF version 6? Do you have the CT10nlo PDF set installed? Is it found?\n";
    die "fnlo-tk-cppscalestest: Table evaluation unit test failed, please fix!\n";
}

print "fnlo-tk-cppscalestest: Table evaluation unit test passed.\n";

exit 0;
