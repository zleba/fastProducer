#!/usr/bin/env perl
#
# Unit test to catenate two fastNLO tables
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
my $tabl = "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_stat_";
my $tabn = "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_stat_";
my @tabs = ("InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_stat_0000", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_stat_0001", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-born-2jet_stat_0002", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_stat_0100", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_stat_0101", "InclusiveNJets_fnr0001midpHT_I723509_v23_fix-hhc-nlo-2jet_stat_0102" );
my @tabd = ("NJets_born-2jet_stat_0000.tab.gz", "NJets_born-2jet_stat_0001.tab.gz", "NJets_born-2jet_stat_0002.tab.gz", "NJets_nlo-2jet_stat_0100.tab.gz", "NJets_nlo-2jet_stat_0101.tab.gz", "NJets_nlo-2jet_stat_0102.tab.gz");
my @tabfls;
my @tabgzs;
foreach my $tab ( @tabs ) {
    push @tabfls, ${tab}.".tab";
    push @tabgzs, ${tab}.".tab".".gz";
}
foreach my $file ( "statlotest.log", "statnlotest.log", "statlodiff.log", "statnlodiff.log", @tabfls, @tabgzs ) {
    if ( -e $file ) {
        system("rm -f $file");
    }
}

# Prepare test setup
# Copy gzipped default tables from data/check storage
for ( my $i=0; $i <= $#tabd ; $i++ ) {
    my $cmd = "cp -f ${src}/../data/check/$tabd[$i] $tabgzs[$i]";
    print "Executing command: $cmd\n";
    my $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-stattest: Copying test table $tabd[$i] failed: $ret, aborted!\n";}
    $cmd = "gunzip $tabgzs[$i]";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    if ( $ret ) {die "fnlo-tk-stattest: Ungzipping test table $tabgzs[$i] failed: $ret, aborted!\n";}
}

# Statistical evaluation
my $cmd = "../src/fnlo-tk-statunc $tabl | grep \"Special info\" -B1 -A34 > statlotest.log";
print "Executing command: $cmd\n";
my $ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-stattest: Statistical evaluation of LO test tables ${tabl} failed: $ret, aborted!\n";}
$cmd = "../src/fnlo-tk-statunc $tabn | grep \"Special info\" -B1 -A34 > statnlotest.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-stattest: Statistical evaluation of NLO test tables ${tabn} failed: $ret, aborted!\n";}

# Determine difference to default statistical uncertainties
$cmd = "diff ${src}/../data/check/NJets_born-2jet_stat.log statlotest.log > statlodiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-stattest: Result comparison with LO 'diff' failed: $ret, aborted!\n";}
$cmd = "diff ${src}/../data/check/NJets_nlo-2jet_stat.log statnlotest.log > statnlodiff.log";
print "Executing command: $cmd\n";
$ret = system("$cmd");
if ( $ret ) {die "fnlo-tk-stattest: Result comparison with NLO 'diff' failed: $ret, aborted!\n";}

# The diff.log must be empty
if ( ! -z "statlodiff.log" ) {
    print "fnlo-tk-stattest: Statistical evaluation of LO test tables differs from default:\n";
    $cmd = "cat statlodiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-tk-stattest: Statistical evaluation unit test failed, please fix!\n";
}
if ( ! -z "statnlodiff.log" ) {
    print "fnlo-tk-stattest: Statistical evaluation of NLO test tables differs from default:\n";
    $cmd = "cat statnlodiff.log";
    print "Executing command: $cmd\n";
    $ret = system("$cmd");
    die "fnlo-tk-stattest: Statistical evaluation unit test failed, please fix!\n";
}

print "fnlo-tk-stattest: Statistical evaluation unit test passed.\n";

exit 0;
