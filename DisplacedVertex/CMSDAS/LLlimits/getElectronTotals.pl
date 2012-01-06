#!/usr/bin/perl -w

use strict;

my $filename = "data/masses_backgroundMC_Electrons.txt";

open INFILE, $filename or die "Couldn't open $filename: $!\n";

my $totev = 0;
my $toterr2 = 0;
my $peak5ev = 0;
my $peak15ev = 0;
my $outsideev = 0;

while (my $line = <INFILE>) {
    chomp $line;

    my @fields = split(" ", $line);
    my $mass = $fields[0];
    my $weight = $fields[1];

    if ($mass > 15.0) {
	$totev += $weight;
	$toterr2 += $weight*$weight;

	if ($mass > 87.5 && $mass < 92.5) { $peak5ev += $weight; }
	if ($mass > 82.5 && $mass < 97.5) { $peak15ev += $weight; }
	else { $outsideev += $weight; }
	
    }
}

my $toterr = sqrt($toterr2);
print "Total: $totev +/- $toterr2\n";
print "Within 5 GeV of Z peak: $peak5ev\n";
print "Within 15 GeV of Z peak: $peak15ev\n";
print "Outside 15 GeV of Z peak: $outsideev\n";
