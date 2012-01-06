#!/usr/bin/perl -w

use strict;

my $filename = "CountingLimits/cls_centralHighLow_newpoints.txt";

open INFILE, $filename or die "Couldn't open $filename: $!\n";

my %results;

while (my $line = <INFILE>) {
    chomp $line;

    $line =~ /^(..)_MH(\d+)_MX(\d+)_ct(\d+.\d).*upper limit: ([0-9.]+)/;

    my $type = $1;
    my $key = $2."_".$3."_".$4;
    my $ul = $5;

    my $mx = $3;
    my $index = 1;
    $index = 2 if ($line =~ /high/);
    $index = 3 if ($line =~ /low/);

    $results{$type}{$key}[0] = $mx;
    $results{$type}{$key}[$index] = $ul;
}

foreach my $type (keys %results) {
    foreach my $key (keys %{ $results{$type}}) {
	my $prefix;
	if ($type eq "el") { $prefix = "ElectronsCLs/"; }
	elsif ($type eq "mu") { $prefix = "MuonsCLs/"; }
	else { die "Unrecognized type $type\n"; }

	my $outfile = $prefix.$key.".txt";

	open OUTFILE, ">$outfile" or die "Couldn't open $outfile: $!\n";

	print OUTFILE join(" ", @ {$results{$type}{$key}}), "\n";

	close OUTFILE;
    }
}
	
