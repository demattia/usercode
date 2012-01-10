#!/usr/bin/perl -w

# This script parses the efficiency files from Kristian
# and produces three parts of output:
#  * arrays of efficiencies vs. mass suitable for pasting
#    into makeEfficiencyFits.C or suchlike
#  * arrays of efficiencies vs. lifetime suitable for pasting
#    into a script to fit vs. lifetime
#  * values of efficiency factor suitable for pasting into
#    runLifetimeJobs.pl
#
# "efficiency factor" is simply the combination of the two efficiencies
# and BR: effFactor = 2*epsilon_1*BR*(1-BR) + 2*epsilon_2*BR^2

use strict;

my $indir = "limits";

# for everything all at once
my %sets = ("e1" => "dielectron1",
	    "e2" => "dielectron2",
	    "mu1" => "dimuon1",
	    "mu2" => "dimuon2");

my @ctaufacts = ("0.033", "0.050", "0.100", "0.333", "0.500", "2.000", "3.000", "10.000", "20.000", "30.000");

foreach my $set (keys %sets) {
    foreach my $ctaufact (@ctaufacts) {
	$sets{$set.$ctaufact} = $sets{$set}."_ctaufact_".$ctaufact;
    }
}
    

# note: no pileup uncertainties for ctau reweighted samples,
# use the unweighted samples

my @uncertainties = ("statistical", "pileup");

my %masses = (200 => ["020", "050"],
	      400 => ["005", "020", "050", "150"],
	      1000 => ["020", "050", "150", "350"]);

my %nominal_lifetimes = ("1000_350" => 35,
			 "1000_150" => 10,
			 "1000_050" => 4,
			 "1000_020" => 1.5,
			 "400_150" => 40,
			 "400_050" => 8,
			 "400_020" => 4,
			 "200_050" => 20,
			 "200_020" => 7,
			 "120_050" => 50,
			 "120_020" => 13);

my %eff;
my %poserr2;
my %negerr2;
my %efferr;

# write out the xvalues first, since they're constant
foreach my $hmass (sort {$b <=> $a} keys %masses) {
    my @xmasses = @{$masses{$hmass}};
    my $npts = scalar(@xmasses);
    my @effvals;
    my @errvals;
    
    print "float eff_${hmass}_x[$npts] = {".
	join(", ", @xmasses)."};\n";
    print "float ex_${hmass}[$npts] = {". join(", ", ("0") x $npts)."};\n";
}
print "\n";

foreach my $set (sort keys %sets) {
    my $filetag = $sets{$set};

    my $infile = "$indir/efficiencies_${filetag}.txt";
    open INFILE, $infile or die "Couldn't open $infile: $!\n";

    while (my $line = <INFILE>) {
	chomp $line;
	my @fields = split(" ", $line);
	foreach my $hmass (keys %masses) {
	    my @xmasses = @{$masses{$hmass}};
	    foreach my $xmass (@xmasses) {
		
		my $pat = "${hmass}_${xmass}";
		if ($line =~ /$pat/) {
		    $eff{$set}{$hmass}{$xmass} = $fields[1];
		}
	    }
	}
    }
    close INFILE;

    foreach my $uncertainty (@uncertainties) {
	my $filename = "$indir/efficiencies_${filetag}_${uncertainty}_uncertainty.txt";

	# have to use the unweighted pileup uncertainty since weighted pileup uncertainties
	# don't exist
	$filename =~ s/ctaufact_\d+\.\d\d\d_// if ($uncertainty eq "pileup");

	open INFILE, $filename or die "Couldn't open $filename: $!\n";
	
	while (my $line = <INFILE>) {
	    chomp $line;
	my @fields = split(" ", $line);
	    foreach my $hmass (keys %masses) {
		my @xmasses = @{$masses{$hmass}};
		foreach my $xmass (@xmasses) {

		    my $pat = "${hmass}_${xmass}";
		    if ($line =~ /$pat/) {
			$poserr2{$set}{$hmass}{$xmass} += $fields[1]*$fields[1];
			$negerr2{$set}{$hmass}{$xmass} += $fields[2]*$fields[2];
		    }
		}
	    }
	}
	close INFILE;
    } # loop over uncertainties

    # now write it out
    foreach my $hmass (sort {$b <=> $a} keys %masses) {
	my @xmasses = @{$masses{$hmass}};
	my $npts = scalar(@xmasses);
	my @effvals;
	my @errvals;

	foreach my $xmass (@xmasses) {
	    push @effvals, $eff{$set}{$hmass}{$xmass};

	    my $poserr = sqrt($poserr2{$set}{$hmass}{$xmass});
	    my $negerr = sqrt($negerr2{$set}{$hmass}{$xmass});
	    my $symerr = sprintf("%.4f", ($poserr + $negerr)/2);
	    $efferr{$set}{$hmass}{$xmass} = $symerr;
	    push @errvals, $symerr;

	} # loop over x mass

	print "float eff${set}_${hmass}_y[$npts] = {".
	    join(", ", @effvals)."};\n";
	print "float eff${set}_${hmass}_ey[$npts] = {".
	    join(", ", @errvals)."};\n";
	print "\n";

    } # loop over h mass
} # loop over sets

# Finally produce the variables for the fits vs. lifetime.

#my %lifetimetags = (1 => "", 3 => "up", 1/3 => "down")
my %lifetimetags = (1 => "");
foreach my $ctaufact (@ctaufacts) {
    $lifetimetags{$ctaufact} = $ctaufact;
}
my $nlifetimes = scalar(keys %lifetimetags);

print "\n// Results by lifetime\n\n";

foreach my $set ("e1", "e2", "mu1", "mu2") {
	
    foreach my $hmass (keys %masses) {
	my @xmasses = @{$masses{$hmass}};
	foreach my $xmass (@xmasses) {
	    
	    my $pat = "${hmass}_${xmass}";
	    next if (!$nominal_lifetimes{$pat});
	    
	    my @lifetimes;
	    my @effvals;
	    my @errvals;
	    foreach my $lifetime (sort keys %lifetimetags) {
		my $modset = $set . $lifetimetags{$lifetime};

		push @lifetimes, ($lifetime * $nominal_lifetimes{$pat});
		push @effvals, $eff{$modset}{$hmass}{$xmass};
		push @errvals, $efferr{$modset}{$hmass}{$xmass};
		
	    } # lifetime loop

	    # print it out

	    print "float lifetimes_${set}_${pat}_y[$nlifetimes] = {".
		join(", ", @lifetimes)."};\n";
	    print "float eff_lifetime${set}_${pat}_y[$nlifetimes] = {".
		join(", ", @effvals)."};\n";
	    print "float eff_lifetime${set}_${pat}_ey[$nlifetimes] = {".
		join(", ", @errvals)."};\n\n";

	} # xmass loop
    } # hmass loop
} # top loop

# One final loop!
# This is to calculate the effective efficiencies for the known samples.
    
my $br = 0.01;

print "\n// Effective efficiency values\n\n";

foreach my $set ("e", "mu") {
    print "$set:\n";
    foreach my $hmass (keys %masses) {
	my @xmasses = @{$masses{$hmass}};
	foreach my $xmass (@xmasses) {
	    
	    my $pat = "${hmass}_${xmass}";
	    next if (!$nominal_lifetimes{$pat});

	    foreach my $lifetime (sort keys %lifetimetags) {
		my $eff1 = $eff{$set."1".$lifetimetags{$lifetime}}{$hmass}{$xmass};
		my $eff2 = $eff{$set."2".$lifetimetags{$lifetime}}{$hmass}{$xmass};
		my $efferr1 = $efferr{$set."1".$lifetimetags{$lifetime}}{$hmass}{$xmass};
		my $efferr2 = $efferr{$set."2".$lifetimetags{$lifetime}}{$hmass}{$xmass};

		my $efffactor = 2*$eff1*$br*(1-$br) + 2*$eff2*$br*$br;
		my $effterm1 = 2*$efferr1*$br*(1-$br);
		my $effterm2 = 2*$efferr2*$br*$br;
		my $efferr = sqrt($effterm1*$effterm1 + $effterm2*$effterm2);

		my $trackerr = 0.2*$efffactor;

		my $toterr = sqrt($efferr*$efferr + $trackerr*$trackerr);
		
		# print "percent error = ", sprintf("%.2f", ($toterr/$efffactor)), "\n";
		my $lt = $lifetime * $nominal_lifetimes{$pat};

		print "\"$pat"."_$lt\" => \[$efffactor, $toterr\],\n";
	    } # lifetime loop
	} # x mass loop
    } # h mass loop
} # e/mu loop

# Repeat for "Jeremy file"
foreach my $set ("e", "mu") {
    foreach my $hmass (keys %masses) {
	my @xmasses = @{$masses{$hmass}};
	foreach my $xmass (@xmasses) {
	    
	    my $pat = "${hmass}_${xmass}";
	    next if (!$nominal_lifetimes{$pat});

	    foreach my $lifetime (sort keys %lifetimetags) {
		my $eff1 = $eff{$set."1".$lifetimetags{$lifetime}}{$hmass}{$xmass};
		my $eff2 = $eff{$set."2".$lifetimetags{$lifetime}}{$hmass}{$xmass};

		my $efffactor = 2*$eff1*$br*(1-$br) + 2*$eff2*$br*$br;
		my $lt = sprintf("%.1f", $lifetime * $nominal_lifetimes{$pat});

		my $modset = $set;
		$modset = "el" if ($set eq "e");

		my $modxmass = $xmass;
		$modxmass =~ s/^0//;

		print $modset . "_MH${hmass}_MX${modxmass}_ct$lt => $efffactor\n";
	    } # lifetime loop
	} # x mass loop
    } # h mass loop
} # e/mu loop
