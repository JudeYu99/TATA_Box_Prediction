#!usr/bin/perl -w
use strict;

open(IN,"$ARGV[0]")||die "not open the file\n";
open(OUT,">$ARGV[1]")|| die "not make the file\n";

### Add positive or negative label to each predicted TATA box. ###

# awk -F'\t' '{if($3=="gene") print $1,$4,$5,$7}' Sc_YJM993.gff > Sc_gene.gff
# perl PositiveNegativeLabel.pl TATA_distance.gff Dis50.gff


while (<IN>) {
	my @lines = split("\t");
	chomp($lines[-1]);
	if ($lines[-1] <= 50) {
		push (@lines, "1");
	} else {
		push (@lines, "0");
	}

	my $res = join("\t", @lines);

	#print "@lines\n";
	print OUT "$res\n";
}


close IN;
close OUT;

