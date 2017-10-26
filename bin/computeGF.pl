#!/bin/env perl
use strict;
my $version = '171019';
####
#
# This software computes a genotype fingerprint for a single personal genome.
# The method is described in Glusman et al., https://doi.org/10.1101/130807
# 
# Copyright 2017 by Gustavo Glusman, Institute for Systems Biology, Seattle, WA, USA.
# It is provided by the Institute for Systems Biology as open source software,
# free for non-commercial use.
# See the accompanying LICENSE file for information about the governing license.
#
####
#
# The expected input format is 23andMe's genotype format: rsid chrom pos genotype.
# The first parameter is an 'id' for the job; output files will use this as base. You can include in it a path to where you want the output files to be located.
# The second parameter is the input file (the genome as VCF, RCF, var or masterVar).
# The third (optional) parameter is the format of the input file: 'vcf', 'rcf', 'var' or 'masterVar'. Defaults to 'vcf'.
# The fourth (optional) parameter is the fingerprint size. Multiple sizes can be specified, comma-delimited. Defaults to including several sizes.
# The fifth (optional) parameter is the distance between consecutive SNVs that are considered 'too close'. Default is 20.
# The sixth (optional) parameter is a bed file specifying regions of interest to be included in the analysis. For example, one could specify the definition of exome segments to compute an exome-compatible fingerprint from whole-genome data. This is available only for VCF and RCF input.
#
####
#
# Examples of usage:
#   computeGF.pl myGenome genotypes/myGenome.gz
#   computeGF.pl fingerprints/anotherGenome genotypes/aGenome.gz 1000,2000
#
####

my($id, $file, $L) = @ARGV;
my @Ls = (500, 1000, 5000);
@Ls = split /,/, $L if $L;	
$file ||= $id;

# Sanitize inputs.
die if $file =~ /[\s\;]/;

# Preparation.
my $rawFileExt = 'out';
my $normFileExt = 'outn';

my(%count, $alleleCount);
my $cat = 'cat';
if ($file =~ /\.gz$/) {
	$cat = 'gunzip -c';
} elsif ($file =~ /\.bz2$/) {
	$cat = 'bzcat';
}

# Read 1000g frequencies.
my %tgf;
open TGF, "gunzip -c genotype-fingerprints/data/1000g.freq.gz |";
while (<TGF>) {
	chomp;
	my($rsid, %f) = split /\t/;
	$tgf{$rsid} = \%f;
}
close TGF;

# Process input file.
open INF, "$cat $file |";
while (<INF>) {
	chomp;
	my($rsid, $chrom, $pos, $gt) = split /\t/;
	
	# Focus the analysis on autosomes only, excluding sex chromosomes, mitochondrial chromosome, alternative haplotypes, etc.
	next unless $chrom =~ /^(chr)?\d+$/;
	
	# Pay attention only to SNVs with rsids.
	$gt =~ s/\W//g;
	$gt = uc $gt;
	next unless $gt =~ /^[ACGT][ACGT]$/;
	next unless defined $tgf{$rsid};
	next unless $rsid =~ /^rs\d+/;
	my $nid = substr($rsid,2);
	
	# Add counts
	my %c;
	foreach my $i (0,1) { $c{substr($gt, $i, 1)}++ }
	foreach my $L (@Ls) {
		my $bin = $nid % $L;
		foreach my $al (keys %c) {
			$count{$L}{$al}[$bin] += $c{$al};
		}
		foreach my $tgal (keys %{$tgf{$rsid}}) {
			$count{$L}{$tgal}[$bin] -= $tgf{$rsid}{$tgal};
		}
	}
	$alleleCount += 2;
		
}
close INF;

my @keys = sort keys %{$count{$Ls[0]}};

my @headers = (
	['#software-version', $version],
	['#source', $file],
	['#alleleCount', $alleleCount],
	['#vectorLengths', join("\t", @Ls)],
	['#created', `date`],
	);
my $header = join("\n", map {join("\t", @{$_})} @headers);

# Output main fingerprint table.
open OUTF, ">$id.$rawFileExt";
print OUTF $header;
foreach my $vl (@Ls) {
	foreach my $key (@keys) {
		print OUTF join("\t", $vl, $key, @{$count{$vl}{$key}}), "\n";
	}
}
close OUTF;

# Normalize the fingerprints.
foreach my $vl (@Ls) {
	# Normalize fingerprint per reduced distance.
	foreach my $col (0..$vl-1) {
		my @v = ();
		push @v, $count{$vl}{$_}[$col] foreach @keys;
		my($avg, $std) = avgstd(\@v);
		$std ||= 1;
		$count{$vl}{$_}[$col] = ($count{$vl}{$_}[$col]-$avg)/$std foreach @keys;
	}
	# Normalize fingerprint per SNV pair key.
	foreach my $sig (@keys) {
		my($avg, $std) = avgstd($count{$vl}{$sig});
		$std ||= 1;
		$_ = ($_-$avg)/$std foreach @{$count{$vl}{$sig}};
	}
}

# Output normalized fingerprint table.
open OUTF, ">$id.$normFileExt";
print OUTF $header;
foreach my $vl (@Ls) {
	foreach my $key (@keys) {
		print OUTF join("\t", $vl, $key, map {sprintf("%.4f", $_)} @{$count{$vl}{$key}}), "\n";
	}
}
close OUTF;

# Compress output files.
`gzip -f $id.$rawFileExt; gzip -f $id.$normFileExt`;


###
sub avgstd {
	my($values) = @_;
	my($sum, $devsqsum);

	my $n = scalar @$values;
	return unless $n>1;
	$sum += $_ foreach @$values;
	my $avg = $sum / $n;
	$devsqsum += ($_-$avg)**2 foreach @$values;
	my $std = sqrt($devsqsum/($n-1));
	return $avg, $std;
}
