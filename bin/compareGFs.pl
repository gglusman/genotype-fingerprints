#!/bin/env perl
use strict;
my $version = '171019';
####
#
# This software compares two genotype fingerprints.
# The method is described in Glusman et al., https://doi.org/10.1101/208025
# 
# Copyright 2017 by Gustavo Glusman, Institute for Systems Biology, Seattle, WA, USA.
# It is provided by the Institute for Systems Biology as open source software,
# free for non-commercial use.
# See the accompanying LICENSE file for information about the governing license.
#
####
#
# The first and second parameters are the files holding the genotype fingerprints to be compared.
# The third (optional) parameter is the fingerprint size to be used. Multiple sizes can be specified, comma-delimited. Defaults to comparing all possible fingerprint sizes, depending on what sizes are available for each of the genomes.
#
# The output displays the number of SNPs present in each of the two genotypes ('q_SNPs' and 't_SNPs' for 'query' and 'target'), the similarity between their binary fingerprints, and their Spearman correlation for each fingerprint size considered.
#
####
#
# Examples of usage:
#   compareGFs.pl fingerprints/genotype1.outn.gz fingerprints/genotype2.outn.gz
#   compareGFs.pl fingerprints/genotype1.outn.gz fingerprints/genotype2.outn.gz 1000,2000
#
####

my($f0, $f1, $L) = @ARGV;
my %todo;
$todo{$_}++ foreach (split /,/, $L);

my($dmf0, $snps0) = readGF($f0);
my($dmf1, $snps1) = readGF($f1);
my $folded0 = foldGF($dmf0);
my $folded1 = foldGF($dmf1);

my @h = qw/q_SNPs t_SNPs/;
my @f = ($snps0 || 'NA', $snps1 || 'NA');

foreach my $vl (sort {$a<=>$b} keys %$dmf0) {
	next unless defined $dmf1->{$vl};
	next if %todo && !$todo{$vl};
	push @h, "L=$vl";
	#my $pearson = correlation($folded0->{$L}, $folded1->{$L});
	my $spearman = correlation(ranks($folded0->{$vl}), ranks($folded1->{$vl}));
	push @f, sprintf("%.4f", $spearman);
}

print join("\t", @h), "\n";
print join("\t", @f), "\n";


###
sub readGF {
	my($file) = @_;
	my(%v, $alleles);
	
	if ($file =~ /\.gz$/) {
		open F, "gunzip -c $file |";
	} else {
		open F, $file;
	}
	while (<F>) {
		chomp;
		my($vl, $key, @v) = split /\t/;
		if ($vl =~ /^#/) {
			$alleles = $key if $vl eq '#alleleCount';
			next;
		}
		next if %todo && !$todo{$vl};
		$v{$vl}{$key} = \@v;
	}
	close F;
	return \%v, $alleles/2;
}

sub foldGF {
	my($v) = @_;
	my %f;
	
	foreach my $vl (keys %$v) {
		push @{$f{$vl}}, @{$v->{$vl}{$_}} foreach sort keys %{$v->{$vl}};
	}
	return \%f;
}

sub correlation {
	my($set1ref, $set2ref, $avg1, $std1, $avg2, $std2) = @_;
	my($n, $i, $corr);

	$n = scalar @{$set1ref};
	return unless $n && ($n == scalar @{$set2ref});
	($avg1, $std1) = avgstd($set1ref) unless $std1;
	($avg2, $std2) = avgstd($set2ref) unless $std2;
	foreach $i (0..$n-1) {
		$corr += ($$set1ref[$i]-$avg1)*($$set2ref[$i]-$avg2);
	}
	return $corr/($n-1)/$std1/$std2;
}

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

sub ranks {
	my($v) = @_;
	
	my @ranks;
	my $i;
	foreach (sort {$v->[$a] <=> $v->[$b]} (0..$#$v)) { $ranks[$_] = $i++ }
	return \@ranks;
}

sub compareBinaries {
	my($b0, $b1) = @_;
	return unless $b0 && $b1;
	my $c;
	my @v0 = split //, $b0;
	my @v1 = split //, $b1;
	foreach my $i (0..$#v0) {
		$c++ if $v0[$i] == $v1[$i];
	}
	return ($c/scalar @v0)**2;
}

