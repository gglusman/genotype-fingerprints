#!/bin/env perl
use strict;
my $version = '180618';
####
#
# This software computes genotype fingerprints for multiple individuals as represented in a PLINK file.
# The method is described in Glusman et al., https://doi.org/10.1101/208025
# 
# Copyright 2017 by Gustavo Glusman, Institute for Systems Biology, Seattle, WA, USA.
# It is provided by the Institute for Systems Biology as open source software,
# free for non-commercial use.
# See the accompanying LICENSE file for information about the governing license.
#
####
#
# The expected input format is Plink bim/fam/bed.
# The first parameter is an 'id' for the job; output files will use this as base. You can include in it a path to where you want the output files to be located.
# The second parameter is the input file (i.e., the base of the Plink files, without the .bed).
# The third (optional) parameter is the fingerprint size. Multiple sizes can be specified, comma-delimited. Defaults to including several sizes.
#
####
#
# Examples of usage:
#   plinkGF.pl outDir plinkfiles/myDataSet
#   plinkGF.pl arbitrary_path/outDir genomes/plink 1000,2000
#
####

my($outdir, $file, $L) = @ARGV;
my @Ls = (500, 1000, 5000);
@Ls = split /,/, $L if $L;	
#$file ||= $id;
die unless $file;

# Sanitize inputs.
die if $outdir =~ /[\s\;]/;
die if $file =~ /[\s\;]/;

# Preparation.
my $rawFileExt = 'out';
my $normFileExt = 'outn';
my $chunkSize = 1000;
my $plink = `which plink`;
chomp $plink;
die "Cannot find plink\n" unless $plink;
my $tmp = "tmp$$";

# Read 1000g frequencies.
my %tgf;
open TGF, "gunzip -c genotype-fingerprints/data/1000g.freq.gz |";
while (<TGF>) {
	chomp;
	my($rsid, %f) = split /\t/;
	$tgf{$rsid} = \%f;
}
close TGF;

# Read list of rsids from .bim file.
my %rsids;
open BIM, "cut -f1,2 $file.bim |";
while (<BIM>) {
	chomp;
	my($chrom, $rsid) = split /\t/;
	push @{$rsids{$chrom}}, $rsid if $chrom>=1 && $chrom<=22;
}
close BIM;

my(%count, %alleleCount);

# Encode genotypes from --list transformed .bed file, per chromosome, in chunks.
foreach my $chrom (sort {$a<=>$b} keys %rsids) {
	my $rsids = $rsids{$chrom};
	my $n = scalar @{$rsids};
	for (my $chunk=0;$chunk<$n;$chunk+=$chunkSize) {
		my $startRsid = $rsids->[$chunk];
		my $endRsid = $rsids->[$chunk+$chunkSize-1] || $rsids->[-1];
		my $now = `date`;
		chomp $now;
		print join("\t", $now, $chrom, $chunk, $startRsid, $endRsid), "\n";
		
		`nice $plink --bfile $file --list --noweb --out $tmp --from $startRsid --to $endRsid`;
		open LST, "$tmp.list";
		while (<LST>) {
			chomp;
			my($chrom, $rsid, $gt, @ids) = split;
			next unless $rsid =~ /^rs\d+/;
			next unless defined $tgf{$rsid};
			$gt = uc $gt;
			next unless $gt =~ /^[ACGT][ACGT]$/;
			my $nid = substr($rsid,2);
			for (my $i=0;$i<$#ids;$i+=2) {
				my $fam = $ids[$i];
				my $id = $ids[$i+1];
				$id = join("-", $fam, $id) unless $fam eq $id;
				my %c;
				foreach my $i (0,1) { $c{substr($gt, $i, 1)}++ }
				foreach my $L (@Ls) {
					my $bin = $nid % $L;
					foreach my $al (keys %c) {
						$count{$id}{$L}{$al}[$bin] += $c{$al};
					}
					foreach my $tgal (keys %{$tgf{$rsid}}) {
						$count{$id}{$L}{$tgal}[$bin] -= $tgf{$rsid}{$tgal};
					}
				}
				$alleleCount{$id} += 2;
			}
		}
		close LST;
	}
}
unlink "$tmp.list";
unlink "$tmp.log";

# Output
mkdir $outdir;
foreach my $id (keys %count) {
	my @keys = sort keys %{$count{$id}{$Ls[0]}};
	my @headers = (
		['#software-version', $version],
		['#source', $file],
		['#alleleCount', $alleleCount{$id}],
		['#vectorLengths', join("\t", @Ls)],
		['#created', `date`],
		);
	my $header = join("\n", map {join("\t", @{$_})} @headers);
	
	# Output main fingerprint table.
	open OUTF, ">$outdir/$id.$rawFileExt";
	print OUTF $header;
	foreach my $vl (@Ls) {
		foreach my $key (@keys) {
			print OUTF join("\t", $vl, $key, @{$count{$id}{$vl}{$key}}), "\n";
		}
	}
	close OUTF;
	
	# Normalize the fingerprints.
	foreach my $vl (@Ls) {
		# Normalize fingerprint per reduced distance.
		foreach my $col (0..$vl-1) {
			my @v = ();
			push @v, $count{$id}{$vl}{$_}[$col] foreach @keys;
			my($avg, $std) = avgstd(\@v);
			$std ||= 1;
			$count{$id}{$vl}{$_}[$col] = ($count{$id}{$vl}{$_}[$col]-$avg)/$std foreach @keys;
		}
		# Normalize fingerprint per SNV pair key.
		foreach my $sig (@keys) {
			my($avg, $std) = avgstd($count{$id}{$vl}{$sig});
			$std ||= 1;
			$_ = ($_-$avg)/$std foreach @{$count{$id}{$vl}{$sig}};
		}
	}
	
	# Output normalized fingerprint table.
	open OUTF, ">$outdir/$id.$normFileExt";
	print OUTF $header;
	foreach my $vl (@Ls) {
		foreach my $key (@keys) {
			print OUTF join("\t", $vl, $key, map {sprintf("%.4f", $_)} @{$count{$id}{$vl}{$key}}), "\n";
		}
	}
	close OUTF;
	
	# Compress output files.
	`gzip -f $outdir/$id.$rawFileExt; gzip -f $outdir/$id.$normFileExt`;
}

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
