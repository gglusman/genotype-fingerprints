#!/bin/env perl
use strict;
my $version = '171019';
use FindBin qw($Bin);
####
#
# This software compares two sets of genotype fingerprints.
# The method is described in Glusman et al., https://doi.org/10.1101/208025
# 
# Copyright 2017 by Gustavo Glusman, Institute for Systems Biology, Seattle, WA, USA.
# It is provided by the Institute for Systems Biology as open source software,
# free for non-commercial use.
# See the accompanying LICENSE file for information about the governing license.
#
####
#
# The first parameters is the query set. It can include one or more fingerprints in serialized format. If a non-serialized fingerprint is used as query, it is automatically serialized using the fingerprint size of the target set.
# The second parameter is the target set: a serialized collection of fingerprints. If absent, all-against-all comparisons are performed in the query set.
#
####
#
# Examples of usage:
#   searchGFs.pl fingerprints/genotype1.outn.gz target.fp
#   searchGFs.pl query.fp target.fp
#   searchGFs.pl self-query.fp
#
####

my $fpc = "$Bin/fpc";
die "Cannot find fpc (the search engine)\n" unless -s $fpc;
my($query, $target) = @ARGV;
unless ($query) {
	print "Usage: searchGFs.pl fingerprints/genome1.outn.gz target.fp\n";
	print "       searchGFs.pl query-set target-set\n";
	print "       searchGFs.pl query-set  (will perform all comparisons within the query data set)\n";
	exit;
}


my($qnames, $tnames, $L);

if (!$target) {
	#will use same as query
} elsif ($target =~ /(.+)\.fp/) {
	$tnames = readIds("$1.id");
} elsif (-e "$target.id") {
	$tnames = readIds("$target.id");
	$target = "$target.fp";
} else {
	die "Couldn't interpret $target as target\n";
}



if ($query =~ /(.+)\.fp/) {
	$qnames = readIds("$1.id");
} elsif (-e "$query.id") {
	$qnames = readIds("$query.id");
	$query = "$query.fp";
} else {
	`./serializeGFs.pl tmp$$ $L $query`;
	$qnames = readIds("tmp$$.id");
	unless (@$qnames) {
		unlink "tmp$$.fp";
		unlink "tmp$$.id";
		die "Couldn't interpret $query as query\n";
	}
	$query = "tmp$$.fp";
}

if ($target) {
	open FPC, "$fpc $query $target |";
} else {
	open FPC, "$fpc $query |";
	$tnames = $qnames;
}
while (<FPC>) {
	chomp;
	my($q, $t, $c) = split /\t/;
	print join("\t", $qnames->[$q-1], $tnames->[$t-1], $c), "\n";
}
close FPC;

unlink "tmp$$.fp";
unlink "tmp$$.id";



sub readIds {
	my($idfile) = @_;
	my @names;
	open ID, "$idfile";
	while (<ID>) {
		if (/^#/) {
			$L = $1 if /^#L\t(\d+)/;
			next;
		}
		chomp;
		my(undef, $id, $file) = split /\t/;
		push @names, $id;
	}
	close ID;
	return \@names;
}


