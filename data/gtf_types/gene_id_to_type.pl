#!/usr/bin/perl

use warnings;
use strict;

use Data::Printer;


print join("\t",qw(gene_id gene_biotype))."\n";
while (<>) {
    next if /^#/;
    chomp;
    my @f = split /\t/;
    next unless $f[2] eq 'gene';
    my %t = map {/(\S+)\s+"([^"]+)"/ && ($1,$2)} split /;/,$f[8];
    print join("\t",$t{gene_id},$t{gene_biotype})."\n";
}
