#!/usr/bin/perl

use warnings;
use strict;

use IO::Uncompress::Gunzip;
use IO::File;


## converts gtf files into a protein_id->gene_id mapping
my %protein_ids;
print "protein_id\tgene_id\tgene_short_name\n";
for my $file (@ARGV) {
    my $fh = IO::Uncompress::Gunzip->new($file);
    while (<$fh>) {
        chomp;
        my ($gene_id) = $_ =~ /gene_id\s+"([^"]+)"/ or next;
        my ($protein_id) = $_ =~ /protein_id\s+"([^"]+)"/ or next;
        my ($gene_name) = $_ =~ /gene_name\s+"([^"]+)"/;
        $gene_name //= '';
        next if exists $protein_ids{$protein_id};
        $protein_ids{$protein_id} =  1;
        print "$protein_id\t$gene_id\t$gene_name\n";
    }
}

