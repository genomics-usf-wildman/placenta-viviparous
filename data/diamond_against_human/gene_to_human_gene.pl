#!/usr/bin/perl

use warnings;
use strict;

use IO::File;
use Data::Dumper;

my ($prot_to_gene,@alignment_files) = @ARGV;

sub open_file {
    my $fh = IO::File->new($_[0],'r') or
        die "Unable to open $_[0] for reading: $!";
    return $fh;
}

my $p2g_fh = open_file($prot_to_gene);

my %p2g;
while (<$p2g_fh>) {
    next unless /^ENS/;
    chomp;
    my ($prot,$gene,$name) = split /\t/;
    $p2g{$prot} = "$gene\t$name";
}

sub print_match {
    my ($gene,$prot,$human_prot) = @_;
    print join("\t",$gene,$prot,$p2g{$human_prot},$human_prot)."\n";
}

print join("\t","gene_id","protein_id","human_gene_id","human_gene_symbol","human_protein")."\n";
for my $alignment_fh (map {open_file($_)} @alignment_files) {
    my $previous_id_cont = '';
    my $previous_gene_id;
    my $previous_prot_id;
    my $previous_human_protein;
    my $previous_bit_score = 0;
    while (<$alignment_fh>) {
        chomp;
        my ($id_cont,$human_protein,$identity_per_len,
            $len,$mismatches,$gaps,$query_begin,$query_end,
            $subject_begin,$subject_end,$evalue,$bitscore) = split /\t/;
        ## we only want the top scoring match
        next if $id_cont eq $previous_id_cont;
        $previous_id_cont = $id_cont;
        my ($prot_id,$gene_id) = split /_/,$id_cont;
        if (defined $previous_gene_id and
            $gene_id ne $previous_gene_id) {
            print_match($previous_gene_id,$previous_prot_id,
                        $previous_human_protein);
            $previous_gene_id=undef;
            $previous_prot_id=undef;
            $previous_human_protein=undef;
            $previous_bit_score = 0;
        }
        if (not defined $previous_gene_id or
            $previous_bit_score < $bitscore
           ) {
            $previous_gene_id = $gene_id;
            $previous_prot_id = $prot_id;
            $previous_human_protein = $human_protein;
            $previous_bit_score = $bitscore;
        }
    }
    if (defined $previous_gene_id) {
        print_match($previous_gene_id,$previous_prot_id,
                    $previous_human_protein);
    }
}
