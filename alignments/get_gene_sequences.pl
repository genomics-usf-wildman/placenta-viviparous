#!/usr/bin/perl

use warnings;
use strict;

use Getopt::Long;

## this might need to be changed on someone else's machine
use lib "/home/don/projects/ensembl/ensembl/modules";
use Bio::EnsEMBL::Registry;
use HTTP::Tiny;
use JSON;

my %options;
$options{flank} = 1e5;
GetOptions(\%options,
           'promoter_only|promoter-only!',
           'flank=i',
           'debug|d+','help|h|?','man|m');



our $http = HTTP::Tiny->new();

Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org',
                                              -user => 'anonymous'
                                             );


our %species_adaptors;


while (<>) {
    chomp;
    my ($id,$type) = $_ =~ /((ENS[A-Z]+)\d+)/;
    if (not defined $id) {
        print STDERR "Unknown id format $id";
        next;
    }
    my $orig = $_;
    $orig =~ s/^>//;
    # figure out spe
    my $adaptor = get_adaptor($type,$id);
    if (not defined $adaptor) {
        print STDERR "Unable to get adaptor for id:$id type:$type\n";
        next;
    }
    my $gene = $adaptor->{gene}->fetch_by_stable_id($id);
    # brilliant job, ensemble. Use camel case sometimes, not other
    # times.
    my $slice = $adaptor->{slice}->fetch_by_Feature($gene);
    if ($options{flank} > 0) {
        if ($options{promoter_only}) {
            if ($gene->strand() == -1) {
                # remove all but flank bases of the sequence, and add on the flank to the promotor
                $slice = $slice->expand(($slice->start()-$slice->end()+
                                         $options{flank}),
                                        $options{flank});
            } else {
                $slice = $slice->expand($options{flank},
                                        ($slice->start()-$slice->end()+
                                         $options{flank}));
            }
        } else {
            $slice = $slice->expand($options{flank},$options{flank});
        }
    }
    my $seq = $slice->seq();
    $seq =~ s/(.{100})/$1\n/g;
    $seq =~ s/\n+$//s;
    print ">$id|".$gene->external_name().'|'.$gene->species()."\n";
    print STDERR "got id:$id\n";
    print $seq;
    print "\n";
}

sub get_adaptor {
    my ($species,$id) = @_;
    if (not exists $species_adaptors{$species}) {
        my $ensembl_species;
        my $response = $http->get('http://rest.ensembl.org/lookup/id/'.$id,
                                 {headers => {'Content-type'=>'application/json'}
                                 }
                                 );

        return undef unless $response->{success};
        my $data = decode_json($response->{content});
        $ensembl_species = $data->{species};
        $species_adaptors{$species} =
           { species => $ensembl_species,
             slice => Bio::EnsEMBL::Registry->get_adaptor($ensembl_species, "core",
                                                          "Slice" ),
             gene  => Bio::EnsEMBL::Registry->get_adaptor($ensembl_species, "core",
                                                          "gene" ),
           };
        if (not defined $species_adaptors{$species}{slice} or
            not defined $species_adaptors{$species}{gene}) {
            return undef;
        }
    }
    return $species_adaptors{$species};
}
