#!/usr/bin/perl
# go_gene_association.pl connects genes with go terms
# and is released under the terms of the GNU GPL version 3, or any
# later version, at your option. See the file README and COPYING for
# more information.
# Copyright 2014 by Don Armstrong <don@donarmstrong.com>.


use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

=head1 NAME

go_gene_association.pl - connects genes with go terms

=head1 SYNOPSIS

go_gene_association.pl [options]

 Options:
   --debug, -d debugging level (Default 0)
   --help, -h display this help
   --man, -m display manual

=head1 OPTIONS

=over

=item B<--debug, -d>

Debug verbosity. (Default 0)

=item B<--help, -h>

Display brief usage information.

=item B<--man, -m>

Display this manual.

=back

=head1 EXAMPLES

go_gene_association.pl

=cut

use IO::File;
use Term::ProgressBar::IO;

use vars qw($DEBUG);

my %options = (debug           => 0,
               help            => 0,
               man             => 0,
              );

GetOptions(\%options,
           'debug|d+','help|h|?','man|m');

pod2usage() if $options{help};
pod2usage({verbose=>2}) if $options{man};

$DEBUG = $options{debug};

my @USAGE_ERRORS;
if (@ARGV != 2) {
    push @USAGE_ERRORS,"You must provide a go-basic.obo and a gene_association.goa file";
}

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;

my $gb = IO::File->new($ARGV[0]) or
    die "Unable to open $ARGV[0] for reading: $!";

my $p = Term::ProgressBar::IO->new($gb);

my %go_terms;
my $in_term;
my $current_term;
while (<$gb>) {
    chomp;
    if (/^\[Term\]/) {
        $in_term = 1;
        next;
    }
    next unless $in_term;
    my ($key,$value) = $_ =~ m/([^:]+)\:\s+(.+?)\s*$/;
    next unless defined $key;
    if ($key eq 'id') {
        $current_term = $value;
        $go_terms{$current_term} = {};
    }
    if ($key eq 'alt_id') {
        $go_terms{$value} = $go_terms{$current_term} if
            not exists $go_terms{$value};
    }
    if (not exists $go_terms{$current_term}{$key}) {
        $go_terms{$current_term}{$key} =
            $value;
    } else {
        if (ref($go_terms{$current_term}{$key})) {
            push @{$go_terms{$current_term}{$key}},
                $value;
        } else {
            $go_terms{$current_term}{$key} =
                [$go_terms{$current_term}{$key},
                 $value,
                ];
        }
    }
    $p->update();
}

my $genes = IO::File->new($ARGV[1]) or
    die "Unable to open $ARGV[1] for reading: $!";

while (<$genes>) {
    s/(GO:\d+)/$1\t$go_terms{$1}{name}\t$go_terms{$1}{namespace}/;
    print $_;
}


__END__
