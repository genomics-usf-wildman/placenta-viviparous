#!/usr/bin/perl
# oma-group-to-long.pl turns an oma file into a long format suitable for use in R
# and is released under the terms of the GNU GPL version 3, or any
# later version, at your option. See the file README and COPYING for
# more information.
# Copyright 2014 by Don Armstrong <don@donarmstrong.com>.


use warnings;
use strict;

use Getopt::Long;
use Pod::Usage;

=head1 NAME

oma-group-to-long.pl - turns an oma file into a long format suitable for use in R

=head1 SYNOPSIS

oma-group-to-long.pl [options]

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

oma-group-to-long.pl

=cut


use vars qw($DEBUG);

use IO::Uncompress::Gunzip;
use IO::File;
use Term::ProgressBar::IO;

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
if (not @ARGV==2) {
    push @USAGE_ERRORS,"You must provide exactly one oma group file and one oma ensemble mapping file";
}

pod2usage(join("\n",@USAGE_ERRORS)) if @USAGE_ERRORS;

my $mapping = IO::File->new($ARGV[1]) or
    die "Unable to open mapping file for reading: $!";

my %mapping;
while (<$mapping>) {
    next if /^#/;
    chomp;
    my ($group,$ensembl) = split /\t/;
    next unless defined $group;
    push @{$mapping{$group}},$ensembl;
}
close($mapping);

my $fh = IO::File->new($ARGV[0]) or
    die "Unable to open file for reading: $!";

my $p = Term::ProgressBar::IO->new($fh);
print "group_num\tfingerprint\toma_entry\tensembl\n";
while (<$fh>) {
    next if /^#/;
    chomp;
    my ($group,$fingerprint,@oma_entries) = split /\t/;
    if ($fingerprint eq "n/a") {
        $fingerprint = "NA";
    }
    next unless @oma_entries;
    for my $entry (@oma_entries) {
        next unless defined $mapping{$entry};
        for my $ensembl (@{$mapping{$entry}}) {
            print join("\t",$group,$fingerprint,$entry,$ensembl),"\n";
        }
    }
    $p->update();
}



__END__
