#!/usr/bin/perl

use warnings;
use strict;
use IPC::System::Simple qw(capture);

while (<>) {
    chomp;
    my $srx = $_;
    my ($srx_short) = $srx =~ m/(SRX\d{3})/;
    my $rsync = capture("rsync",
                        "--list-only",
                        "rsync://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/".
                        "SRX/$srx_short/$srx/",
                       );
    print $_."\n" for $rsync =~ /(SRR\d+)/g;
}
