#!/usr/bin/env perl

##
## parse_RM_out_diversity.pl
## reubwn 2020
##

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS
  Parse RepeatMasker *.out file to report TE diversity at family and superfamily level.

OPTIONS [*] = required
  -i|--infile [FILE]*  : mapping file with genome sizes for each file
  -h|--help            : this message
\n";

## input
my (
  $RM_infile,
  $help, $debug
);
## defaults
my $outsuffix = "diversity";

GetOptions (
  'i|infile=s'  => \$RM_infile,
  'h|help'      => \$help,
  'debug'       => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ( $RM_infile );

print STDERR "[####] TE-EVOLUTION parse_RM_out_diversity.pl\n";
print STDERR "[####] " . `date`;

(my $sample = $RM_infile) =~ s/\.out.*//;
my %results_hash;

open(my $IN, $RM_infile) or die $!;
while (my $line=<$IN>) {
  next if $. < 4; ## skip first 3 lines
  $line =~ s/^\s+|\s+$//g; ## trim that leading whitespace
  my @F = split (m/\s+/, $line);

  ## get class
  my @a = split (m/\//, $F[10]);
  my $class = $a[0];
  ## get superfamily
  ## think this is more consistently provided by the prefix to the left of the '-' in col 10
  my $superfamily = $1 if $F[9] =~ m/(\w+)\-/;
  ## get family
  my $family = $1 if $F[9] =~ m/\-(\w+)\_/;

  if ($class =~ m/DNA|RC|LTR|LINE|SINE/) { ## skip unknown, simple etc.
    $results_hash{$sample}{$class}{superfamily}++ if ( $superfamily );
    $results_hash{$sample}{$class}{family}++ if ( $family );

    # if ( ($superfamily) && ($family) ) {
    #   # print join ("\t", $class,$superfamily,$family) . "\n";
    #   $results_hash{$sample}{$class}{$superfamily}{$family}++
    # }
  }
}
close $IN;

foreach my $sample (nsort keys %results_hash) {
  print "$sample\t";
  foreach my $class (nsort keys %{$results_hash{$sample}}) {
    print $$results_hash{$sample}{$class}{superfamily} . "\n";
  }
}


__END__
