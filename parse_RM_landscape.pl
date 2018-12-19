#!/usr/bin/env perl

##
## parse_RM_landscape.pl
## reubwn 2018
##

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);;
use Sort::Naturally;
use List::Util qw( min max );
use Data::Dumper;

my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -i|--infile     [FILE]   : RepeatMasker landscape file (HTML) [*]
  -o|--out        [STRING] : Prefix for outfiles [landscape]
  -n|--nocondense          : Print full results [no]
  -h|--help                : this message

OUTPUTS
\n";

## input
my (
  $html_infile,
  $no_condense,
  $help, $debug
);
## defaults
my $outprefix = "landscape";

GetOptions (
  'i|infile=s' => \$html_infile,
  'n|nocondense' => \$no_condense,
  'o|outprefix:s' => \$outprefix,
  'h|help' => \$help,
  'debug' => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ($html_infile);
print STDERR "[####] TE-EVOLUTION parse_RM_landscape.pl\n";
print STDERR "[####] " . `date`;

## STUFF
my (%html_hash, %te_full_hash, %te_class_hash, %print_hash);
my ($i,$j) = (0,0);

## parse HTML file
open (my $HTML, $html_infile) or die $!;
while (<$HTML>) {
  chomp;
  my $class;
  if ($_ =~ m/data.addColumn\('number',\s'(\w+\/*.*)'\);/) { ## get TEs present
    # print "$1\n";
    $te_full_hash{$i} = $1;
    ($class = $1) =~ s/\/.*//; ## trim family
    $te_class_hash{$i} = $class; ## put into hash
    $i++;
  } elsif ($_ =~ m/\['\d{1,2}',\s(.+)\],/) { ## get coverage of each TE
    my @values = split (/,\s/, $1);
    for my $k (0..$#values) {
      $html_hash{$j}{$te_class_hash{$k}} += $values[$k];
    }
    $j++;
  }
}
close $HTML;
print STDERR "[####] Number of TE families: ".scalar(keys %te_full_hash)."\n";

## Dumper
print Dumper(\%html_hash) if $debug;

## open outfile
open (my $OUT, ">$html_infile.class_landscape") or die $!;

## print condensed data
print $OUT join ("\t", "Divergence", "Unknown", "Other", "DNA", "RC", "LTR", "LINE", "SINE" ) . "\n";
foreach my $col (sort {$a<=>$b} keys %html_hash) {
  print $OUT "$col\t";
  print $OUT join ("\t", ($html_hash{$col}{Unknown} ? $html_hash{$col}{Unknown} : 0),
                          ($html_hash{$col}{Other} ? $html_hash{$col}{Other} : 0),
                          ($html_hash{$col}{DNA} ? $html_hash{$col}{DNA} : 0),
                          ($html_hash{$col}{RC} ? $html_hash{$col}{RC} : 0),
                          ($html_hash{$col}{LTR} ? $html_hash{$col}{LTR} : 0),
                          ($html_hash{$col}{LINE} ? $html_hash{$col}{LINE} : 0),
                          ($html_hash{$col}{SINE} ? $html_hash{$col}{SINE} : 0)
                  ) . "\n";
}
close $OUT;

__END__
