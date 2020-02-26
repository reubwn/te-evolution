#!/usr/bin/env perl

##
## parse_RM_landscape.pl
## reubwn 2018
##

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -i|--infile     [FILE]   : landscape file (HTML) [*]
  -o|--out        [STRING] : prefix for outfiles [landscape]
  -c|--condense            : condense results to class level
  -p|--ple                 : do NOT modify 'LINE/Penelope' to 'Penelope'
  -h|--help                : this message

OUTPUTS
\n";

## input
my (
  $html_infile,
  $condense,
  $ple,
  $help, $debug
);
## defaults
my $outprefix = "landscape";

GetOptions (
  'i|infile=s' => \$html_infile,
  'o|outprefix:s' => \$outprefix,
  'c|condense' => \$condense,
  'p|ple' => \$ple,
  'h|help' => \$help,
  'debug' => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ($html_infile);
print STDERR "[####] TE-EVOLUTION parse_RM_landscape.pl\n";
print STDERR "[####] " . `date`;

## STUFF
my (%html_hash, %te_full_hash, %te_hash);
my ($i,$j) = (0,0);

## parse HTML file
open (my $HTML, $html_infile) or die $!;
while (<$HTML>) {
  chomp;
  my $repeat;
  if ($_ =~ m/data.addColumn\('number',\s'(\w+\/*.*)'\);/) { ## get TEs present
    # print "$1\n";
    $te_full_hash{$i} = $1;
    if ($condense) {
      ($repeat = $1) =~ s/\/.*//; ## trim to class
    } else {
      $repeat = $1; ## keep full class/family name
    }
    $te_hash{$i} = $repeat; ## key=index order; val= repeat name
    $i++;
  } elsif ($_ =~ m/\['\d{1,2}',\s(.+)\],/) { ## get coverage of each TE
    my @values = split (/,\s/, $1);
    for my $k (0..$#values) {
      $html_hash{$j}{$te_hash{$k}} += $values[$k]; ## HoH; key= Kimura divergence [0..50]; val= %{key= repeat name at index $k; val= SUM(TE span for each key)}
    }
    $j++;
  }
}
close $HTML;
print STDERR "[####] Number of TE families: ".scalar(keys %te_full_hash)."\n";

## modify LINE/Penelope to Penelope post hoc
unless ( $ple ) {
  foreach (keys %te_hash) {
    if ($te_hash{$_} eq 'LINE/Penelope') {
      $te_hash{$_} = 'Penelope';
    }
  }
  foreach (keys %te_full_hash) {
    if ($te_full_hash{$_} eq 'LINE/Penelope') {
      $te_full_hash{$_} = 'Penelope';
    }
  }
  foreach (keys %html_hash) {
    $html_hash{$_}{'Penelope'} = delete ( $html_hash{$_}{'LINE/Penelope'} );
  }
}

## Dumper
print Dumper(\%te_hash) if $debug;

## open outfile
open (my $OUT, ">$html_infile.landscape.txt") or die $!;

## condensed to Class
if ( $condense ) {
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

} else { ## print full table
  my %nr = map {$_ => 1} values %te_hash;
  print $OUT join ("\t", "\"Divergence\"", (map { qq/"$_"/ } nsort keys %nr) ) . "\n";
  foreach my $col (sort {$a<=>$b} keys %html_hash) {
    print $OUT "$col\t";
    foreach (nsort keys %{$html_hash{$col}}) {
      print $OUT "$html_hash{$col}{$_}\t";
    }
    print $OUT "\n";
  }
}
close $OUT;

__END__
