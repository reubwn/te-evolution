#!/usr/bin/env perl

##
## parse_RM_landscape.pl
## reubwn 2018
##

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -i|--infile [FILE] : RepeatMasker out file (*.out) [*]
  -f|--fasta  [FILE] : Genome fasta file
  -n|--nocondense    : Classify based on subfamilies [family]
  -h|--help          : this message

OUTPUTS
\n";

## input
my (
  $RM_infile,
  $fasta_infile,
  $no_condense,
  $help, $debug
);
## defaults
my $outprefix = "plot_RM";

GetOptions (
  'i|infile=s' => \$RM_infile,
  'f|fasta:s' => \$fasta_infile,
  'h|help' => \$help,
  'debug' => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ($RM_infile);
print STDERR "[####] TE-EVOLUTION parse_RM_out.pl\n";
print STDERR "[####] " . `date`;

## STUFF
my %print_hash;
my $genome_length;

if ($fasta_infile) {
  open (my $FASTA, $fasta_infile) or die $!;
  while (<$FASTA>) {
    if ($_ =~ m/^\>/) {
      next;
    } else {
      chomp;
      $genome_length += length ($_);
    }
  }
  print STDERR "[####] Genome length: $genome_length bp\n";
}

open (my $IN, $RM_infile) or die $!;
while (<$IN>) {
  chomp;
  my $line = trim ($_); ## remove leading/trailing whitespace
  if ($line !~ m/^\d+/) {
    next;
  } else {
    my @F = split (/\s+/, $line);
    my $rep = $F[10];
    unless ($no_condense) {
      $rep =~ s/-.*//;
      $rep =~ s/\?//g;
    }
    $print_hash{$rep} += ( $F[6] - $F[5] + 1 );
  }
}
close $IN;

foreach (sort { "\L$a" cmp "\L$b" || $print_hash{$b} <=> $print_hash{$a} } keys %print_hash) {
  print STDOUT join ("\t", $_, $print_hash{$_});
  if ($fasta_infile) {
    print "\t" . (($print_hash{$_}/$genome_length)*100) . "\n";
  } else {
    print "\n";
  }
}

############ SUBS

sub trim {
  my $s = shift;
  $s =~ s/^\s+|\s+$//g;
  return $s
}

__END__
