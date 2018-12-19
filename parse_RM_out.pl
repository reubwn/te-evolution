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
  -i|--infile [FILE] : RepeatMasker out file (*.out)
  -f|--fasta  [FILE] : Genome fasta file
  -g|--glob          : Will glob all *.out files and attempt to process them
  -n|--nocondense    : Classify based on subfamilies [family]
  -h|--help          : this message

OUTPUTS
\n";

## input
my (
  $RM_infile,
  $fasta_infile,
  $glob,
  $no_condense,
  $help, $debug
);
## defaults
my $outprefix = "plot_RM";

GetOptions (
  'i|infile=s' => \$RM_infile,
  'f|fasta:s' => \$fasta_infile,
  'g|glob' => \$glob,
  'h|help' => \$help,
  'debug' => \$debug
);
## help and usage
die $usage if $help;
print STDERR "[####] TE-EVOLUTION parse_RM_out.pl\n";
print STDERR "[####] " . `date`;

###############
## single files
###############
if ($RM_infile) {

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

######################
## glob multiple files
######################
} elsif ($glob) {

  ## STUFF
  my %files_hash;
  my %repeats_hash;
  my %genome_lengths_hash;

  my @files = glob "*.out";
  print STDERR "[####] Number of RM files: ".@files."\n";

  foreach my $file (@files) {
    print STDERR "[####] Working on file: $file\n";
    ## STUFF
    my %print_hash;
    (my $fasta_infile = $file) =~ s/\.out//;
    my $genome_length;

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
    $genome_lengths_hash{$file} = $genome_length;

    open (my $IN, $file) or die $!;
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
        $repeats_hash{$rep}++;
      }
    }
    close $IN;
    $files_hash{$file} = \%print_hash;
  }

  print STDOUT join ("\t", "\t", nsort keys %files_hash);
  foreach my $repeat (nsort keys %repeats_hash) {
    print STDOUT "$repeat\t";
    foreach my $file (nsort keys %files_hash) {
      my %print_hash = %{ $files_hash{$file} };
      if ($print_hash{$repeat} ) {
        print STDOUT (($print_hash{$repeat}/$genome_lengths_hash{$file})*100) . "\t";
      } else {
        print STDOUT "0\t";
      }
    }
    print STDOUT "\n";
  }
}



############ SUBS

sub trim {
  my $s = shift;
  $s =~ s/^\s+|\s+$//g;
  return $s
}

__END__
