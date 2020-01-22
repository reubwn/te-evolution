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
  -i|--infile  [FILE]* : RepeatMasker out file (*.out)
  -f|--fasta   [FILE]* : Genome fasta file
  -m|--mapping [FILE]  : Mapping file with genome sizes for each file
  -g|--glob            : Will glob all *.out files and attempt to process them
  -n|--nocondense      : Classify based on subfamilies [family]
  -h|--help            : this message

OUTPUTS
\n";

## input
my (
  $RM_infile,
  $fasta_infile,
  $mapping_infile,
  $glob,
  $no_condense,
  $help, $debug
);
## defaults
my $outprefix = "repeats_tab";

GetOptions (
  'i|infile=s' => \$RM_infile,
  'f|fasta:s' => \$fasta_infile,
  'm|mapping:s' => \$mapping_infile,
  'g|glob' => \$glob,
  'h|help' => \$help,
  'debug' => \$debug
);
## help and usage
die $usage if $help;
print STDERR "[####] TE-EVOLUTION parse_RM_out.pl\n";
print STDERR "[####] " . `date`;

## custom order to print results
my @sort_order = qw (
 Unknown Other
 DNA DNA/Academ DNA/CMC DNA/Crypton DNA/Dada DNA/Ginger DNA/IS3EU DNA/IS5 DNA/Kolobok DNA/MULE DNA/Maverick DNA/Merlin DNA/P DNA/PIF DNA/PiggyBac DNA/Sola DNA/TcMar DNA/Zator DNA/hAT
 RC RC/Helitron
 LTR LTR/Copia LTR/DIRS LTR/Gypsy LTR/Juno LTR/Mag LTR/Ngaro LTR/Pao LTR/TelKA
 LINE/Penelope
 LINE LINE/CR1 LINE/CRE LINE/Dong LINE/I LINE/L1 LINE/L2 LINE/Proto2 LINE/R1 LINE/R2 LINE/RTE LINE/Rex LINE/Soliton
 SINE SINE/5S SINE/I SINE/R2 SINE/tRNA
 tRNA rRNA Satellite Simple_repeat Low_complexity ARTEFACT RNA Retroposon repeat
);
my $count = 0;
my %position_of;
$position_of{$_} = $count++ for @sort_order;

## open outfile
open (my $OUTFILE, ">".$outprefix.".txt") or die $!;

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
    close $FASTA;
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
    print $OUTFILE join ("\t", $_, $print_hash{$_});
    if ($fasta_infile) {
      print $OUTFILE "\t" . (($print_hash{$_}/$genome_length)*100) . "\n";
    } else {
      print $OUTFILE "\n";
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

  if ( $mapping_infile ) {
    open (my $MAP, $mapping_infile) or die $!;
    while (<$MAP>) {
      chomp;
      my @F = split (m/\s+/, $_);
      ## filename (*.out) in col1; genome span in col2
      $genome_lengths_hash{$F[0]} = $F[1];
    }
    close $MAP;
  }

  foreach my $file (@files) {
    print STDERR "[####] Working on file: $file\n";
    ## STUFF
    my %print_hash;
    (my $fasta_infile = $file) =~ s/\.out//;
    my $genome_length;

    ## open fasta only if mapping file not provided
    unless ( $mapping_infile ) {
      open (my $FASTA, $fasta_infile) or die $!;
      while (<$FASTA>) {
        if ($_ =~ m/^\>/) {
          next;
        } else {
          chomp;
          $genome_length += length ($_);
        }
      }
      close $FASTA;

      print STDERR "[####] Genome length: $genome_length bp\n";
      $genome_lengths_hash{$file} = $genome_length;
  } else {
    print STDERR "[####] Genome length: $genome_lengths_hash{$file} bp\n";
  }

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

  print $OUTFILE join ("\t", "repeat", nsort keys %files_hash) . "\n";
  foreach my $repeat (nsort keys %repeats_hash) {
    print $OUTFILE "$repeat\t";
    foreach my $file (nsort keys %files_hash) {
      my %print_hash = %{ $files_hash{$file} };
      if ($print_hash{$repeat} ) {
        print $OUTFILE (($print_hash{$repeat}/$genome_lengths_hash{$file})*100) . "\t";
      } else {
        print $OUTFILE "0\t";
      }
    }
    print $OUTFILE "\n";
  }
}

close $OUTFILE;

############ SUBS

sub trim {
  my $s = shift;
  $s =~ s/^\s+|\s+$//g;
  return $s
}

sub by_order {
#  my @a = split m/\//, $a;
#  my @b = split m/\//, $b;
  $position_of{$a} <=> $position_of{$b};
}


__END__
