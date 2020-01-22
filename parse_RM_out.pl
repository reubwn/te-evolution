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
  Globs *.out (or *.out.gz) RM files and summarizes counts across all input files.
  Requires a genome span mapping file in 'filename [TAB] genome_size' format.

OPTIONS [*] = required
  -m|--mapping [FILE]  : mapping file with genome sizes for each file
  -z|--gz              : gunzip any *.out.gz files prior to glob
  -n|--nocondense      : classify based on subfamilies [family]
  -h|--help            : this message
\n";

## input
my (
  $mapping_infile,
  $gzip,
  $no_condense,
  $help, $debug
);
## defaults
my $outprefix = "repeats_tab";

GetOptions (
  'm|mapping=s' => \$mapping_infile,
  'z|gz'        => \$gzip,
  'h|help'      => \$help,
  'debug'       => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ( $mapping_infile );

print STDERR "[####] TE-EVOLUTION parse_RM_out.pl\n";
print STDERR "[####] " . `date`;

## open outfile
open (my $OUTFILE, ">".$outprefix.".txt") or die $!;

## gunzip the files
if ( $gzip ) {
  my @gzipped_files = glob "*.out.gz";
  print STDERR "[INFO] Gunzip on ".scalar(@gzipped_files)." files\n";
  if (system ("parallel --help &>/dev/null") != 0){
    foreach my $file (@gzipped_files) {
      `gunzip $file`;
    }
  } else {
    `ls *.out.gz | parallel gunzip {}`;
  }
}

## stuff
my %files_hash;
my %repeats_hash;
my %genome_lengths_hash;

my @files = glob "*.out";
print STDERR "[INFO] Number of RM files: ".@files."\n";

## open mapping file to get genome sizes
open (my $MAP, $mapping_infile) or die $!;
while (<$MAP>) {
  chomp;
  my @F = split (m/\s+/, $_);
  ## filename (*.out) in col1; genome span in col2
  $genome_lengths_hash{$F[0]} = $F[1];
}
close $MAP;

foreach my $file (@files) {
  print STDERR "[INFO] Working on file: $file\n";
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

    print STDERR "[INFO] Genome length: $genome_length bp\n";
    $genome_lengths_hash{$file} = $genome_length;
} else {
  print STDERR "[INFO] Genome length: $genome_lengths_hash{$file} bp\n";
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

close $OUTFILE;

## gunzip the files
if ( $gzip ) {
  print STDERR "[INFO] Gzip on ".scalar(@files)." files\n";
  if (system ("parallel --help &>/dev/null") != 0){
    foreach my $file (@files) {
      `gzip $file`;
    }
  } else {
    `ls *.out | parallel gzip {}`;
  }
}

print STDERR "[INFO] Finished! " . `date`;

############ SUBS

sub trim {
  my $s = shift;
  $s =~ s/^\s+|\s+$//g;
  return $s
}

__END__
