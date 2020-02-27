#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use File::Basename;
use Sort::Naturally;
use Getopt::Long qw(:config no_ignore_case);

## TODO

my $usage = "
SYNOPSIS
  Extract sequences for full-length LTRs

OPTIONS [*] = required
  -g|--gff     [FILE]   : GFF file of full-length LTRs (GFF) [*]
  -f|--fasta   [FILE]   : fasta file of genome(s) (fasta) [*]
  -e|--feature [STRING] : select GFF features from 3rd field
  -o|--out     [STRING] : output prefix for outfiles ['extraction']
  -h|--help             : this message
\n";

## input
my ( $gff_file, $fasta_file, $gff_feature, $help, $debug );
## defaults
my $outprefix = "extraction";

GetOptions (
  'i|input=s' => \$gff_file,
  'f|fasta=s' => \$fasta_file,
  'o|out:s' => \$outprefix,
  'e|feature:s' => \$gff_feature,
  'h|help' => \$help,
  'debug' => \$debug
  );
## help and usage
die $usage if $help;
die $usage unless ( $gff_file && $fasta_file );

print STDERR "[####] TE-EVOLUTION extract_LTRs.pl\n";
print STDERR "[####] " . `date`;

my (%scaffolds_hash);

## parse fasta
my $in = Bio::SeqIO -> new ( -file => $fasta_file, -format => "fasta" );
while (my $seq_obj = $in->next_seq) {
  $scaffolds_hash{$seq_obj->display_id()} = $seq_obj;
}
print STDERR "[INFO] Number of sequences in '$fasta_file': ".commify(scalar(keys %scaffolds_hash))."\n";

my $OUT;
if ( $gff_feature ) {
  print STDERR "[INFO] Extracting regions matching GFF feature '$gff_feature'\n";
  ## open outfile
  open ($OUT, ">$outprefix.$gff_feature.fasta") or die $!;
} else {
  print STDERR "[INFO] Extracting all regions from GFF\n";
  ## open outfile
  open ($OUT, ">$outprefix.all_features.fasta") or die $!;
}

## parse GFF3
open (my $GFF, $gff_file) or die $!;
while (my $line = <$GFF>) {
  chomp $line;
  next if ($line =~ m/^\#/);
  my @F = split (m/\s+/, $line);
  if ( $gff_feature ) {
    if ( $F[2] eq $gff_feature ) {
      print $OUT ">$F[0]:$F[3]..$F[4]:$F[2]\n";
      print $OUT $scaffolds_hash{$F[0]} -> subseq($F[3],$F[4]) . "\n";
    } else {
      next;
    }
  } else {
    print $OUT ">$F[0]:$F[3]..$F[4]:$F[2]\n";
    print $OUT $scaffolds_hash{$F[0]} -> subseq($F[3],$F[4]) . "\n";
  }
}
close $OUT;

############### SUBS

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

__END__
