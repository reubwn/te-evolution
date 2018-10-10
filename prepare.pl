#!/usr/bin/env perl

##
## prepare.pl
## reubwn 2018
##

use strict;
use warnings;
use Getopt::Long;
use Sort::Naturally;
use Data::Dumper;

my $usage = "
SYNOPSIS

OPTIONS
  -c|--collinearity [FILE]   : MCScanX collinearity file (reformatted!)
  -s|--score        [FILE]   : MCScanX score file
  -t|--te1          [FILE]   : TE annotation file (species 1) (GFF3)
  -T|--te2          [FILE]   : TE annotation file (species 1) (GFF3)
  -g|--gff1         [FILE]   : Genes annotation file (species 1) (GFF3)
  -G|--gff2         [FILE]   : Genes annotation file (species 2) (GFF3)
  -n|--nnns         [FILE]   : NNNs annotation file (GFF3) [!!TODO]
  -v|--coverage     [FILE]   : coverage annotation file [!!TODO]
  -k|--ks           [FLOAT]  : Ks threshold to filter LCBs
  -o|--out          [STRING] : Output filename
  -h|--help                  : this message

OUTPUTS
\n";

## input
my (
  $tes_infile,
  $mcscanx_infile,
  $gff_infile,
  $score_infile,
  $nnns_infile,
  $coverage_infile,
  $help, $debug
);
## defaults
my $ks_threshold = "0.2";
my $outprefix = "prepare";

GetOptions (
  't|tes=s' => \$tes_infile,
  'c|collinearity=s' => \$mcscanx_infile,
  'g|gff=s' => \$gff_infile,
  's|score:s' => \$score_infile,
  'n|nnns:s' => $nnns_infile,
  'c|coverage:s' => $coverage_infile,
  'k|ks:s' => $ks_threshold,
  'o|outprefix:s' => \$outprefix,
  'h|help' => \$help,
  'd|debug' => \$debug
);

die $usage if $help;
die $usage unless ($tes_infile && $mcscanx_infile && $gff_infile);

## stuff
my (
  %scores_hash,
  %collinearity_hash,
  %genes_hash,
  %tes_hash,
  %nnns_hash,
  %coverage_hash
);

## parse scores file
my ($total,$skip) = (0,0);
## open file
open (my $SCORES, $score_infile) or die $!;
while (<$SCORES>) {
  if ($. == 1) {
    next;
  } else {
    my @F = split (m/\s+/, $_);
    if ($F[11] <= $ks_threshold) {
      $scores_hash{$F[0]} = { 'chrom1' => $F[1], 'chrom2' => $F[2], 'ks_avg' => $F[11] }; ##this will include LCBs within species...
    } else {
      $skip++;
    }
  }
  $total++;
}
close $SCORES;
print STDERR "[INFO] Scores file: $score_infile\n";
print STDERR "[INFO] Number of LCBs passing Ks threshold (<= $ks_threshold): ".(scalar(keys %scores_hash))."\n";

## parse collinearity file
open (my $COLL, $mcscanx_infile) or die $!;
while (<$COLL>) {
  if ($_ =~ m/^#/) {
    next;
  } else {
    my @F = split (m/\s+/, $_);
    if ($scores_hash{$F[0]}) { ##analyse blocks with Ks <= $ks_threshold
      push ( @{ $collinearity_hash{$F[0]}{'genes1'} }, $F[2]); ##anon array of member genes from species1
      push ( @{ $collinearity_hash{$F[0]}{'genes2'} }, $F[3]);
    } else {
      next;
    }
  }
}
close $COLL;
print STDERR "[INFO] MCScanX collinearity file: $mcscanx_infile\n";

## open ideogram and cytoband files
open (my $IDEOGRAM, ">".$outprefix."_ideogram.txt") or die $!;
print $IDEOGRAM join ("\t", "chr", "start", "end", "\n");
open (my $CYTOBANDS, ">".$outprefix."_cytobands.txt") or die $!;
print $CYTOBANDS join ("\t", "chr", "start", "end", "name", "gieStain", "\n");

## process $collinearity_hash
foreach (sort {$a<=>$b} keys %collinearity_hash) {
  ## get start and end genes in LCB array
  my $start = ${ $collinearity_hash{$_}{'genes1'} }[0];
  my $end = ${ $collinearity_hash{$_}{'genes1'} }[-1];
  ## slice from GFF to get coordinates
  `perl -e 'while (<>) {print if (/\Q$start\E/../\Q$end\E/)}' $gff_infile | grep mRNA > tmp1`;
  my %ideogram;
  open (my $TMP, "tmp1") or die $!;
  while (<$TMP>) {
    my @F = split (m/\s+/, $_);
    ## regex to get the 'ID' field from GFF as name, assumes '|' delim with Augustus style gene names eg: 'ABC|g1.t1' etc
    if ($F[8] =~ m/ID\=(\w+\|g\d+\.t\d+);.+/) {
      print $CYTOBANDS join ("\t", $F[0], $F[3], $F[4], $1, "\n");
    } else {
      ## otherwise set name to 9th field in GFF
      print $CYTOBANDS join ("\t", $F[0], $F[3], $F[4], $F[8], "\n");
    }
    ## get all coords of all genes in region
    push (@{$ideogram{$F[0]}}, $F[3]);
    push (@{$ideogram{$F[0]}}, $F[4]);
  }
  ## print regions to ideogram file, setting start and end coords as 1st bp and last bp of genes in that region...
  foreach (nsort keys %ideogram) {
    print $IDEOGRAM join ("\t", $_, ${$ideogram{$_}}[0], ${$ideogram{$_}}[-1], "\n");
  }
}
close $IDEOGRAM;
close $CYTOBANDS;


__END__
