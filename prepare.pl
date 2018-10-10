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
  -g|--genes        [FILE]   : MCScanX genes annotation file
  -t|--te1          [FILE]   : TE annotation file (species 1) (GFF3)
  -T|--te2          [FILE]   : TE annotation file (species 1) (GFF3)
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
  $collinearity_infile,
  $genes_infile,
  $score_infile,
  $nnns_infile,
  $coverage_infile,
  $help, $debug
);
## defaults
my $ks_threshold = 0.2;
my $outprefix = "prepare";

GetOptions (
  't|tes=s' => \$tes_infile,
  'c|collinearity=s' => \$collinearity_infile,
  'g|genes=s' => \$genes_infile,
  's|score:s' => \$score_infile,
  'n|nnns:s' => \$nnns_infile,
  'v|coverage:s' => \$coverage_infile,
  'k|ks:f' => \$ks_threshold,
  'o|outprefix:s' => \$outprefix,
  'h|help' => \$help,
  'd|debug' => \$debug
);

die $usage if $help;
die $usage unless ($tes_infile && $collinearity_infile && $genes_infile);

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
      $scores_hash{$F[0]} = { 'chrom1' => $F[1], 'chrom2' => $F[2], 'orientation' => $F[6], 'ks_avg' => $F[11] }; ##this will include LCBs within species...
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
print STDERR "[INFO] MCScanX collinearity file: $collinearity_infile\n";
open (my $COLL, $collinearity_infile) or die $!;
while (<$COLL>) {
  if ($_ =~ m/^#/) {
    next;
  } else {
    my @F = split (m/\s+/, $_);
    if ($scores_hash{$F[0]}) { ##analyse blocks with Ks <= $ks_threshold
      $genes_hash{$F[2]}++; ##note genes that form LCB
      $genes_hash{$F[3]}++;
      push ( @{ $collinearity_hash{$F[0]}{'genes1'} }, $F[2]); ##anon array of member genes from species1
      push ( @{ $collinearity_hash{$F[0]}{'genes2'} }, $F[3]);
    } else {
      next;
    }
  }
}
close $COLL;
print STDERR "[INFO] Number of genes in LCBs: ".scalar(keys %genes_hash)."\n";

## open ideogram and cytoband files
open (my $IDEOGRAM, ">".$outprefix."_ideogram.txt") or die $!;
print $IDEOGRAM join ("\t", "chr", "start", "end", "\n");
open (my $CYTOBANDS, ">".$outprefix."_cytobands.txt") or die $!;
print $CYTOBANDS join ("\t", "chr", "start", "end", "name", "gieStain", "\n");

## process $collinearity_hash
foreach my $block (sort {$a<=>$b} keys %collinearity_hash) {
  print STDERR "[INFO] Block number: $block\n";
  print STDERR "[INFO] LCB orientation: $scores_hash{$block}{'orientation'}\n";

  ## CODE BLOCK
  ## print consecutively each chr in LCB
  ## first do for 'genes1'
  ## =====================
  ## get start and end genes in LCB array
  my $start1 = ${ $collinearity_hash{$block}{'genes1'} }[0];
  my $end1 = ${ $collinearity_hash{$block}{'genes1'} }[-1];
  print STDERR "[INFO] Start1: $start1\n[INFO] End1: $end1\n";
  ## slice from MCScanX genes file to get coordinates
  if ($scores_hash{$block}{'orientation'} eq "plus") {
    `perl -e 'while (<>) {print if (/\Q$start1\E/../\Q$end1\E/)}' $genes_infile > tmp1`;
  } else {
    ## switch orientation for minus strand LCBs
    `perl -e 'while (<>) {print if (/\Q$end1\E/../\Q$start1\E/)}' $genes_infile > tmp1`;
  }
  ## parse tmp file to get LCB coords as ideogram and gene coords as cytobands
  my %ideogram;
  open (my $TMP1, "tmp1") or die $!;
  while (<$TMP1>) {
    my @F = split (m/\s+/, $_);
    print $CYTOBANDS join ("\t", join("_","LCB$block",$F[0]), -1e+9, 1e+9, "background", "gneg", "\n") if $. == 1; ##print arbitrarily large blank cytoband for each block for visualisation
    ## print genes to cytobands file
    if ($genes_hash{$F[1]}) { ##gene is part of LCB
      print $CYTOBANDS join ("\t", join("_","LCB$block",$F[0]), $F[2], $F[3], $F[1], "stalk", "\n");
    } else {
      print $CYTOBANDS join ("\t", join("_","LCB$block",$F[0]), $F[2], $F[3], $F[1], "gpos25", "\n");
    }

    ## get all coords of all genes in region
    push (@{$ideogram{$F[0]}}, $F[2]);
    push (@{$ideogram{$F[0]}}, $F[3]);
  }
  close $TMP1;

  ## then do for 'genes2'
  ## ====================
  my $start2 = ${ $collinearity_hash{$block}{'genes2'} }[0];
  my $end2 = ${ $collinearity_hash{$block}{'genes2'} }[-1];
  print STDERR "[INFO] Start2: $start2\n[INFO] End2: $end2\n";
  if ($scores_hash{$block}{'orientation'} eq "plus") {
    `perl -e 'while (<>) {print if (/\Q$start2\E/../\Q$end2\E/)}' $genes_infile > tmp2`;
  } else {
    `perl -e 'while (<>) {print if (/\Q$end2\E/../\Q$start2\E/)}' $genes_infile > tmp2`;
  }
  open (my $TMP2, "tmp2") or die $!;
  while (<$TMP2>) {
    my @F = split (m/\s+/, $_);
    print $CYTOBANDS join ("\t", join("_","LCB$block",$F[0]), -1e+9, 1e+9, "background", "gneg", "\n") if $. == 1; ##print arbitrarily large blank cytoband for each block for visualisation
    if ($genes_hash{$F[1]}) {
      print $CYTOBANDS join ("\t", join("_","LCB$block",$F[0]), $F[2], $F[3], $F[1], "stalk", "\n");
    } else {
      print $CYTOBANDS join ("\t", join("_","LCB$block",$F[0]), $F[2], $F[3], $F[1], "gpos25", "\n");
    }
    push (@{$ideogram{$F[0]}}, $F[2]);
    push (@{$ideogram{$F[0]}}, $F[3]);
  }
  close $TMP2;
  foreach (nsort keys %ideogram) {
    print $IDEOGRAM join ("\t", join("_","LCB$block",$_), ${$ideogram{$_}}[0], ${$ideogram{$_}}[-1], "\n");
  }

}
close $IDEOGRAM;
close $CYTOBANDS;


__END__
