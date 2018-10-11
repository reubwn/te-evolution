#!/usr/bin/env perl

##
## prepare.pl
## reubwn 2018
##

use strict;
use warnings;

use Getopt::Long;
use Sort::Naturally;
use List::Util qw( min max );
use Data::Dumper;

use Bio::Seq;
use Bio::SeqIO;

my $usage = "
SYNOPSIS

OPTIONS
  -c|--collinearity [FILE]   : MCScanX collinearity file (reformatted!) [*]
  -s|--score        [FILE]   : MCScanX score file [*]
  -g|--genes        [FILE]   : MCScanX genes annotation file [*]
  -f|--fasta        [FILE]   : Genome fasta file [*]
  -t|--tes1         [FILE]   : TE annotation file (species 1) (GFF3) [*]
  -T|--tes2         [FILE]   : TE annotation file (species 1) (GFF3)
  -n|--nnns         [FILE]   : NNNs annotation file (GFF3) [!!TODO]
  -v|--coverage     [FILE]   : coverage annotation file [!!TODO]
  -k|--ks           [FLOAT]  : Ks threshold to filter LCBs [0.2]
  -d|--find         [STRING] : Search string to grep TE id from GFF file [Family]
  -o|--out          [STRING] : Output filename
  -h|--help                  : this message
[*] Required input

OUTPUTS
\n";

## input
my (
  $collinearity_infile,
  $score_infile,
  $genes_infile,
  $genome_infile,
  $tes_infile1,
  $tes_infile2,
  $nnns_infile,
  $coverage_infile,
  $help, $debug
);
## defaults
my $ks_threshold = 0.2;
my $outprefix = "prepare";
my $find = "Family";

GetOptions (
  'c|collinearity=s' => \$collinearity_infile,
  's|score=s' => \$score_infile,
  'g|genes=s' => \$genes_infile,
  'f|fasta=s' => \$genome_infile,
  't|tes1=s' => \$tes_infile1,
  'T|tes2:s' => \$tes_infile2,
  'n|nnns:s' => \$nnns_infile,
  'v|coverage:s' => \$coverage_infile,
  'k|ks:f' => \$ks_threshold,
  'd|find:s' => \$find,
  'o|outprefix:s' => \$outprefix,
  'h|help' => \$help,
  'd|debug' => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ($collinearity_infile && $score_infile && $genes_infile && $genome_infile && $tes_infile1);

## stuff
my (
  %scores_hash,
  %collinearity_hash,
  %genes_hash,
  %genome_hash,
  %tes_hash,
  %nnns_hash,
  %coverage_hash
);

## parse scores file
my ($total,$skip) = (0,0);
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

## parse genome fasta file
my $in = Bio::SeqIO -> new( -file => $genome_infile, -format => "fasta" );
while ( my $seqobj = $in->next_seq() ) {
  $genome_hash{ $seqobj->display_name() } = $seqobj->length();
}

## make ideogram genome file
open (my $IDEOGRAM_GENOME, ">".$outprefix."_ideogram.genome.txt") or die $!;
open (my $CYTOBANDS_GENOME, ">".$outprefix."_cytobands.genome.txt") or die $!;
print $CYTOBANDS_GENOME join ("\t", "chr", "start", "end", "name", "gieStain") . "\n";
print $IDEOGRAM_GENOME join ("\t", "chr", "start", "end") . "\n";
foreach (nsort keys %genome_hash) {
  print $IDEOGRAM_GENOME join ("\t", $_, "1", $genome_hash{$_}) . "\n";
}
close $IDEOGRAM_GENOME;

## open cytobands and ideogram LCBs file
open (my $IDEOGRAM_LCB, ">".$outprefix."_ideogram.LCBs.txt") or die $!;
open (my $CYTOBANDS_LCB, ">".$outprefix."_cytobands.LCBs.txt") or die $!;
open (my $REPEATS_LCB, ">".$outprefix."_repeats.LCBs.txt") or die $!;
print $IDEOGRAM_LCB join ("\t", "chr", "start", "end", "name") . "\n";
print $CYTOBANDS_LCB join ("\t", "chr", "start", "end", "name", "gieStain") . "\n";
print $REPEATS_LCB join ("\t", "chr", "start", "end", "strand", "name") . "\n";


## process $collinearity_hash
foreach my $block (sort {$a<=>$b} keys %collinearity_hash) {
  print STDERR "[INFO] Block number: $block\n";
  print STDERR "[INFO] LCB orientation: $scores_hash{$block}{'orientation'}\n";

  ## MAIN BLOCK
  ## print consecutively each chr in LCB
  ## first do for 'genes1'
  ## =====================
  ## get start and end genes in LCB array
  my $start1 = ${ $collinearity_hash{$block}{'genes1'} }[0];
  my $end1 = ${ $collinearity_hash{$block}{'genes1'} }[-1];
  print STDERR "[INFO] Start1: $start1\n[INFO] End1: $end1\n";
  ## slice from MCScanX genes file to get coordinates
  `perl -e 'while (<>) {print if (/\Q$start1\E/../\Q$end1\E/)}' $genes_infile > tmp1`;
  ## parse tmp file to get LCB coords as ideogram and gene coords as cytobands
  my %ideogram;
  open (my $TMP1, "tmp1") or die $!;
  while (<$TMP1>) {
    my @F = split (m/\s+/, $_);
    if ($. == 1) { ## this prints arbitrarily large blank cytoband for each block for prettier visualisation
      print $CYTOBANDS_LCB join ("\t", "LCB#$block:1", -1e+5, 1e+5, "background", "gneg") . "\n";
    }
    ## print genes to cytobands file
    if ($genes_hash{$F[1]}) { ##gene is part of LCB
      print $CYTOBANDS_LCB join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "stalk") . "\n";
    } else {
      print $CYTOBANDS_LCB join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "gpos25") . "\n";
    }
    ## get all coords of all genes in region
    $ideogram{"LCB#$block:1"}{chrom} = $F[0]; ##key=LCB##:1; val=chrom
    push ( @{ $ideogram{"LCB#$block:1"}{coords} }, $F[2] ); ##key=LCB##:1; val=@{array of start-end coordinates}
    push ( @{ $ideogram{"LCB#$block:1"}{coords} }, $F[3] );
  }
  close $TMP1;
  ## get TEs that intersect with LCB region using bedtools
  my ($min, $max) = ( (min @{ $ideogram{"LCB#$block:1"}{coords}), (max @{ $ideogram{"LCB#$block:1"}{coords}) );
  my $text = `printf "$ideogram{LCB#$block:1}{chrom}\t(min @{ $ideogram{LCB#$block:1}{coords})\t(max @{ $ideogram{LCB#$block:1}{coords})" | bedtools intersect -a $tes_infile1 -b stdin -wa | perl -lane 'if(m/$find\=([\w]+)(\;.+)*/){print join("\t","LCB#$block:1",$F[3],$F[4],$F[2],$1)}'`;
  print $REPEATS_LCB "$text";

  ## then do for 'genes2'
  ## ====================
  my $start2 = ${ $collinearity_hash{$block}{'genes2'} }[0];
  my $end2 = ${ $collinearity_hash{$block}{'genes2'} }[-1];
  ## check orientation of LCB on second strand! (genes1 is always +)
  if ($scores_hash{$block}{'orientation'} eq "plus") {
    print STDERR "[INFO] Start2: $start2\n[INFO] End2: $end2\n";
    `perl -e 'while (<>) {print if (/\Q$start2\E/../\Q$end2\E/)}' $genes_infile > tmp2`;
  } else {
    print STDERR "[INFO] Start2: $end2\n[INFO] End2: $start2\n"; ##switcheroo
    `perl -e 'while (<>) {print if (/\Q$end2\E/../\Q$start2\E/)}' $genes_infile > tmp2`;
  }
  open (my $TMP2, "tmp2") or die $!;
  while (<$TMP2>) {
    my @F = split (m/\s+/, $_);
    if ($. == 1) {
      print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", -1e+5, 1e+5, "background", "gneg") . "\n";
    }
    if ($genes_hash{$F[1]}) {
      print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "stalk") . "\n";
    } else {
      print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "gpos25") . "\n";
    }
    $ideogram{"LCB#$block:2"}{chrom} = $F[0]; ##
    push ( @{ $ideogram{"LCB#$block:2"}{coords} }, $F[2] );
    push ( @{ $ideogram{"LCB#$block:2"}{coords} }, $F[3] );
  }
  close $TMP2;
  ## check there are 2 chroms in %ideogram
  if (scalar(keys %ideogram) == 2) {
    ## now print the ideogram for each LCB
    foreach (nsort keys %ideogram) {
      print $IDEOGRAM_LCB join ("\t", $_, (min @{$ideogram{$_}{coords}}), (max @{$ideogram{$_}{coords}}), $ideogram{$_}{chrom}) . "\n";
    }
  } else {
    print STDERR "[WARN] LCB#$block does not have two chromosomes\n";
  }
}
## close fhs
close $IDEOGRAM_GENOME;
close $CYTOBANDS_GENOME;
close $IDEOGRAM_LCB;
close $CYTOBANDS_LCB;
close $REPEATS_LCB;


__END__
