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
  -1|--tes1         [FILE]   : TE annotation file (species 1) (GFF3) [*]
  -2|--tes2         [FILE]   : TE annotation file (species 2) (GFF3) [*]
  -n|--nnns         [FILE]   : NNNs annotation file (GFF3) [!!TODO]
  -v|--coverage     [FILE]   : coverage annotation file [!!TODO]
  -k|--ks           [FLOAT]  : Ks threshold to filter LCBs [0.2]
  -d|--find         [STRING] : Search string to grep TE id from GFF file [Family]
  -o|--out          [STRING] : Prefix for outfiles
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
  '1|tes1=s' => \$tes_infile1,
  '2|tes2=s' => \$tes_infile2,
  'n|nnns:s' => \$nnns_infile,
  'v|coverage:s' => \$coverage_infile,
  'k|ks:f' => \$ks_threshold,
  'd|find:s' => \$find,
  'o|outprefix:s' => \$outprefix,
  'h|help' => \$help,
  'debug' => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ($collinearity_infile && $score_infile && $genes_infile && $genome_infile && $tes_infile1 && $tes_infile2);

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
print STDERR "[INFO] Number of LCBs passing Ks threshold (<=$ks_threshold): ".commify(scalar(keys %scores_hash))."\n";

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
print STDERR "[INFO] Number of genes in LCBs: ".commify(scalar(keys %genes_hash))."\n";

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
  print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'})"; $|=1;

  ## MAIN BLOCK
  ## print consecutively each chr in LCB
  ## first do for 'genes1'
  ## =====================
  ## get start and end genes in LCB array
  my $start1 = ${ $collinearity_hash{$block}{'genes1'} }[0];
  my $end1 = ${ $collinearity_hash{$block}{'genes1'} }[-1];
  print STDERR "\n[INFO] Start1: $start1\n[INFO] End1: $end1\n" if $debug;
  ## slice from MCScanX genes file to get coordinates
  # `perl -e 'while (<>) {print if (/\Q$start1\E/../\Q$end1\E/)}' $genes_infile > tmp1`;
  ## parse tmp file to get LCB coords as ideogram and gene coords as cytobands
  my %ideogram;
  open (my $TMP1, "perl -e 'while (<>) {print if (/\Q$start1\E/../\Q$end1\E/)}' $genes_infile |") or die $!;
  while (<$TMP1>) {
    my @F = split (m/\s+/, $_);
    if ($. == 1) { ## this prints arbitrarily large blank cytoband for each block for prettier visualisation
      print $CYTOBANDS_LCB join ("\t", "LCB#$block:1", -1e+9, 1e+9, "background", "gneg") . "\n";
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
  my ($chrom1, $min1, $max1) = ($ideogram{"LCB#$block:1"}{chrom}, (min @{ $ideogram{"LCB#$block:1"}{coords} }), (max @{ $ideogram{"LCB#$block:1"}{coords} }) );
  open (my $BED1, "printf '$chrom1\t$min1\t$max1' | bedtools intersect -a $tes_infile1 -b stdin -wa |") or die $!;
  while (<$BED1>) {
    if (m/$find\=([\w]+)(\;.+)*/) {
      my @F = split (/\s+/, $_);
      print $REPEATS_LCB join("\t", "LCB#$block:1", $F[3], $F[4], $F[2], $1) . "\n";
    }
  }
  close $BED1;

  ## then do for 'genes2'
  ## ====================
  my $start2 = ${ $collinearity_hash{$block}{'genes2'} }[0];
  my $end2 = ${ $collinearity_hash{$block}{'genes2'} }[-1];
  ## check orientation of LCB on second strand! (genes1 is always +)
  my $TMP2;
  if ($scores_hash{$block}{'orientation'} eq "plus") {
    print STDERR "[INFO] Start2: $start2\n[INFO] End2: $end2\n" if $debug;
    open ($TMP2, "perl -e 'while (<>) {print if (/\Q$start2\E/../\Q$end2\E/)}' $genes_infile |") or die $!;
  } else {
    print STDERR "[INFO] Start2: $end2\n[INFO] End2: $start2\n" if $debug; ##switcheroo
    open ($TMP2, "perl -e 'while (<>) {print if (/\Q$end2\E/../\Q$start2\E/)}' $genes_infile |") or die $!;
  }
  while (<$TMP2>) {
    my @F = split (m/\s+/, $_);
    if ($. == 1) {
      print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", -1e+9, 1e+9, "background", "gneg") . "\n";
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
  my ($chrom2, $min2, $max2) = ($ideogram{"LCB#$block:2"}{chrom}, (min @{ $ideogram{"LCB#$block:2"}{coords} }), (max @{ $ideogram{"LCB#$block:2"}{coords} }) );
  open (my $BED2, "printf '$chrom2\t$min2\t$max2' | bedtools intersect -a $tes_infile2 -b stdin -wa |") or die $!;
  while (<$BED2>) {
    if (m/$find\=([\w]+)(\;.+)*/) {
      my @F = split (/\s+/, $_);
      print $REPEATS_LCB join("\t", "LCB#$block:2", $F[3], $F[4], $F[2], $1) . "\n";
    }
  }
  close $BED2;

  ## now print LCB ideogram file
  ## check there are 2 chroms in %ideogram
  if (scalar(keys %ideogram) == 2) {
    ## now print the ideogram for each LCB
    foreach (nsort keys %ideogram) {
      print $IDEOGRAM_LCB join ("\t", $_, (min @{$ideogram{$_}{coords}}), (max @{$ideogram{$_}{coords}}), $ideogram{$_}{chrom}) . "\n";
    }
  } else {
    print STDERR "[WARN] LCB#$block does not seem to have two chromosomes?\n";
  }
}
## close fhs
close $IDEOGRAM_GENOME;
close $CYTOBANDS_GENOME;
close $IDEOGRAM_LCB;
close $CYTOBANDS_LCB;
close $REPEATS_LCB;

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
