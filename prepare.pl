#!/usr/bin/env perl

##
## prepare.pl
## reubwn 2018
##

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case);;
use Sort::Naturally;
use List::Util qw( min max );
use Data::Dumper;

use Bio::Seq;
use Bio::SeqIO;

my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -o|--out          [STRING] : Prefix for outfiles
  -c|--collinearity [FILE]   : MCScanX collinearity file (reformatted!) [*]
  -s|--score        [FILE]   : MCScanX score file [*]
  -a|--annot        [FILE]   : MCScanX genes annotation file (GFF) [*]
  -t|--te1          [FILE]   : TE annotation file (species1) (GFF3) [*]
  -T|--te2          [FILE]   : TE annotation file (species2) (GFF3) [*]
  -g|--fasta1       [FILE]   : Genome fasta file (species1)
  -G|--fasta2       [FILE]   : Genome fasta file (species2)
  -b|--bam1         [FILE]   : Coverage file [BAM] (species1)
  -B|--bam2         [FILE]   : Coverage file [BAM] (species2)
  -p|--split1       [FILE]   : Split reads coverage file [BAM] (species1)
  -P|--split2       [FILE]   : Split reads coverage file [BAM] (species2)
  -d|--disc1        [FILE]   : Discordant reads coverage file [BAM] (species1)
  -D|--disc2        [FILE]   : Discordant reads coverage file [BAM] (species2)
  -k|--ks           [FLOAT]  : Ks threshold to filter LCBs [0.2]
  -N|--gaps         [INT]    : Size of scaffold gap to annotate (>=10 Ns)
  -f|--find         [STRING] : Search string to grep TE id from GFF file ['Class']
  -l|--nolegend              : Don't plot legend [FALSE]
  -n|--names                 : Plot repeat names [FALSE]
  -h|--help                  : this message

OUTPUTS
\n";

## input
my (
  $collinearity_infile,
  $score_infile,
  $genes_infile,
  $genome1_infile,
  $genome2_infile,
  $tes1_infile,
  $tes2_infile,
  $coverage1_infile,
  $coverage2_infile,
  $split1_infile,
  $split2_infile,
  $disc1_infile,
  $disc2_infile,
  $plot_legend_off,
  $plot_names,
  $help, $keep, $debug
);
## defaults
my $ks_threshold = 0.2;
my $gaps_threshold = 10;
my $outprefix = "prepare";
my $find = "Class"; ## this will be grepped from the TE annotation file like: if (m/$find\=([\w]+)(\;.+)*/) { ... }
my $numticks = 5;

GetOptions (
  'c|collinearity=s' => \$collinearity_infile,
  's|score=s' => \$score_infile,
  'a|annot=s' => \$genes_infile,
  't|te1=s' => \$tes1_infile,
  'T|te2=s' => \$tes2_infile,
  'g|fasta1:s' => \$genome1_infile,
  'G|fasta2:s' => \$genome2_infile,
  'b|bam1:s' => \$coverage1_infile,
  'B|bam2:s' => \$coverage2_infile,
  'p|split1:s' => \$split1_infile,
  'P|split2:s' => \$split2_infile,
  'd|disc1:s' => \$disc1_infile,
  'D|disc2:s' => \$disc2_infile,
  'k|ks:f' => \$ks_threshold,
  'N|gaps:i' => \$gaps_threshold,
  'f|find:s' => \$find,
  'o|outprefix:s' => \$outprefix,
  'l|legend' => \$plot_legend_off,
  'n|names' => \$plot_names,
  'm|numticks:i' => \$numticks,
  'j|keep' => \$keep,
  'h|help' => \$help,
  'debug' => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ($collinearity_infile && $score_infile && $genes_infile && $tes1_infile && $tes2_infile);
print STDERR "[####] TE-EVOLUTION prepare.pl\n";
print STDERR "[####] " . `date`;

## look for bedtools
if ( can_run('bedtools') ) {
  print STDERR "[INFO] Found bedtools! (".can_run('bedtools').")\n";
} else {
  die "[ERROR] Bedtools is not installed or not discoverable in \$PATH\n\n";
}

## look for seqtk
if ($genome1_infile || $genome2_infile) {
  if ( can_run('seqtk') ) {
    print STDERR "[INFO] Found seqtk! (".can_run('seqtk').")\n";
  } else {
    die "[ERROR] Seqtk is not installed or not discoverable in \$PATH\n\n";
  }
}

## stuff
my ( %scores_hash, %collinearity_hash, %genes_hash, %genome_hash, %tes_hash, %nnns_hash, %coverage_hash );

## open outdir
my $results_dir = $outprefix."_data";
if (-d $results_dir) {
  if ($keep) {
    die "[ERROR] Dir $results_dir already exists!\n";
  } else {
    system ("rm -r $results_dir && mkdir $results_dir");
  }
} else {
  system ("mkdir $results_dir");
}
print STDERR "[INFO] Prefix is '$outprefix'\n";
print STDERR "[INFO] Writing output to: '$results_dir/'\n";
print STDERR "[INFO] MCScanX collinearity file: $collinearity_infile\n";
print STDERR "[INFO] MCScanX scores file: $score_infile\n";

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
print STDERR "[INFO] Number of LCBs passing Ks threshold (<=$ks_threshold): ".commify(scalar(keys %scores_hash))."\n";

## parse collinearity file
open (my $COLL, $collinearity_infile) or die $!;
while (<$COLL>) {
  if ($_ =~ m/^#/) {
    next;
  } else {
    my @F = split (m/\s+/, $_);
    if ($scores_hash{$F[0]}) { ##analyse blocks with Ks <= $ks_threshold
      push ( @{ $genes_hash{$F[2]} }, $F[0] ); ##note genes that form LCB; key=gene; val=@{blocks}
      push ( @{ $genes_hash{$F[3]} }, $F[0] );
      push ( @{ $collinearity_hash{$F[0]}{'genes1'} }, $F[2]); ##anon array of member genes from species1
      push ( @{ $collinearity_hash{$F[0]}{'genes2'} }, $F[3]);
    } else {
      next;
    }
  }
}
close $COLL;
print STDERR "[INFO] Number of genes in LCBs: ".commify(scalar(keys %genes_hash))."\n";

## make Ns annotation files if genomes are provided
if ($genome1_infile) {
  if (system("seqtk cutN -gp1000000 -n10 $genome1_infile > $genome1_infile.$outprefix.gaps_annot.txt") != 0) {
    die "[ERROR] Seqtk command ran with some errors\n";
  } else {
    print STDERR "[INFO] Seqtk: annotated gaps (Ns) in $genome1_infile\n";
  }
}
if ($genome2_infile) {
  if (system("seqtk cutN -gp1000000 -n10 $genome2_infile > $genome2_infile.$outprefix.gaps_annot.txt") != 0) {
    die "[ERROR] Seqtk command ran with some errors\n";
  } else {
    print STDERR "[INFO] Seqtk: annotated gaps (Ns) in $genome2_infile\n";
  }
}

## make R file for easy plotting
open (my $R, ">$outprefix.R") or die $!;
my $date = `date`;
print $R "## $date\n";
print $R "## libraries\n";
print $R "library(karyoploteR)\n";
print $R "library(viridis)\n";
print $R "## path to working dir\n";
print $R "setwd(getwd())\n";
print $R "## graphics\n";
# print $R "par(mfrow=c(1,1))\n";
print $R "par(family=\"Ubuntu Light\",ps=12, las=1)\n";
print $R "par(mar=c(1,1,1,1), oma=c(1,1,1,10))\n";
print $R "par(tcl=-0.25)\n";
print $R "par(mgp=c(2, 0.6, 0))\n";
print $R "hists<-viridis(4, alpha=1, option='A')\n";
## define some colors for plotting
my %repeat_colors;
if ($find eq "Class") {
  %repeat_colors = ( DNA => '1', RC => '2', LINE => '3', SINE => '4', LTR => '5', Other => '6');
  print $R "cols<-c(viridis(5,alpha=1),'black')\n\n";
} elsif ($find eq "Family") {
  %repeat_colors = (); ##TODO
}

## MAIN BLOCK
## print consecutively each chr in LCB
## process $collinearity_hash
foreach my $block (sort {$a<=>$b} keys %collinearity_hash) {
  print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'})"; $|=1;

  ## open cytobands and ideogram LCBs file
  open (my $IDEOGRAM, ">$results_dir/$outprefix"."_"."LCB\#$block.ideogram") or die $!;
  open (my $CYTOBANDS, ">$results_dir/$outprefix"."_"."LCB\#$block.cytobands") or die $!;
  open (my $LINK_STARTS, ">$results_dir/$outprefix"."_"."LCB\#$block.starts") or die $!;
  open (my $LINK_ENDS, ">$results_dir/$outprefix"."_"."LCB\#$block.ends") or die $!;
  open (my $REPEATS, ">$results_dir/$outprefix"."_"."LCB\#$block.repeats") or die $!;
  open (my $COV, ">$results_dir/$outprefix"."_"."LCB\#$block.coverage") or die $!;
  open (my $SPLIT, ">$results_dir/$outprefix"."_"."LCB\#$block.split") or die $!;
  open (my $DISC, ">$results_dir/$outprefix"."_"."LCB\#$block.disc") or die $!;
  print $IDEOGRAM join ("\t", "chr", "start", "end", "name") . "\n";
  print $CYTOBANDS join ("\t", "chr", "start", "end", "name", "gieStain") . "\n";
  print $LINK_STARTS join ("\t", "chr", "start", "end", "strand", "name") . "\n";
  print $LINK_ENDS join ("\t", "chr", "start", "end", "strand", "name") . "\n";
  print $REPEATS join ("\t", "chr", "start", "end", "strand", "name") . "\n";
  # print $COV join ("\t", "chr", "start", "end", "strand", "name") . "\n";
  # print $SPLIT join ("\t", "chr", "start", "end", "strand", "name") . "\n";
  # print $DISC join ("\t", "chr", "start", "end", "strand", "name") . "\n";
  print $COV join ("\t", "chr", "x", "y") . "\n";
  print $SPLIT join ("\t", "chr", "x", "y") . "\n";
  print $DISC join ("\t", "chr", "x", "y") . "\n";

  ## first do for 'genes1'
  ## =====================
  my (@cytobands_chr1, @link_starts, @repeats_chr1, @repeatClass_chr1, @N_chr1, @repeatClass_both); ##stuff

  ## LCB GENES
  ## =========
  ## get first and last genes in LCB for chr1
  my $firstGene_chr1 = ${ $collinearity_hash{$block}{'genes1'} }[0];
  my $lastGene_chr1 = ${ $collinearity_hash{$block}{'genes1'} }[-1];
  print STDERR "\n[DEBUG] LCB#$block:1 ::: $firstGene_chr1-$lastGene_chr1\n" if $debug; ##debug
  ## slice from MCScanX genes file to get coordinates
  my %ideogram;
  open (my $TMP1, "perl -e 'while (<>) {print if (/\Q$firstGene_chr1\E/../\Q$lastGene_chr1\E/)}' $genes_infile |") or die $!;
  while (<$TMP1>) {
    my @F = split (m/\s+/, $_);
    my ($cytoband, $link_start);
    ## print genes to cytobands file
    if ($genes_hash{$F[1]}) { ## gene is part of ANY LCB
      if ( grep { $block eq $_ } @{$genes_hash{$F[1]}} ) { ## gene is part of THIS SPECIFIC LCB (and also possibly others)
        $cytoband = join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "stalk");
        $link_start = join ("\t", "LCB#$block:1", $F[2], $F[3], "+", $F[1]);
      } else {
        $cytoband = join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "gpos75"); ## gene is a member of some other LCB...
      }
    } else {
      $cytoband = join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "gpos25");
    }
    push (@cytobands_chr1, $cytoband);
    push (@link_starts, $link_start) if ($link_start); ##skip non-linking genes
    ## get all coords of all genes in region
    $ideogram{"LCB#$block:1"}{chrom} = $F[0]; ##key=LCB##:1; val=chrom
    push ( @{ $ideogram{"LCB#$block:1"}{coords} }, $F[2] ); ##key=LCB##:1; val=@{array of start-end coordinates}
    push ( @{ $ideogram{"LCB#$block:1"}{coords} }, $F[3] );
  }
  close $TMP1;
  ## get min max and range of coords
  my ($chrom_chr1, $min_chr1, $max_chr1) = ($ideogram{"LCB#$block:1"}{chrom}, (min @{ $ideogram{"LCB#$block:1"}{coords} }), (max @{ $ideogram{"LCB#$block:1"}{coords} }) );
  my $range_chr1 = ($max_chr1 - $min_chr1);
  ## print cytobands and links files
  print $CYTOBANDS join ("\t", "LCB#$block:1", $min_chr1, $max_chr1, "background", "gneg") . "\n"; ## print background blank band
  print $CYTOBANDS join ("\n", @cytobands_chr1) . "\n";
  print $LINK_STARTS join ("\n", @link_starts) . "\n";

  ## REPEATS
  ## =======
  ## get TEs that intersect with LCB region using bedtools
  print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr1: REPEATS]"; $|=1;
  print STDERR "[DEBUG] Bedtools command: printf '$chrom_chr1\t$min_chr1\t$max_chr1' | bedtools intersect -a $tes1_infile -b stdin -wa |\n" if $debug;
  open (CMD, "printf '$chrom_chr1\t$min_chr1\t$max_chr1' | bedtools intersect -a $tes1_infile -b stdin -wa |") or die $!;
  while (<CMD>) {
    if (m/$find\=([\w]+)(\;.+)*/) {
      my @F = split (/\s+/, $_);
      push (@repeatClass_chr1, $1);
      push (@repeats_chr1, join("\t", "LCB#$block:1", $F[3], $F[4], $F[6], $1)); ##could just print it here for chr1
    }
  }
  close CMD;
  print $REPEATS join ("\n", @repeats_chr1) . "\n";

  ## GAPS
  ## ====
  ## get NNNs that intersect with LCB region using bedtools and seqtk
  if ($genome1_infile) {
    print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr1: GAPS]"; $|=1;
    open (CMD, "printf '$chrom_chr1\t$min_chr1\t$max_chr1' | bedtools intersect -a $genome1_infile.$outprefix.gaps_annot.txt -b stdin -wa |") or die $!;
    while (<CMD>) {
      my @F = split (/\s+/, $_);
      print $CYTOBANDS join ("\t", "LCB#$block:1", $F[1], $F[2], "gap", "acen") . "\n"; ## just print it directly
      # push (@N_chr1, join("\t", "LCB#$block:1", $F[1], $F[2], "gap", "acen"));
    }
    close CMD;
    # print $CYTOBANDS join ("\n", @N_chr1) . "\n";
  }

  ## COVERAGE
  ## ========
  ## get overall coverage information from CIS reads
  if ($coverage1_infile) {
    print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr1: COVERAGE]"; $|=1;
    open (CMD, "printf '$chrom_chr1\t$min_chr1\t$max_chr1' | bedtools intersect -abam $coverage1_infile -b stdin -wa | bedtools bamtobed -i stdin |") or die $!;
    while (<CMD>) {
      my @F = split (/\s+/, $_);
      print $COV join ("\t", "LCB#$block:1", $F[1], $F[2], $F[5], $F[3]) . "\n";
    }
    close CMD;
    # open (CMD, "perl -e 'while (<>) {print if (/$chrom_chr1\\s+$min_chr1\\s+\\d+\n/../$chrom_chr1\\s+$max_chr1\\s+\\d+\n/)}' $coverage1_infile |") or die $!;
    # while (<CMD>) {
    #   my @F = split (/\s+/, $_);
    #   print $COV join ("\t", "LCB#$block:1", $F[1], $F[2], $F[5], $F[3]) . "\n";
    # }
    # close CMD;
  }
  ## get coverage information from TRANS split reads (i.e. reads from species 2 mapped to species 1)
  if ($split1_infile) {
    print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr1: SPLIT]"; $|=1;
    open (CMD, "printf '$chrom_chr1\t$min_chr1\t$max_chr1' | bedtools intersect -abam $split1_infile -b stdin -wa | bedtools bamtobed -i stdin |") or die $!;
    while (<CMD>) {
      my @F = split (/\s+/, $_);
      print $SPLIT join ("\t", "LCB#$block:1", $F[1], $F[2], $F[5], $F[3]) . "\n";
    }
    close CMD;
  }
  ## get coverage information from TRANS discordant reads
  if ($disc1_infile) {
    print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr1: DISCORDANT]"; $|=1;
    open (CMD, "printf '$chrom_chr1\t$min_chr1\t$max_chr1' | bedtools intersect -abam $disc1_infile -b stdin -wa | bedtools bamtobed -i stdin |") or die $!;
    while (<CMD>) {
      my @F = split (/\s+/, $_);
      print $DISC join ("\t", "LCB#$block:1", $F[1], $F[2], $F[5], $F[3]) . "\n";
    }
    close CMD;
  }

  ## then do for 'genes2'
  ## ====================
  my (@cytobands_chr2, @link_ends, @repeats_chr2, @repeatClass_chr2, @N_chr2, @cov_chr2, @covSplit_chr2, @covDisc_chr2); ##stuff

  ## LCB GENES
  ## =========
  ## get first and last genes in LCB for chr2
  my $firstGene_chr2 = ${ $collinearity_hash{$block}{'genes2'} }[0];
  my $lastGene_chr2 = ${ $collinearity_hash{$block}{'genes2'} }[-1];
  ## check orientation of LCB on second strand!
  my $TMP2;
  if ($scores_hash{$block}{'orientation'} eq "plus") {
    open ($TMP2, "perl -e 'while (<>) {print if (/\Q$firstGene_chr2\E/../\Q$lastGene_chr2\E/)}' $genes_infile |") or die $!;
    print STDERR "[DEBUG] LCB#$block:2 ::: $firstGene_chr2-$lastGene_chr2\n" if $debug;
  } else { ##switcheroo!
    ## note this collects gene order as they appear in GFF
    open ($TMP2, "perl -e 'while (<>) {print if (/\Q$lastGene_chr2\E/../\Q$firstGene_chr2\E/)}' $genes_infile |") or die $!;
    print STDERR "[DEBUG] LCB#$block:2 ::: $lastGene_chr2-$firstGene_chr2\n" if $debug;
  }
  while (<$TMP2>) {
    my @F = split (m/\s+/, $_);
    my ($cytoband, $link_end);
    if ($genes_hash{$F[1]}) { ## gene is part of ANY LCB
      if ( grep { $block eq $_ } @{$genes_hash{$F[1]}} ) { ## gene is part of THIS SPECIFIC LCB
        $cytoband = join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "stalk");
        $link_end = join ("\t", "LCB#$block:2", $F[2], $F[3], "+", $F[1]);
      } else {
        $cytoband = join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "gpos75"); ## gene is a member of some other LCB...
      }
    } else { ##gene is not part of any LBC
      $cytoband = join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "gpos25");
    }
    push (@cytobands_chr2, $cytoband);
    push (@link_ends, $link_end) if ($link_end); ##skip non-linking genes
    $ideogram{"LCB#$block:2"}{chrom} = $F[0]; ##
    push ( @{ $ideogram{"LCB#$block:2"}{coords} }, $F[2] );
    push ( @{ $ideogram{"LCB#$block:2"}{coords} }, $F[3] );
  }
  close $TMP2;
  ## calculate min, max coords and range
  my ($chrom_chr2, $min_chr2, $max_chr2) = ($ideogram{"LCB#$block:2"}{chrom}, (min @{ $ideogram{"LCB#$block:2"}{coords} }), (max @{ $ideogram{"LCB#$block:2"}{coords} }) );
  my $range_chr2 = ($max_chr2 - $min_chr2);

  ## REPEATS
  ## =======
  ## get TEs that intersect with LCB region using bedtools
  print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr2: REPEATS]"; $|=1;
  print STDERR "[DEBUG] Bedtools command: printf '$chrom_chr2\t$min_chr2\t$max_chr2' | bedtools intersect -a $tes2_infile -b stdin -wa |\n" if $debug;
  open (my $BED2, "printf '$chrom_chr2\t$min_chr2\t$max_chr2' | bedtools intersect -a $tes2_infile -b stdin -wa |") or die $!;
  my %repeat_count;
  while (<$BED2>) {
    if (m/$find\=([\w]+)(\;.+)*/) {
      my @F = split (/\s+/, $_);
      push (@repeatClass_chr2, $1);
      push (@repeats_chr2, join("\t", "LCB#$block:2", $F[3], $F[4], $F[6], $1));
    }
  }
  close $BED2;

  ## GAPS
  ## ====
  ## get NNNs that intersect with LCB region using bedtools and seqtk
  if ($genome2_infile) {
    print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr2: GAPS]"; $|=1;
    open (CMD, "printf '$chrom_chr2\t$min_chr2\t$max_chr2' | bedtools intersect -a $genome2_infile.$outprefix.gaps_annot.txt -b stdin -wa |") or die $!;
    while (<CMD>) {
      my @F = split (/\s+/, $_);
      push (@N_chr2, join("\t", "LCB#$block:2", $F[1], $F[2], "gap", "acen"));
    }
    close CMD;
  }

  ## COVERAGE
  ## ========
  ## get overall coverage information from CIS reads
  if ($coverage2_infile) {
    print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr2: COVERAGE]"; $|=1;
    open (CMD, "printf '$chrom_chr2\t$min_chr2\t$max_chr2' | bedtools intersect -abam $coverage2_infile -b stdin -wa | bedtools bamtobed -i stdin |") or die $!;
    while (<CMD>) {
      my @F = split (/\s+/, $_);
      push (@cov_chr2, join ("\t", "LCB#$block:2", $F[1], $F[2], $F[5], $F[3]));
    }
    close CMD;
  }
  ## get coverage information from TRANS split reads (i.e. reads from species 2 mapped to species 1)
  if ($split2_infile) {
    print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr2: SPLIT]"; $|=1;
    open (CMD, "printf '$chrom_chr2\t$min_chr2\t$max_chr2' | bedtools intersect -abam $split2_infile -b stdin -wa | bedtools bamtobed -i stdin |") or die $!;
    while (<CMD>) {
      my @F = split (/\s+/, $_);
      push (@covSplit_chr2, join ("\t", "LCB#$block:2", $F[1], $F[2], $F[5], $F[3]));
    }
    close CMD;
  }
  ## get coverage information from TRANS discordant reads
  if ($disc2_infile) {
    print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'}) [$chrom_chr2: DISC]"; $|=1;
    open (CMD, "printf '$chrom_chr2\t$min_chr2\t$max_chr2' | bedtools intersect -abam $disc2_infile -b stdin -wa | bedtools bamtobed -i stdin |") or die $!;
    while (<CMD>) {
      my @F = split (/\s+/, $_);
      push (@covDisc_chr2, join ("\t", "LCB#$block:2", $F[1], $F[2], $F[5], $F[3]));
    }
    close CMD;
  }

  ## PRINT DEPENDING ON LCB ORIENTATION
  ## ==================================
  ## print cytobands, links and repeat files depending on LCB orientation
  if ($scores_hash{$block}{'orientation'} eq "plus") {
    print $CYTOBANDS join ("\t", "LCB#$block:2", $min_chr2, $max_chr2, "background", "gneg") . "\n"; ## print background blank band
    print $CYTOBANDS join ("\n", @cytobands_chr2) . "\n";
    print $CYTOBANDS join ("\n", @N_chr2) . "\n" if ($genome1_infile);
    print $LINK_ENDS join ("\n", @link_ends) . "\n";
    print $REPEATS join ("\n", @repeats_chr2) . "\n";
    print $COV join ("\n", @cov_chr2) . "\n" if ($coverage2_infile);
    print $SPLIT join ("\n", @covSplit_chr2) . "\n" if ($split2_infile);
    print $DISC join ("\n", @covDisc_chr2) . "\n" if ($disc2_infile);
    @repeatClass_both = (@repeatClass_chr1, @repeatClass_chr2);
  } else {
    ## reverse order so that LCB is always visualised in FF orientation
    print $CYTOBANDS join ("\t", "LCB#$block:2", $min_chr2, $max_chr2, "background", "gneg") . "\n"; ## print background blank band
    print $CYTOBANDS join ("\n", @{ revcomp_lcb(\@cytobands_chr2, $min_chr2, $max_chr2) }) . "\n";
    print $CYTOBANDS join ("\n", @{ revcomp_lcb(\@N_chr2, $min_chr2, $max_chr2) }) . "\n" if ($genome2_infile);
    print $LINK_ENDS join ("\n", @{ revcomp_lcb(\@link_ends, $min_chr2, $max_chr2) }) . "\n";
    print $REPEATS join ("\n", @{ revcomp_lcb(\@repeats_chr2, $min_chr2, $max_chr2) }) . "\n";
    print $COV join ("\n", @{ revcomp_lcb(\@cov_chr2, $min_chr2, $max_chr2) }) . "\n" if ($coverage2_infile);
    print $SPLIT join ("\n", @{ revcomp_lcb(\@covSplit_chr2, $min_chr2, $max_chr2) }) . "\n" if ($split2_infile);
    print $DISC join ("\n", @{ revcomp_lcb(\@covDisc_chr2, $min_chr2, $max_chr2) }) . "\n" if ($disc2_infile);
    @repeatClass_both = (@repeatClass_chr1, (reverse @repeatClass_chr2));
  }

  ## IDEOGRAM (BACKBONE)
  ## ===================
  ## print the ideogram for each LCB, checking there are 2 chroms in %ideogram
  if (scalar(keys %ideogram) == 2) {
    foreach (nsort keys %ideogram) {
      print $IDEOGRAM join ("\t", $_, (min @{$ideogram{$_}{coords}}), (max @{$ideogram{$_}{coords}}), $ideogram{$_}{chrom}) . "\n";
      print STDERR "[DEBUG] $_ ::: $ideogram{$_}{chrom} ".(min @{$ideogram{$_}{coords}})."-".(max @{$ideogram{$_}{coords}})."\n" if ($debug); ##debug
    }
  } else {
    print STDERR "[WARN] LCB#$block does not seem to have two chromosomes?\n";
  }

  ## close FHs
  close $IDEOGRAM;
  close $CYTOBANDS;
  close $REPEATS;
  close $LINK_STARTS;
  close $LINK_ENDS;
  close $COV;
  close $SPLIT;
  close $DISC;

  ## R STUFF
  ## =======
  ## print commands to R file
  print $R "{\n";
  print $R "\t## load data for LCB\#$block\n";
  print $R "\tgenome<-toGRanges('$results_dir/$outprefix\_LCB\#$block.ideogram')\n";
  print $R "\tcytobands<-toGRanges('$results_dir/$outprefix\_LCB\#$block.cytobands')\n";
  print $R "\trepeats<-toGRanges('$results_dir/$outprefix\_LCB\#$block.repeats')\n";
  print $R "\tcoverage<-toGRanges('$results_dir/$outprefix\_LCB\#$block.coverage')\n";
	print $R "\tsplit<-toGRanges('$results_dir/$outprefix\_LCB\#$block.split')\n";
	print $R "\tdisc<-toGRanges('$results_dir/$outprefix\_LCB\#$block.disc')\n";
  # print $R "\trepeats<-toGRanges('$results_dir/$outprefix\_LCB\#$block.repeats')\n";
  print $R "\tstarts<-toGRanges('$results_dir/$outprefix\_LCB\#$block.starts')\n";
  print $R "\tends<-toGRanges('$results_dir/$outprefix\_LCB\#$block.ends')\n";
  print $R "\t## plot data for LCB\#$block\n";
  my $strand = $scores_hash{$block}{'orientation'};
  print $R "\tkp <- plotKaryotype(genome=genome, cytobands=cytobands, plot.type=2, labels.plotter=NULL, main=\'LCB\#$block ($strand)\')\n"; #paste(genome\$name[1],\" / \",genome\$name[2],\" \($strand\)\",sep=\"\")
  print $R "\tkpPlotLinks(kp, data=starts, data2=ends, col='grey95', border='grey95', data.panel=1, y=-0.36)\n";
  print $R "\tkpAddCytobands(kp)\n";
  my $tick_dist = ($range_chr1+$range_chr2/2) / $numticks; ## calculate appropriate inter-tickmark distance
  print $R "\tkpAddBaseNumbers(kp,tick.dist=$tick_dist, add.units=T, cex = 0.8)\n";
  print $R "\tmtext(genome\$name[[1]], side=2, outer=T, at=0.698, adj=0, line=0, cex=1)\n";
  print $R "\tmtext(genome\$name[[2]], side=2, outer=T, at=0.286, adj=0, line=0, cex=1)\n";
  my $r_colors_string;
  if (($find eq "Class") or ($find eq "Family")) {
    my @r_colors;
    foreach (@repeatClass_both) {
      if ($repeat_colors{$_}) {
        # print STDERR "$_\n";
        push(@r_colors, "cols[$repeat_colors{$_}]");
      } else {
        push(@r_colors, "cols[$repeat_colors{Other}]");
      }
    }
    $r_colors_string = join (",", @r_colors);
  } else {
    $r_colors_string = "black";
  }
  print $R "\tkpPlotRegions(kp, data=repeats, r0=0, r1=0.1, avoid.overlapping=F, col=c($r_colors_string), border=c($r_colors_string), lwd=2)\n";
  print $R "\tkpPlotCoverage(kp, data=coverage, r0=0.2, r1=0.48, col=hists[1])\n";
	print $R "\tkpPlotCoverage(kp, data=split, r0=0.5, r1=0.78, col=hists[2])\n";
	print $R "\tkpPlotCoverage(kp, data=disc, r0=0.8, r1=1.08, col=hists[3])\n";
  print $R "\tlegend(1,0.9, title='Chromosome', c('collinear gene (this LCB)','collinear gene (different LCB)','non-collinear gene','assembly gap (\u2265 10 Ns)'), pch=15, col=c('#647fa4','#828282','#c8c8c8','#d92f27'), xjust=0, yjust=1, bg='grey95', box.col='grey50',cex=0.75, pt.cex=2, xpd=NA)\n";
  print $R "\tlegend(1,0.75, title='Repeats', c('DNA','Helitron','LINE','SINE','LTR','Other'), pch=15, col=cols, xjust=0, yjust=1, bg='grey95', box.col='grey50',cex=0.75, pt.cex=2, xpd=NA)\n";
  print $R "\tkpPlotNames(kp, data=repeats, y0=0.1, y1=0.1, labels=repeats\$name,cex=0.5)\n" if $plot_names;
  print $R "}\n\n";
}
## close FHs
close $R;

######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

sub can_run {
  my $tool_name = $_[0];
  my $tool_path = `which $tool_name`;
  $tool_path =~ s/\n//g;
  return $tool_path;
}

sub check_bam { ## this doesn't work...
  my $chrom = $_[0];
  my $bamfile = $_[1];
  my %check = ($chrom => 0);
  open (CMD, "samtools view -H $bamfile |");
  while (<CMD>) {
    my @F = split (m/\s+/, $_);
    $check{$chrom}++;
  }
  close CMD;
  if ($check{$chrom} == 0) {
    return 1;
  } else {
    return 0;
  }
}

sub revcomp_lcb {
  my @reversed = reverse @{$_[0]};
  my $lcb_start = $_[1];
  my $lcb_end = $_[2];
  my @revcomp;
  foreach (@reversed) {
    my @F = split (m/\s+/, $_);
    my ($rc_start, $rc_end) = (($lcb_end - $F[2] + $lcb_start), ($lcb_end - $F[1] + $lcb_start));
    splice (@F, 1, 1, $rc_start);
    splice (@F, 2, 1, $rc_end);
    my $new_coords = join ("\t", @F);
    push (@revcomp, $new_coords);
  }
  return \@revcomp;
}

__END__
