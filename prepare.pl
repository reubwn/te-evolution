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
  $help, $keep, $debug
);
## defaults
my $ks_threshold = 0.2;
my $outprefix = "prepare";
my $find = "Family";
my $numticks = 5;

GetOptions (
  'c|collinearity=s' => \$collinearity_infile,
  's|score=s' => \$score_infile,
  'g|genes=s' => \$genes_infile,
  'f|fasta=s' => \$genome_infile,
  '1|tes1=s' => \$tes_infile1,
  '2|tes2=s' => \$tes_infile2,
  'N|Ns:s' => \$nnns_infile,
  'v|coverage:s' => \$coverage_infile,
  'k|ks:f' => \$ks_threshold,
  'd|find:s' => \$find,
  'b|numticks:i' => \$numticks,
  'o|outprefix:s' => \$outprefix,
  'a|keep' => \$keep,
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
print STDERR "[INFO] Files written to: $results_dir/\n";

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
open (my $IDEOGRAM_GENOME, ">$results_dir/$outprefix"."_"."genome.ideogram") or die $!;
open (my $CYTOBANDS_GENOME, ">$results_dir/$outprefix"."_"."genome.cytobands") or die $!;
print $CYTOBANDS_GENOME join ("\t", "chr", "start", "end", "name", "gieStain") . "\n";
print $IDEOGRAM_GENOME join ("\t", "chr", "start", "end") . "\n";
foreach (nsort keys %genome_hash) {
  print $IDEOGRAM_GENOME join ("\t", $_, "1", $genome_hash{$_}) . "\n";
}
close $IDEOGRAM_GENOME;

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
print $R "par(mfrow=c(1,1))\n";
print $R "par(family=\"Ubuntu Light\",ps=12, las=1)\n";
print $R "par(mar=c(2,2,0,2), oma=c(4,4,0,4))\n";
print $R "par(tcl=-0.25)\n";
print $R "par(mgp=c(2, 0.6, 0))\n";
print $R "cols<-viridis(5,alpha=1)\n\n";

## MAIN BLOCK
## print consecutively each chr in LCB
## process $collinearity_hash
foreach my $block (sort {$a<=>$b} keys %collinearity_hash) {
  print STDERR "\r[INFO] Block number: $block ($scores_hash{$block}{'orientation'})"; $|=1;

  ## open cytobands and ideogram LCBs file
  open (my $IDEOGRAM_LCB, ">$results_dir/$outprefix"."_"."LCB\#$block.ideogram") or die $!;
  open (my $CYTOBANDS_LCB, ">$results_dir/$outprefix"."_"."LCB\#$block.cytobands") or die $!;
  open (my $LINKS_S_LCB, ">$results_dir/$outprefix"."_"."LCB\#$block.starts") or die $!;
  open (my $LINKS_E_LCB, ">$results_dir/$outprefix"."_"."LCB\#$block.ends") or die $!;
  open (my $REPEATS_LCB, ">$results_dir/$outprefix"."_"."LCB\#$block.repeats") or die $!;
  print $IDEOGRAM_LCB join ("\t", "chr", "start", "end", "name") . "\n";
  print $CYTOBANDS_LCB join ("\t", "chr", "start", "end", "name", "gieStain") . "\n";
  print $LINKS_S_LCB join ("\t", "chr", "start", "end", "strand", "name") . "\n";
  print $LINKS_E_LCB join ("\t", "chr", "start", "end", "strand", "name") . "\n";
  print $REPEATS_LCB join ("\t", "chr", "start", "end", "strand", "name") . "\n";

  ## first do for 'genes1'
  ## =====================
  ## get start and end genes in LCB array
  my $start1 = ${ $collinearity_hash{$block}{'genes1'} }[0];
  my $end1 = ${ $collinearity_hash{$block}{'genes1'} }[-1];
  my (@cytobands1_array, @linkstarts_array);
  ## debug
  print STDERR "\n[DEBUG] LCB#$block:1 ::: $start1-$end1\n" if $debug;

  ## slice from MCScanX genes file to get coordinates
  my %ideogram;
  open (my $TMP1, "perl -e 'while (<>) {print if (/\Q$start1\E/../\Q$end1\E/)}' $genes_infile |") or die $!;
  while (<$TMP1>) {
    my @F = split (m/\s+/, $_);
    my ($cytoband, $link_start);
    # if ($. == 1) { ## this prints arbitrarily large blank cytoband for each block for prettier visualisation
    #   print $CYTOBANDS_LCB join ("\t", "LCB#$block:1", -1e+9, 1e+9, "background", "gneg") . "\n";
    # }
    ## print genes to cytobands file
    if ($genes_hash{$F[1]}) { ##gene is part of LCB
      $cytoband = join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "stalk");
      $link_start = join ("\t", "LCB#$block:1", $F[2], $F[3], "+", $F[1]);
      # print $CYTOBANDS_LCB join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "stalk") . "\n";
      # print $LINKS_S_LCB join ("\t", "LCB#$block:1", $F[2], $F[3], "+", $F[1]) . "\n";
    } else {
      $cytoband = join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "gpos25");
      # print $CYTOBANDS_LCB join ("\t", "LCB#$block:1", $F[2], $F[3], $F[1], "gpos25") . "\n";
    }
    push (@cytobands1_array, $cytoband);
    push (@linkstarts_array, $link_start) if ($link_start); ##skip non-linking genes
    ## get all coords of all genes in region
    $ideogram{"LCB#$block:1"}{chrom} = $F[0]; ##key=LCB##:1; val=chrom
    push ( @{ $ideogram{"LCB#$block:1"}{coords} }, $F[2] ); ##key=LCB##:1; val=@{array of start-end coordinates}
    push ( @{ $ideogram{"LCB#$block:1"}{coords} }, $F[3] );
  }
  close $TMP1;

  ## get min max and range of coords
  my ($chrom1, $min1, $max1) = ($ideogram{"LCB#$block:1"}{chrom}, (min @{ $ideogram{"LCB#$block:1"}{coords} }), (max @{ $ideogram{"LCB#$block:1"}{coords} }) );
  my $range1 = ($max1 - $min1);

  ## print cytobands and links files
  print $CYTOBANDS_LCB join ("\t", "LCB#$block:1", $min1, $max1, "background", "gneg") . "\n"; ## print background blank band
  print $CYTOBANDS_LCB join ("\n", @cytobands1_array);
  print $LINKS_S_LCB join ("\n", @linkstarts_array);

  ## get TEs that intersect with LCB region using bedtools
  print STDERR "[DEBUG] Bedtools command: printf '$chrom1\t$min1\t$max1' | bedtools intersect -a $tes_infile1 -b stdin -wa |\n" if $debug;
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
    open ($TMP2, "perl -e 'while (<>) {print if (/\Q$start2\E/../\Q$end2\E/)}' $genes_infile |") or die $!;
    print STDERR "[DEBUG] LCB#$block:2 ::: $start2-$end2\n" if $debug;
  } else {
    ## switcheroo
    ## note this collects gene order as they appear in GFF
    open ($TMP2, "perl -e 'while (<>) {print if (/\Q$end2\E/../\Q$start2\E/)}' $genes_infile |") or die $!;
    print STDERR "[DEBUG] LCB#$block:2 ::: $end2-$start2\n" if $debug;
  }
  my (@cytobands2_array, @linkends_array); ##
  while (<$TMP2>) {
    my @F = split (m/\s+/, $_);
    # if ($. == 1) {
    #   print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", -1e+9, 1e+9, "background", "gneg") . "\n";
    # }
    my ($cytoband, $link_end);
    if ($genes_hash{$F[1]}) { ##gene is part of LBC
      $cytoband = join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "stalk");
      $link_end = join ("\t", "LCB#$block:2", $F[2], $F[3], "+", $F[1]); ##the strand '+' here is a bit arbitrary...
      # print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "stalk") . "\n";
      # print $LINKS_E_LCB join ("\t", "LCB#$block:2", $F[2], $F[3], "+", $F[1]) . "\n"; ##the strand '+' here is a bit arbitrary...
    } else { ##gene is not part of LBC
      $cytoband = join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "gpos25");
      # print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", $F[2], $F[3], $F[1], "gpos25") . "\n";
    }
    push (@cytobands2_array, $cytoband);
    push (@linkends_array, $link_end) if ($link_end); ##skip non-linking genes
    $ideogram{"LCB#$block:2"}{chrom} = $F[0]; ##
    push ( @{ $ideogram{"LCB#$block:2"}{coords} }, $F[2] );
    push ( @{ $ideogram{"LCB#$block:2"}{coords} }, $F[3] );
  }
  close $TMP2;
  ## calculate min, max coords and range
  my ($chrom2, $min2, $max2) = ($ideogram{"LCB#$block:2"}{chrom}, (min @{ $ideogram{"LCB#$block:2"}{coords} }), (max @{ $ideogram{"LCB#$block:2"}{coords} }) );
  my $range2 = ($max2 - $min2);

  ## print cytobands and links file depending on LCB orientation
  if ($scores_hash{$block}{'orientation'} eq "plus") {
    print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", $min2, $max2, "background", "gneg") . "\n"; ## print background blank band
    print $CYTOBANDS_LCB join ("\n", @cytobands2_array);
    print $LINKS_E_LCB join ("\n", @linkends_array);
  } else { ## reverse order so that LCB is always visualised in FF orientation
    @cytobands2_array = reverse @cytobands2_array;
    print $CYTOBANDS_LCB join ("\t", "LCB#$block:2", $min2, $max2, "background", "gneg") . "\n"; ## print background blank band
    print $CYTOBANDS_LCB join ("\n", @cytobands2_array);
    @linkends_array = reverse @linkends_array;
    print $LINKS_E_LCB join ("\n", @linkends_array);
  }

  ## search for intersection with TE file and pull out repeats that occur in LCB region
  print STDERR "[DEBUG] Bedtools command: printf '$chrom2\t$min2\t$max2' | bedtools intersect -a $tes_infile2 -b stdin -wa |\n" if $debug;
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
      ## debug:
      if ($debug) {
        print STDERR "[DEBUG] $_ ::: $ideogram{$_}{chrom} ".(min @{$ideogram{$_}{coords}})."-".(max @{$ideogram{$_}{coords}})."\n";
      }
    }
  } else {
    print STDERR "[WARN] LCB#$block does not seem to have two chromosomes?\n";
  }
  ## close FHs
  close $IDEOGRAM_LCB;
  close $CYTOBANDS_LCB;
  close $REPEATS_LCB;
  close $LINKS_S_LCB;
  close $LINKS_E_LCB;

  ## print commands to R file
  print $R "\n## load data for LCB\#$block\n";
  print $R "genome<-toGRanges('$results_dir/$outprefix\_LCB\#$block.ideogram')\n";
  print $R "cytobands<-toGRanges('$results_dir/$outprefix\_LCB\#$block.cytobands')\n";
  print $R "repeats<-toGRanges('$results_dir/$outprefix\_LCB\#$block.repeats')\n";
  print $R "starts<-toGRanges('$results_dir/$outprefix\_LCB\#$block.starts')\n";
  print $R "ends<-toGRanges('$results_dir/$outprefix\_LCB\#$block.ends')\n";
  print $R "## plot data for LCB\#$block\n";
  my $strand = $scores_hash{$block}{'orientation'};
  print $R "kp <- plotKaryotype(genome=genome, cytobands=cytobands, plot.type=2, main=paste(genome\$name[1],\" / \",genome\$name[2],\" \($strand\)\",sep=\"\"))\n";
  print $R "kpPlotLinks(kp, data=starts, data2=ends, col='grey95', border=NA, data.panel=1, y=-0.36)\n";
  print $R "kpAddCytobands(kp)\n";
  my $tick_dist = ($range1+$range2/2) / $numticks; ## calculate appropriate inter-tickmark distance
  print $R "kpAddBaseNumbers(kp,tick.dist=$tick_dist, add.units=T, cex = 0.8)\n";
  print $R "mtext(genome\$name[[1]], side=2, outer=T, at=0.56, adj=0, cex=0.75)\n";
  print $R "mtext(genome\$name[[2]], side=2, outer=T, at=0.2, adj=0, cex=0.75)\n";
  print $R "kpPlotLinks(kp, data=starts, data2=ends, col='grey95')\n";
  print $R "kpPlotRegions(kp, data=repeats, r0=0, r1=0.5, col=cols[1], border=cols[1])\n";
  print $R "kpPlotNames(kp, data=repeats, y0=0.1, y1=0.1, labels=repeats\$name,cex=0.5)\n";
}
## close FHs
close $R;
close $CYTOBANDS_GENOME;


######################## SUBS

sub commify {
  my $text = reverse $_[0];
  $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
  return scalar reverse $text;
}

__END__
