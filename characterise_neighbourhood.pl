#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use Sort::Naturally;

my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -i|--in        [FILE] : list of sequence names to analyse (TXT) [*]
  -s|--scaffolds [FILE] : scaffolds encoding sequences (FASTA) [*]
  -g|--gff       [FILE] : GFF gene location file (GFF) [*]
  -a|--annot     [FILE] : TSV gene annotation file from InterProScan (TSV)
  -t|--telo    [STRING] : telomeric repeat string to search for (default TGTGGG)
  -x|--telox      [INT] : number of repeat units to identify before classifying as a match (default 2)
  -z|--mismatch   [INT] : number of mismatches in --telo to allow (default 1)
  -t|--type    [STRING] : define name/type of sequence e.g. LTR, PLE, etc.
  -h|--help             : this message
\n";

## input
my (
  $sam_infile, $db_string, $names_string, $use_qcovhsp_as_score, $collapse_marginal_scores, $mark, $eyeball,
  $help, $verbose, $debug
  );
## defaults
my $outprefix = "results_mapTEags";
my $overhang_threshold = 50;
my $match_threshold = 50;
my $evalue = "1e-20";
my $threads = 4;

GetOptions (
  'i|sam=s' => \$sam_infile,
  'd|db=s' => \$db_string,
  'n|names:s' => \$names_string,
  'o|out:s' => \$outprefix,
  'g|overhang:i' => \$overhang_threshold,
  'm|match:i' => \$match_threshold,
  'e|evalue:s' => \$evalue,
  't|threads:i' => \$threads,
  'q|qcovhsp' => \$use_qcovhsp_as_score,
  'c|collapse' => \$collapse_marginal_scores,
  'a|mark:s' => \$mark,
  'y|eyeball' => \$eyeball,
  'h|help' => \$help,
  'v|verbose' => \$verbose,
  'debug' => \$debug
  );
## help and usage
die $usage if $help;
die $usage unless ( $sam_infile && $db_string );

print STDERR "[####] TE-EVOLUTION blast_tags.pl\n";
print STDERR "[####] " . `date`;
