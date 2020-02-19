#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Sort::Naturally;
use Getopt::Long qw(:config no_ignore_case);

## TODO


my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -i|--sam     [FILE] : BAM file of reads mapped to LTRs (BAM/SAM) [*]
  -d|--dbs   [STRING] : comma delim list of databases to search (BAM/SAM) [*]
  -n|--names [STRING] : comma delim list of column names to print, instead of full db names [inherits from '--db'] ! must be same order as '--dbs'
  -o|--out   [STRING] : output prefix for outfiles ['results_mapTEags']
  -g|--overhang [INT] : require at least this number of bases as genome tag [50]
  -m|--match    [INT] : require at least this number of bases either side of TE boundary [50]
  -e|--evalue   [STR] : BLASTn evalue [1e-20]
  -t|--threads  [INT] : number of threads for multicore operations [4]
  -q|--qcovhsp        : use blast qcovhsp value for scoring
  -c|--collapse       : collapse marginal scores to zero, rather than qcovhsp value (recommended if -q)
  -a|--mark     [STR] : mark focal individual with special symbol in output [asterisks '**']
  -y|--eyeball        : print delimiter between LTR id's for easier reading
  -h|--help           : this message
  -v|--verbose        : verbose mode
  -d|--debug          : debug mode
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
