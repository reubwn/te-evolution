#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use File::Basename;
use Sort::Naturally;
use Getopt::Long qw(:config no_ignore_case);

## TODO
## Currently this works from a BAM file with the internally mapped reads already stripped away
## i.e. only reads overlapping the ends of the LTR are retained
## It would be better to do this within the script, so you can just feed it a unprocessed BAM
## and the relevant reads are pulled straight from that

## furthermore, it currently just blasts all of the reads vs the assemblies...
## better would be to make it non-redundant - ie. pick a read spanning >50 (-g) into
## genome region; chuck any others that are identical across the same region (SOMEHOW?)

my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -i|--in      [PATH] : path to 'sequence/' dir output from 'create_tags.pl' [*]
  -d|--dbs   [STRING] : comma delim list of genomes to search (fasta) [*]
  -n|--names [STRING] : comma delim list of column names to print, instead of full db names [inherits from '--db'] ! must be same order as '--dbs'
  -o|--out   [STRING] : output prefix for outfiles ['results_blastTEags']
  -e|--evalue   [STR] : BLASTn evalue [1e-5]
  -t|--threads  [INT] : number of threads for multicore operations [4]
  -q|--qcovhsp        : use blast qcovhsp value for scoring
  -c|--collapse       : collapse marginal scores to zero, rather than qcovhsp value (recommended if -q)
  -a|--mark     [STR] : mark focal individual with special symbol in output [asterisks '**']
  -y|--eyeball        : print delimiter between LTR id's for easier reading
  -h|--help           : this message
  -d|--debug          : debug mode
\n";

## input
my (
  $IN_path, $db_string, $names_string, $use_qcovhsp_as_score, $collapse_marginal_scores, $mark, $eyeball,
  $help, $verbose, $debug
  );
## defaults
my $OUT_prefix = "results_blastTEags";
my $evalue = "1e-5";
my $threads = 4;

GetOptions (
  'i|in=s' => \$IN_path,
  'd|db=s' => \$db_string,
  'n|names:s' => \$names_string,
  'o|out:s' => \$OUT_prefix,
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
die $usage unless ( $IN_path && $db_string );
$IN_path =~ s/\/$//;

print STDERR "[####] TE-EVOLUTION blast_tags.pl\n";
print STDERR "[####] " . `date`;
##
my ( %ltr_hash, %top_hits, %final_results );

## check system for required programs
check_progs();
## get dbs to blast against
my @databases_makedb = split( m/\,/, $db_string );
my @databases_blastdb = split( m/\,/, $db_string );
my %databases_names;
if ( $names_string ) { ## get names subs if exists, otherwise just use input db names (makes for nicer table)
  @databases_names{@databases_blastdb} = split( m/\,/, $names_string ); ##key= full db name; val= sub name
} else {
  @databases_names{@databases_blastdb} = split( m/\,/, $db_string ); ##key= full db name; val= full db name
}

## shout about scoring system
if ( $use_qcovhsp_as_score ) {
  print STDERR "[INFO] Scoring system set to 'qcovhsp'\n";
} else {
  print STDERR "[INFO] Scoring system set to '1, 0.8, 0.2, 0'\n";
}

print STDERR "[INFO] Check/make blastdb's for: @databases_makedb\n";
@databases_makedb = @{ check_blastdbs(\@databases_makedb) }; ## check which database files need makeblastdb run on them
make_blastdbs( \@databases_makedb ); ## and then do it

## open BLAST results file
my $blast_file = $OUT_prefix . "_blast.txt";
open (my $BLAST, ">$blast_file") or die $!;
print $BLAST join ("\t", "qacc","sacc","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovhsp","repeat_id","ltr_pos","db","descr","result","score") . "\n";

## glob tag fasta files
my @fasta_files = glob ("$IN_path/*fasta $IN_path/*fnaa $IN_path/*fa");
print STDERR "[INFO] There are ".scalar(@fasta_files)." files in '$IN_path'\n";

## iterate across files and dbs; do BLAST
foreach my $fasta_file (nsort @fasta_files) {

  ## load $ltr_hash
  (my $repeat_id = $fasta_file) =~ s{^.*/|\.[^.]+$}{}g;
  open (my $GREP, "grep '>' $fasta_file |") or die $!;
  while (<$GREP>) {
    chomp;
    (my $ltr_id = $_) =~ s/\>//;
    $ltr_hash{$repeat_id}{L} = $ltr_id if ($ltr_id =~ m/:L:/);
    $ltr_hash{$repeat_id}{R} = $ltr_id if ($ltr_id =~ m/:R:/);
  }

  ## BLAST vs each database in turn
  foreach my $database (@databases_blastdb) {
    print STDERR "[INFO] BLASTing query '$repeat_id' TE-tags versus database '$databases_names{$database}'...\n";

    ## open blast filehandle
    open (my $BLASTCMD, "blastn -task blastn -num_threads $threads -evalue $evalue -query $fasta_file -db $database -outfmt '6 std qcovhsp' |") or die $!;

    ## iterate thru blast results
    LINE: while (my $line = <$BLASTCMD>) {
      chomp $line;
      my ($qacc, $sacc, $pident, $length, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore, $qcovhsp) = split( m/\s+/, $line );

      my $score;
      my $ltr_pos = $qacc =~ m/:L:/ ? "L" : "R";
      if ( $qcovhsp == 100 ) { ## successful BLAST alignment across entire tag
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 1;
        ## annotate BLAST result
        print $BLAST join ("\t", $line,$repeat_id,$databases_names{$database},$ltr_pos,"full","PASS",$score) . "\n";
        next LINE;

      } elsif ( $qcovhsp > 80 ) { ## successful match spanning at least 30 bp over TE/genome boundary
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.8;
        ## annotate BLAST result
        print $BLAST join ("\t", $line,$repeat_id,$databases_names{$database},$ltr_pos,"partial","PASS",$score) . "\n";
        next LINE;

      } elsif ( $qcovhsp > 50 ) { ## marginal match spanning at least 1 bp over TE/genome boundary
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.5;
        $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
        ## annotate BLAST result
        print $BLAST join ("\t", $line,$repeat_id,$databases_names{$database},$ltr_pos,"partial","MARGINAL",$score) . "\n";
        next LINE;

      } else { ## match that does not span TE/genome boundary by any overlap
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0;
        $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
        ## annotate BLAST result
        print $BLAST join ("\t", $line,$repeat_id,$databases_names{$database},$ltr_pos,"partial","FAIL",$score) . "\n";
        next LINE;
      }
    }
    close $BLASTCMD;
  }
}
close $BLAST;

## process annotated blast results
## want to save the 'best' score per query-subject pair only
print STDERR "[INFO] Parsing BLAST results file '$blast_file'...\n";
open (my $ANNOT_BLAST, "sort -k14,14V -k15,15V -k16,16V -k19,19n $blast_file |") or die $!;
while (my $line = <$ANNOT_BLAST >) {
  chomp $line;
  my @F = split ( m/\t/, $line );
  if ( scalar(@F) == 19 ) {
    $top_hits{$F[13]}{$F[14]}{$F[15]}{qcovhsp} = $F[12];
    $top_hits{$F[13]}{$F[14]}{$F[15]}{mismatches} = $F[4];
  }
}
close $ANNOT_BLAST;

print STDERR Dumper(\%top_hits);

## print condensed results
## to show presence / absence of tags across all databases
my $table_file = $OUT_prefix . "_table.txt";
open (my $TAB, ">$table_file") or die $!;
if ( $mark ) {
  print $TAB "repeat_id\tposition\tqacc\t";
  foreach ( nsort values %databases_names ) {
    if ( $_ eq $mark ) { ## convoluted way of printing a couple of asterisks, but heyho
      print $TAB "**$_\t";
    } else {
      print $TAB "$_\t";
    }
  }
  print $TAB "\n";
} else {
  print $TAB join ( "\t", "repeat_id", nsort values %databases_names ) . "\n"; ## print header
}

foreach my $repeat_id ( nsort keys %ltr_hash ) {
  ## iterate thru ALL ltrs, not just the ones with blast hits
  print $TAB $repeat_id; ## print ltr_id
  foreach my $database ( nsort values %databases_names ) {
    my $final_score;
    if ( ($top_hits{$repeat_id}{$database}{L}{qcovhsp}) && ($top_hits{$repeat_id}{$database}{R}{qcovhsp}) ) { ## if hit exists for both L && R
      $final_score = (($top_hits{$repeat_id}{$database}{L}{qcovhsp} - $top_hits{$repeat_id}{$database}{L}{mismatches}) + ($top_hits{$repeat_id}{$database}{R}{qcovhsp} - $top_hits{$repeat_id}{$database}{R}{mismatches})) / 200;
      print STDERR "Ever get in here?\n";
    } elsif ( ($top_hits{$repeat_id}{$database}{L}{qcovhsp}) && !($top_hits{$repeat_id}{$database}{R}{qcovhsp}) ) {
      $final_score = ($top_hits{$repeat_id}{$database}{L}{qcovhsp} - $top_hits{$repeat_id}{$database}{L}{mismatches}) / 200;
      # print $TAB "\t$final_score";
    } elsif ( !($top_hits{$repeat_id}{$database}{L}{qcovhsp}) && ($top_hits{$repeat_id}{$database}{R}{qcovhsp}) ) {
      $final_score = ($top_hits{$repeat_id}{$database}{R}{qcovhsp} - $top_hits{$repeat_id}{$database}{R}{mismatches}) / 200;
      # print $TAB "\t$final_score";
    } else { ## else print 0
      $final_score = 0;
      # print $TAB "\t0";
    }
    print $TAB "\t$final_score";
  }
  print $TAB "\n";
  print $TAB "###\n" if ( $eyeball );

  # ## ~~~~~~~~~~~~~
  # ## then righties
  # ## ~~~~~~~~~~~~~
  # foreach my $name ( nsort @{ $ltr_hash{$ltr}{right_names} } ) {
  #   print $TAB "$ltr\tR\t$name"; ## print ltr_id and position factors
  #   foreach my $db ( nsort values %databases_names ) {
  #     if ( $results{$name}{$db} ) { ## if exists, print score
  #       print $TAB "\t$results{$name}{$db}";
  #     } else { ## else print 0
  #       print $TAB "\t0";
  #     }
  #   }
  #   print $TAB "\n";
  # }
  # print $TAB "###\n" if ( $eyeball );
}
close $TAB;

print STDERR "[####] Done!\n";
print STDERR "[####] " . `date`;

####################
#################### SUBS
####################

sub check_progs {
  chomp( my $makeblastdb_path = `which makeblastdb` );
  chomp( my $blastn_path = `which blastn` );
  chomp( my $samtools_path = `which samtools` );
  if (!( $makeblastdb_path )) {
    die "[ERROR] Cannot find makeblastdb in \$PATH\n";
  } else {
    print STDERR "[INFO] Found makeblastdb at $makeblastdb_path\n";
  }
  if (!( $blastn_path )) {
    die "[ERROR] Cannot find blastn in \$PATH\n";
  } else {
    print STDERR "[INFO] Found blastn at $blastn_path\n";
  }
  if (!( $samtools_path )) {
    die "[ERROR] Cannot find samtools in \$PATH\n";
  } else {
    print STDERR "[INFO] Found samtools at $samtools_path\n"; ##
  }
}

sub parallel {
  chomp( my $parallel_path = `which parallel` );
  my $result = 1;
  if ( $parallel_path ) {
    $result = 0;
  }
  return $result;
}

sub check_blastdbs {
  my @in = @{ $_[0] };
  my @out;
  for my $i ( 0 .. $#in ) {
    ## file not exist
    if (! -f $in[$i]) {
      die "[ERROR] File $in[$i] does not exist! $!\n\n";
    }
    ## file is gzipped
    if ($in[$i] =~ m/gz$/) {
      die "[ERROR] Please gunzip your fasta file: $in[$i]\n\n";
    }
    ## blastdb already made, remove from @out
    if ( (-f "$in[$i].nhr") && (-f "$in[$i].nin") && (-f "$in[$i].nsq") ) {
      print STDERR "[INFO] BlastDB already exists for $in[$i]\n";
      # splice( @out, $i, 1 );
    } else {
      push ( @out, $in[$i] );
    }
  }
  return \@out;
}

sub make_blastdbs {
  my @databases = @{ $_[0] };

  ## exit if nothing to do
  return if ( scalar(@databases) == 0 );

  if ( &parallel == 0 ) { ## run in parallel
    if ( system("parallel -j $threads 'makeblastdb -in {} -dbtype nucl' ::: @databases") != 0 ) {
      die "\n[ERROR] Something wrong with makeblastdb command: $!\n\n";
    }
  } else { ## or sequentially
    foreach (@databases) {
      if ( system("makeblastdb -in $_ -dbtype nucl") != 0 ) {
        die "\n[ERROR] Something wrong with makeblastdb command: $!\n\n";
      }
    }
  }
}
