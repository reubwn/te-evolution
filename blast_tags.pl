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

OPTIONS [*] = required
  -i|--input   [FILE] : FoFN list of genomes to search (txt) [*]
  -f|--fasta   [PATH] : path to 'sequence/' dir output from 'create_tags.pl' [*]
  -d|--dbs   [STRING] : comma delim list of genomes to search (fasta) [alternative to -i]
  -n|--names [STRING] : comma delim list of column names to print, instead of full db names; must be same order as '--dbs'
  -o|--out   [STRING] : output prefix for outfiles ['results_blastTEags']
  -e|--evalue   [STR] : BLASTn evalue [1e-5]
  -t|--threads  [INT] : number of threads for multicore operations [4]
  -q|--qcovhsp        : use blast qcovhsp value for scoring
  -c|--collapse       : collapse marginal scores to zero, rather than qcovhsp value
  -a|--mark     [STR] : mark focal individual with special symbol in output [asterisks '**']
  -y|--eyeball        : print delimiter between LTR id's for easier reading
  -h|--help           : this message

DETAILS
  --in is a simple list of paths to fasta files, with sample ID as second column
  eg. '~/path/to/genome_foo.fasta foo'
\n";

## input
my (
  $infile, $fasta_path, $db_string, $names_string, $use_qcovhsp_as_score, $collapse_marginal_scores, $mark, $eyeball,
  $help, $verbose, $debug
  );
## defaults
my $outprefix = "results";
my $evalue = "1e-5";
my $threads = 4;

GetOptions (
  'i|input:s' => \$infile,
  'f|fasta=s' => \$fasta_path,
  'd|db:s' => \$db_string,
  'n|names:s' => \$names_string,
  'o|out:s' => \$outprefix,
  'e|evalue:s' => \$evalue,
  't|threads:i' => \$threads,
  'q|qcovhsp' => \$use_qcovhsp_as_score,
  'c|collapse' => \$collapse_marginal_scores,
  'a|mark:s' => \$mark,
  'y|eyeball' => \$eyeball,
  'h|help' => \$help,
  'verbose' => \$verbose,
  'debug' => \$debug
  );
## help and usage
die $usage if $help;
die $usage unless ( ($infile || $db_string) && $fasta_path );
$fasta_path =~ s/\/$//;

print STDERR "[####] TE-EVOLUTION blast_tags.pl\n";
print STDERR "[####] " . `date`;

## stuff
my ( %ltr_hash, %top_hits, %final_results );

## check system for required programs
check_progs();

## get dbs to blast against
my (@databases_makedb, @databases_blastdb, %databases_names);
if ($infile) {
  print STDERR "[INFO] Using genomes found in '$infile' as input\n";
  open (my $DB, $infile) or die $!;
  while (<$DB>) {
    chomp;
    my @F = split (m/\s+/);
    if (scalar(@F)==2) {
      push @databases_makedb, $F[0];
      push @databases_blastdb, $F[0];
      $databases_names{$F[0]} = $F[1];
    } else {
      push @databases_makedb, $F[0];
      push @databases_blastdb, $F[0];
      $databases_names{$F[0]} = $F[0];
    }
  }
  close $DB;
} elsif ( $db_string ) {
  @databases_makedb = split( m/\,/, $db_string );
  @databases_blastdb = split( m/\,/, $db_string );
  if ( $names_string ) { ## get names subs if exists, otherwise just use input db names (makes for nicer table)
    @databases_names{@databases_blastdb} = split( m/\,/, $names_string ); ##key= full db name; val= sub name
  } else {
    @databases_names{@databases_blastdb} = split( m/\,/, $db_string ); ##key= full db name; val= full db name
  }
}

## shout about scoring system
if ( $use_qcovhsp_as_score ) {
  print STDERR "[INFO] Scoring system set to 'qcovhsp'\n";
} else {
  print STDERR "[INFO] Scoring system set to '1, 0.8, 0.2, 0'\n";
}

print STDERR "[INFO] Check/make blastdb's for: \n\t" . join("\n\t", @databases_makedb) . "\n";
@databases_makedb = @{ check_blastdbs(\@databases_makedb) }; ## check which database files need makeblastdb run on them
make_blastdbs( \@databases_makedb ); ## and then do it

## open BLAST results file
my $blast_file = $outprefix . "_blastTEags_out.txt";
open (my $BLAST, ">$blast_file") or die $!;
print $BLAST join ("\t", "qacc","sacc","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovhsp","repeat_id","ltr_pos","db","descr","result","score") . "\n";

## glob tag fasta files
my @fasta_files = glob ("$fasta_path/*fasta $fasta_path/*fna $fasta_path/*fa");
print STDERR "[INFO] There are ".scalar(@fasta_files)." files in '$fasta_path'\n";

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
        print $BLAST join ("\t", $line,$repeat_id,$databases_names{$database},$ltr_pos,"partial","PASS",$score) . "\n";
        next LINE;

      } elsif ( $qcovhsp > 50 ) { ## marginal match spanning at least 1 bp over TE/genome boundary
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.5;
        $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
        print $BLAST join ("\t", $line,$repeat_id,$databases_names{$database},$ltr_pos,"partial","MARGINAL",$score) . "\n";
        next LINE;

      } else { ## match that does not span TE/genome boundary by any overlap
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0;
        $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
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
## sorting like this ensures the best BLAST hit is last; easier to extract with code below
open (my $ANNOT_BLAST, "sort -k14,14V -k15,15V -k16,16V -k19,19n $blast_file |") or die $!;
while (my $line = <$ANNOT_BLAST>) {
  chomp $line;
  my @F = split ( m/\t/, $line );
  $top_hits{$F[13]}{$F[14]}{$F[15]}{qcovhsp} = $F[12];
  $top_hits{$F[13]}{$F[14]}{$F[15]}{mismatches} = $F[4];
}
close $ANNOT_BLAST;

# print STDERR Dumper(\%top_hits);

## print condensed results
## to show presence / absence of tags across all databases
my $table_file = $outprefix . "_blastTEags_table.txt";
open (my $TAB, ">$table_file") or die $!;
if ( $mark ) {
  print $TAB "repeat_id";
  my $success = 0;
  foreach ( nsort values %databases_names ) {
    if ( $_ eq $mark ) { ## convoluted way of printing a couple of asterisks, but heyho
      print $TAB "\t**$_";
      $success = 1;
    } else {
      print $TAB "\t$_";
    }
  }
  print $TAB "\n";
  if ($success == 0) {
    print STDERR "[WARN] Mark specified for taxon '$mark', but '$mark' does not exist!\n";
  }
} else {
  print $TAB join ( "\t", "repeat_id", nsort values %databases_names ) . "\n"; ## print header
}

foreach my $repeat_id ( nsort keys %ltr_hash ) {
  ## iterate thru ALL ltrs, not just the ones with blast hits
  print $TAB $repeat_id; ## print ltr_id
  foreach my $database ( nsort values %databases_names ) {
    my $final_score;
    if ( ($top_hits{$repeat_id}{$database}{L}) && ($top_hits{$repeat_id}{$database}{R}) ) { ## if hit exists for both L && R
      $final_score = (($top_hits{$repeat_id}{$database}{L}{qcovhsp} - $top_hits{$repeat_id}{$database}{L}{mismatches}) + ($top_hits{$repeat_id}{$database}{R}{qcovhsp} - $top_hits{$repeat_id}{$database}{R}{mismatches})) / 200;

    } elsif ( ($top_hits{$repeat_id}{$database}{L}) && !($top_hits{$repeat_id}{$database}{R}) ) { ## or only L
      $final_score = ($top_hits{$repeat_id}{$database}{L}{qcovhsp} - $top_hits{$repeat_id}{$database}{L}{mismatches}) / 200;

    } elsif ( !($top_hits{$repeat_id}{$database}{L}) && ($top_hits{$repeat_id}{$database}{R}) ) { ## or only R
      $final_score = ($top_hits{$repeat_id}{$database}{R}{qcovhsp} - $top_hits{$repeat_id}{$database}{R}{mismatches}) / 200;

    } else { ## else print 0
      $final_score = 0;

    }
    print $TAB "\t$final_score";
  }
  print $TAB "\n";
  print $TAB "###\n" if ( $eyeball );
}
close $TAB;

print STDERR "[####] Done!\n";
print STDERR "[####] " . `date`;


#################### SUBS

sub check_progs {
  chomp( my $makeblastdb_path = `which makeblastdb` );
  chomp( my $blastn_path = `which blastn` );
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
  my $full_path;
  for my $i ( 0 .. $#in ) {
    ($full_path = $in[$i]) =~ s/^~(\w*)/ ( getpwnam( $1 || $ENV{USER} ))[7] /e; ## to interpret home '~' correctly
    # print STDERR "File: $in[$i]\n";
    # print STDERR "Full path: $full_path\n";
    ## file not exist
    if (! -f "$full_path") {
      die "[ERROR] File '$full_path' does not exist! $!\n\n";
    }
    ## file is gzipped
    if ($full_path =~ m/gz$/) {
      die "[ERROR] Please gunzip your fasta file: $full_path\n\n";
    }
    ## blastdb already made, remove from @out
    if ( (-f "$full_path.nhr") && (-f "$full_path.nin") && (-f "$full_path.nsq") ) {
      print STDERR "[INFO] BlastDB already exists for '$full_path'\n";
      # splice( @out, $i, 1 );
    } else {
      push ( @out, $full_path );
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
