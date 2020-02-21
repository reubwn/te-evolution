#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
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
  -g|--overhang [INT] : require at least this number of bases as genome tag [50]
  -m|--match    [INT] : require at least this number of bases either side of TE boundary [50]
  -e|--evalue   [STR] : BLASTn evalue [1e-20]
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
my $overhang_threshold = 50;
my $match_threshold = 50;
my $evalue = "1e-20";
my $threads = 4;

GetOptions (
  'i|in=s' => \$IN_path,
  'd|db=s' => \$db_string,
  'n|names:s' => \$names_string,
  'o|out:s' => \$OUT_prefix,
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
die $usage unless ( $IN_path && $db_string );

print STDERR "[####] TE-EVOLUTION blast_tags.pl\n";
print STDERR "[####] " . `date`;
##
my ( %ltr_hash, %blast_hash, %lefties_hash, %righties_hash, %results );
my ( $ambiguous ) = ( 0 );

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

# ## parse SAM input file
# ## to generate te-tags with metadata
# open (my $SAM, "samtools view $IN_file |") or die $!;
# while (my $line = <$SAM>) {
#   chomp $line;
#   my @F = split ("\t", $line);
#   ## parse reads and build hash
#   if ( $F[5] =~ m/^(\d+)S/ ) {
#     ## read overhangs left-side of LTR element
#     if ( $1 >= $overhang_threshold ) {
#       ## read has te-tag >= $overhang_threshold
#       $F[0] =~ s/\s.+// if ($F[0] =~ m/\s+/); ## remove trailing text with whitespaces, common in fastq headers
#       push @{ $ltr_hash{$F[2]}{left_names} }, $F[0];
#       push @{ $ltr_hash{$F[2]}{left_cigars} }, $F[5];
#       push @{ $ltr_hash{$F[2]}{left_seqs} }, $F[9];
#
#     }
#   } elsif ( $F[5] =~ m/(\d+)S$/ ) {
#     ## read overhangs right-side of LTR element
#     if ( $1 >= $overhang_threshold ) {
#       ## read has te-tag >= $overhang_threshold
#       $F[0] =~ s/\s.+// if ($F[0] =~ m/\s+/);
#       push @{ $ltr_hash{$F[2]}{right_names} }, $F[0];
#       push @{ $ltr_hash{$F[2]}{right_cigars} }, $F[5];
#       push @{ $ltr_hash{$F[2]}{right_seqs} }, $F[9];
#
#     }
#
#   } else {
#     print STDERR "[INFO] Skipping read $F[0] due to ambiguous cigar $F[5]\n" if ( $verbose );
#     $ambiguous++;
#   }
# }
# close $SAM;
#
# print Dumper \%ltr_hash if ( $debug );
#
# ## make dir structure for sequences
# my $base_dir = $OUT_prefix . "_sequences";
# unlink $base_dir if (-d $base_dir); ## delete it if exists, it's a brutal world
# system mkdir => '-p' => "$base_dir/lefties";
# system mkdir => '-p' => "$base_dir/righties";

## open BLAST results file
my $blast_file = $OUT_prefix . "_blast.txt";
open (my $BLAST, ">$blast_file") or die $!;
print $BLAST join ("\t", "qacc","sacc","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore","qcovhsp","ltr_id","ltr_pos","db","cigar","boundary","descr","result","score") . "\n";

## glob tag fasta files
my @fasta_files = glob ("$IN_path/*fasta $IN_path/*fnaa $IN_path/*fa");
print STDERR "[INFO] There are ".scalar(@fasta_files)." files in '$IN_path'\n";

## iterate across files and dbs; do BLAST
foreach my $fasta_file (@fasta_files) {

  ## BLAST vs each database in turn
  foreach my $database (@databases_blastdb) {
    print STDERR "[INFO] BLASTing query '$fasta_file' TE-tags versus database '$database'...\n";

    ## open blast filehandle
    open (my $IN, "blastn -task megablast -num_threads $threads -evalue $evalue -query $query_file -db $database -outfmt '6 std qcovhsp' |") or die $!;

    ## iterate thru blast results
    while (my $line = <$IN>) {
      chomp $line;
      my ($qacc, $sacc, $pident, $length, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore, $qcovhsp) = split( m/\s+/, $line );

      my $score;
      if ( $qcovhsp == 100 ) { ## successful alignment across entire tag
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 1;
        ## annotate BLAST result
        print $BLAST join ("\t", $line,$query,"L",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"full","PASS",$score) . "\n";

      } elsif ( ($qstart < ($boundaries{$qacc}-$match_threshold)) && ($qend > ($boundaries{$qacc}+$match_threshold)) ) { ## successful match spanning +/- $match_threshold over TE/genome boundary
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.8;
        ## annotate BLAST result
        print $BLAST join ("\t", $line,$query,"L",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","PASS",$score) . "\n";

      } elsif ( ($qstart <= $boundaries{$qacc}) && ($qend => $boundaries{$qacc}) ) { ## successful match spanning TE/genome boundary by any margin
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.2;
        $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
        ## annotate BLAST result
        print $BLAST join ("\t", $line,$query,"L",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","MARGINAL",$score) . "\n";

      } else { ## match that does not span TE/genome boundary by any overlap
        $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0;
        $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
        ## annotate BLAST result
        print $BLAST join ("\t", $line,$query,"L",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","FAIL",$score) . "\n";
      }
    }
    close $IN;
  }
}

# ## print lefties and righties to file for BLASTing; do BLAST
# foreach my $query (nsort keys %ltr_hash) {
#
#   ## ~~~~~~~~~~~~~~~
#   ## process lefties
#   ## ~~~~~~~~~~~~~~~
#   if ( ($ltr_hash{$query}{left_names}) && ($ltr_hash{$query}{left_seqs}) ) {
#     my @names_arr = @{$ltr_hash{$query}{left_names}};
#     my @seqs_arr = @{$ltr_hash{$query}{left_seqs}};
#     my @cigars_arr = @{$ltr_hash{$query}{left_cigars}};
#     my @overhangs_arr = map { m/^(\d+)S/; $1 } @cigars_arr;
#     my %cigars; @cigars{@names_arr} = @cigars_arr; ##key= seqid; val= cigar
#     my %boundaries; @boundaries{@names_arr} = @overhangs_arr; ##key= seqid; val=overhang value from cigar
#     print Dumper \%boundaries if ( $debug );
#
#     ## write to file
#     my $query_file = "$base_dir/lefties/$query.fa";
#     open (my $F, ">$query_file") or die $!;
#     for my $i ( 0 .. $#names_arr ) {
#       print $F ">$names_arr[$i]\n$seqs_arr[$i]\n";
#     }
#     close $F;
#
#     ## BLAST vs each database in turn
#     foreach my $database (@databases_blastdb) {
#       print STDERR "[INFO] BLASTing query '$query' LEFT TE-tags versus database '$database'...\n";
#
#       ## open blast filehandle
#       open (my $IN, "blastn -task megablast -num_threads $threads -evalue $evalue -query $query_file -db $database -outfmt '6 std qcovhsp' |") or die $!;
#
#       ## iterate thru blast results
#       while (my $line = <$IN>) {
#         chomp $line;
#         my ($qacc, $sacc, $pident, $length, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore, $qcovhsp) = split( m/\s+/, $line );
#
#         my $score;
#         if ( $qcovhsp == 100 ) { ## successful alignment across entire tag
#           $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 1;
#           ## annotate BLAST result
#           print $BLAST join ("\t", $line,$query,"L",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"full","PASS",$score) . "\n";
#
#         } elsif ( ($qstart < ($boundaries{$qacc}-$match_threshold)) && ($qend > ($boundaries{$qacc}+$match_threshold)) ) { ## successful match spanning +/- $match_threshold over TE/genome boundary
#           $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.8;
#           ## annotate BLAST result
#           print $BLAST join ("\t", $line,$query,"L",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","PASS",$score) . "\n";
#
#         } elsif ( ($qstart <= $boundaries{$qacc}) && ($qend => $boundaries{$qacc}) ) { ## successful match spanning TE/genome boundary by any margin
#           $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.2;
#           $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
#           ## annotate BLAST result
#           print $BLAST join ("\t", $line,$query,"L",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","MARGINAL",$score) . "\n";
#
#         } else { ## match that does not span TE/genome boundary by any overlap
#           $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0;
#           $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
#           ## annotate BLAST result
#           print $BLAST join ("\t", $line,$query,"L",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","FAIL",$score) . "\n";
#         }
#       }
#       close $IN;
#     }
#   }
#
#   ## ~~~~~~~~~~~~~~~~
#   ## process righties
#   ## ~~~~~~~~~~~~~~~~
#   if ( ($ltr_hash{$query}{right_names}) && ($ltr_hash{$query}{right_seqs}) ) {
#     my @names_arr = @{$ltr_hash{$query}{right_names}};
#     my @seqs_arr = @{$ltr_hash{$query}{right_seqs}};
#     my @cigars_arr = @{$ltr_hash{$query}{right_cigars}};
#     my @overhangs_arr = map { m/(\d+)S$/; $1 } @cigars_arr;
#     my %cigars; @cigars{@names_arr} = @cigars_arr; ##key= seqid; val= cigar
#     my %boundaries; @boundaries{@names_arr} = @overhangs_arr; ##key= seqid; val=overhang value from cigar
#     print Dumper \%boundaries if ( $debug );
#
#     ## write to file
#     my $query_file = "$base_dir/righties/$query.fa";
#     open (my $F, ">$query_file") or die $!;
#     for my $i ( 0 .. $#names_arr ) {
#       print $F ">$names_arr[$i]\n$seqs_arr[$i]\n";
#     }
#     close $F;
#
#     ## BLAST vs each database in turn
#     foreach my $database (@databases_blastdb) {
#       print STDERR "[INFO] BLASTing query '$query' RIGHT TE-tags versus database '$database'...\n";
#       open (my $IN, "blastn -task megablast -num_threads $threads -evalue $evalue -query $query_file -db $database -outfmt '6 std qcovhsp' |") or die $!;
#       while (my $line = <$IN>) {
#         chomp $line;
#         my ($qacc, $sacc, $pident, $length, $mismatch, $gapopen, $qstart, $qend, $sstart, $send, $evalue, $bitscore, $qcovhsp) = split( m/\s+/, $line );
#         my $score;
#         if ( $qcovhsp == 100 ) { ## successful alignment across entire tag
#           $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 1;
#           print $BLAST join ("\t", $line,$query,"R",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"full","PASS",$score) . "\n";
#         } elsif ( ($qstart < ($boundaries{$qacc}-$match_threshold)) && ($qend > ($boundaries{$qacc}+$match_threshold)) ) { ## successful match spanning +/- $match_threshold over TE/genome boundary
#           $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.8;
#           print $BLAST join ("\t", $line,$query,"R",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","PASS",$score) . "\n";
#         } elsif ( ($qstart <= $boundaries{$qacc}) && ($qend => $boundaries{$qacc}) ) { ## successful match spanning TE/genome boundary by any margin
#           $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0.2;
#           $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
#           print $BLAST join ("\t", $line,$query,"R",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","MARGINAL",$score) . "\n";
#         } else { ## match that does not span TE/genome boundary by any overlap
#           $score = ( $use_qcovhsp_as_score ) ? $qcovhsp : 0;
#           $score = 0 if ( $collapse_marginal_scores ); ## collapse marginal calls to score = 0
#           print $BLAST join ("\t", $line,$query,"R",$databases_names{$database},$cigars{$qacc},$boundaries{$qacc},"partial","FAIL",$score) . "\n";
#         }
#       }
#       close $IN;
#     }
#   }
# }
# close $BLAST;
# print Dumper \%results;

## process annotated blast results
## want to save the 'best' score per query-subject pair only
open (my $ANNOT, $blast_file) or die $!;
while (my $line = <$ANNOT>) {
  chomp $line;
  my @F = split ( m/\t/, $line );
  if ( !($results{$F[0]}{$F[15]}) ) { ## first time
    $results{$F[0]}{$F[15]} = $F[20]; ## key= TE-tag name; val= %{ key= database name; val= score }
  } else { ## subsequent
    ## keep highest score
    $results{$F[0]}{$F[15]} = $F[20] if ($F[20] > $results{$F[0]}{$F[15]});
  }
}
close $ANNOT;

## print condensed results
## to show presence / absence of tags across all databases
my $table_file = $OUT_prefix . "_table.txt";
open (my $TAB, ">$table_file") or die $!;
if ( $mark ) {
  print $TAB "ltr_id\tposition\tqacc\t";
  foreach ( nsort values %databases_names ) {
    if ( $_ eq $mark ) { ## convoluted way of printing a couple of asterisks, but heyho
      print $TAB "**$_\t";
    } else {
      print $TAB "$_\t";
    }
  }
  print $TAB "\n";
} else {
  print $TAB join ( "\t", "ltr_id","position","qacc", nsort values %databases_names ) . "\n"; ## print header
}

foreach my $ltr ( nsort keys %ltr_hash ) {
  ## iterate thru ALL ltrs, not just the ones with blast hits
  ## ~~~~~~~~~~
  ## do lefties
  ## ~~~~~~~~~~
  foreach my $name ( nsort @{ $ltr_hash{$ltr}{left_names} } ) {
    print $TAB "$ltr\tL\t$name"; ## print ltr_id and position factors
    foreach my $db ( nsort values %databases_names ) {
      if ( $results{$name}{$db} ) { ## if exists, print score
        print $TAB "\t$results{$name}{$db}";
      } else { ## else print 0
        print $TAB "\t0";
      }
    }
    print $TAB "\n";
  }
  print $TAB "###\n" if ( $eyeball );

  ## ~~~~~~~~~~~~~
  ## then righties
  ## ~~~~~~~~~~~~~
  foreach my $name ( nsort @{ $ltr_hash{$ltr}{right_names} } ) {
    print $TAB "$ltr\tR\t$name"; ## print ltr_id and position factors
    foreach my $db ( nsort values %databases_names ) {
      if ( $results{$name}{$db} ) { ## if exists, print score
        print $TAB "\t$results{$name}{$db}";
      } else { ## else print 0
        print $TAB "\t0";
      }
    }
    print $TAB "\n";
  }
  print $TAB "###\n" if ( $eyeball );
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
