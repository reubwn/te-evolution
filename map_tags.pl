#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Sort::Naturally;
use List::Util qw/sum/;
use Getopt::Long qw(:config no_ignore_case);

## TODO


my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -i|--in      [FILE] : list of SAM/BAM files to search (txt) [*]
  -f|--fasta   [PATH] : path to 'sequence/' dir output from 'create_tags.pl' [*]
  -d|--dbs   [STRING] : comma delim list of databases to search (BAM/SAM) [alternative to -i]
  -n|--names [STRING] : comma delim list of column names to print, instead of full db names; must be same order as '--dbs'
  -o|--out   [STRING] : output prefix for outfiles ['results_mapTEags']
  -t|--threads  [INT] : number of threads for multicore operations [4]
  -a|--mark     [STR] : mark focal individual with special symbol in output [asterisks '**']
  -h|--help           : this message

DETAILS
  --in is a simple list of paths to SAM/BAM files, with sample ID as second column
  eg. '~/path/to/mapping_foo.bam foo'

  Example of mapping command with BBMap:
  > bbmap.sh \
      ref=\$REF in=\$READS1 in2=\$READS2 \
      outm=mapped_\${PREFIX}.sam.gz \
      minid=0.5 local=t nodisk
\n";

## input
my (
  $infile, $fasta_path, $db_string, $names_string, $use_qcovhsp_as_score, $collapse_marginal_scores, $mark, $eyeball,
  $help, $verbose, $debug
  );
## defaults
my $outprefix = "results";
my $overhang_threshold = 50;
my $match_threshold = 50;
my $evalue = "1e-20";
my $threads = 4;

GetOptions (
  'i|in:s' => \$infile,
  'f|fasta=s' => \$fasta_path,
  'd|db:s' => \$db_string,
  'n|names:s' => \$names_string,
  'o|out:s' => \$outprefix,
  'g|overhang:i' => \$overhang_threshold,
  'm|match:i' => \$match_threshold,
  't|threads:i' => \$threads,
  'c|collapse' => \$collapse_marginal_scores,
  'a|mark:s' => \$mark,
  'y|eyeball' => \$eyeball,
  'h|help' => \$help,
  'verbose' => \$verbose,
  'debug' => \$debug
  );
## help and usage
die $usage if $help;
die $usage unless ( $infile || $db_string );

print STDERR "[####] TE-EVOLUTION map_tags.pl\n";
print STDERR "[####] " . `date`;

## stuff
my ( %repeat_hash, %ltr_hash, %sam_hash );

## check system for required programs
check_progs();

## get dbs to blast against
my (@databases_sams, %databases_names);
if ($infile) {
  print STDERR "[INFO] Using SAM files found in '$infile' as input\n";
  open (my $CONFIG, $infile) or die $!;
  while (<$CONFIG>) {
    chomp;
    my @F = split (m/\s+/);
    if (scalar(@F)==2) {
      push @databases_sams, $F[0];
      $databases_names{$F[0]} = $F[1];
    } else {
      push @databases_sams, $F[0];
      $databases_names{$F[0]} = $F[0];
    }
  }
  close $CONFIG;
} elsif ( $db_string ) {
  @databases_sams = split( m/\,/, $db_string );
  if ( $names_string ) { ## get names subs if exists, otherwise just use input db names (makes for nicer table)
    @databases_names{@databases_sams} = split( m/\,/, $names_string ); ##key= full db name; val= sub name
  } else {
    @databases_names{@databases_sams} = split( m/\,/, $db_string ); ##key= full db name; val= full db name
  }
}

## parse seq info for LTR ids
my @fasta_files = glob ("$fasta_path/*fasta $fasta_path/*fna $fasta_path/*fa");
print STDERR "[INFO] There are ".scalar(@fasta_files)." files in '$fasta_path'\n";
foreach my $fasta_file (nsort @fasta_files) {
  my $repeat_id = `basename -s .fa $fasta_file`; chomp $repeat_id;
  open (my $FA, "grep '>' $fasta_file |") or die $!;
  while (<$FA>) {
    chomp;
    (my $ltr_id = $_) =~ s/\>//;
    $repeat_hash{$repeat_id}{left} = $ltr_id if ($ltr_id =~ /:L:/);
    $repeat_hash{$repeat_id}{right} = $ltr_id if ($ltr_id =~ /:R:/);
    $ltr_hash{$ltr_id} = $repeat_id;
  }
  close $FA;
}

# print Dumper(\%ltr_hash);

## iterate across SAM files
foreach my $database ( @databases_sams ) {
  (my $full_path = $database) =~ s/^~(\w*)/ ( getpwnam( $1 || $ENV{USER} ))[7] /e; ## to interpret home '~' correctly
  # my $subject = `basename $full_path`; ## get the
  open (my $SAM, "samtools view $full_path |") or die $!;
  while (my $line = <$SAM>) {
    my @F = split (m/\t/, $line); ## split on tab not whitespace as some readnames have whitespace
    if (scalar(@F) >= 11) { ## SAM spec at least 11 columns
      ## parse CIGAR strings of mapped reads to comupte score
      my @m = ($F[5] =~ m/(\d+)=/g); ## pull out the number of matches '='
      my @x = ($F[5] =~ m/(\d+)X/g); ## pull out the number of mismatches 'X'
      my $matches = ( sum(@m) ) ? sum(@m) : 0;
      my $mismatches = ( sum(@x) ) ? sum(@x) : 0;
      # print STDERR join("\t", sum(@m), $mismatches) . "\n";
      $sam_hash{$database}{$ltr_hash{$F[2]}}{$F[2]}{(($matches+$mismatches)-$mismatches)}++; ## key= name of samfile; val= %{key= TEag; val=%{key= matches; val= count}}
    } else {
      ## no reads have mapped
      $sam_hash{$database} = ();
      print STDERR "NO READS IN $database\n";
    }
  }
  close $SAM;
}

print Dumper(\%sam_hash) if $debug;

## open MAP scores and counts results file
my $scores_file = $outprefix . "_mapTEags_scores.txt";
my $counts_file = $outprefix . "_mapTEags_counts.txt";
open (my $SCORES, ">$scores_file") or die $!;
open (my $COUNTS, ">$counts_file") or die $!;
if ( $mark ) {
  print $SCORES "repeat_id";
  print $COUNTS "repeat_id";
  my $success = 0;
  foreach ( nsort values %databases_names ) {
    if ( $_ eq $mark ) { ## convoluted way of printing a couple of asterisks, but heyho
      print $SCORES "\t**$_";
      print $COUNTS "\t**$_";
      $success = 1;
    } else {
      print $SCORES "\t$_";
      print $COUNTS "\t$_";
    }
  }
  print $SCORES "\n";
  print $COUNTS "\n";
  if ($success == 0) {
    print STDERR "[WARN] Mark specified for taxon '$mark', but '$mark' does not exist!\n";
  }
} else {
  print $SCORES join ( "\t", "repeat_id", nsort values %databases_names ) . "\n"; ## print header
  print $COUNTS join ( "\t", "repeat_id", nsort values %databases_names ) . "\n"; ## print header
}

## sort and print results
foreach my $repeat_id ( nsort keys %repeat_hash ) {
  print $SCORES $repeat_id; ## print ltr_id
  print $COUNTS $repeat_id; ## print ltr_id
  foreach my $database ( nsort keys %sam_hash ) {
    my @scores = qw/0/;
    my @counts = qw/0/;
    foreach my $ltr_id ( nsort keys %{$sam_hash{$database}{$repeat_id}} ) {
      foreach my $score ( (sort {$b<=>$a} keys %{$sam_hash{$database}{$repeat_id}{$ltr_id}})[0] ) { ## top hit!
        print STDOUT join ("\t", "TOPHIT:",$databases_names{$database},$repeat_id,$ltr_id,$score,$sam_hash{$database}{$repeat_id}{$ltr_id}{$score}) . "\n" if ($debug);
        push (@scores, $score);
        push (@counts, $sam_hash{$database}{$repeat_id}{$ltr_id}{$score});
      }
    }
    # my $final_score = ( sum(@scores) ) ? sum(@scores)/200 : 0;
    # print $SCORES "\t$final_score";
    # my $final_count = ( sum(@counts) ) ? sum(@counts) : 0;
    # print $COUNTS "\t$final_count";
    print $SCORES "\t" . sum(@scores)/200;
    print $COUNTS "\t" . sum(@counts);
  }
  print $SCORES "\n";
  print $COUNTS "\n";
}

print STDERR "[####] Done!\n";
print STDERR "[####] " . `date`;

#################### SUBS

sub check_progs {
  chomp( my $samtools_path = `which samtools` );
  if (!( $samtools_path )) {
    die "[ERROR] Cannot find samtools in \$PATH\n";
  } else {
    print STDERR "[INFO] Found samtools at $samtools_path\n";
  }
}

__END__
