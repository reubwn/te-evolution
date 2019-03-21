#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use Sort::Naturally;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "
SYNOPSIS

OPTIONS [*] = required
  -i|--sam     [FILE] : BAM/SAM file of reads mapped to LTRs [*]
  -d|--dbs   [STRING] : comma delim list of genomes to search [*]
  -e|--overhang [INT] : require at least this number of bases as genome tag [50]
  -t|--threads  [INT] : number of threads for multicore operations [4]
  -h|--help           : this message
  -v|--verbose        : verbose mode
  -d|--debug          : debug mode
\n";

## input
my (
  $sam_infile, $db_string,
  $help, $verbose
  );
## defaults
my $overhang_threshold = 50;
my $threads = 4;

GetOptions (
  'i|sam=s' => \$sam_infile,
  'd|db=s' => \$db_string,
  'e|overhang:i' => \$overhang_threshold,
  't|threads:i' => \$threads,
  'h|help' => \$help,
  'v|verbose' => \$verbose
  );
## help and usage
die $usage if $help;
die $usage unless ( $sam_infile && $db_string );

print STDERR "[####] TE-EVOLUTION blast_tags.pl\n";
print STDERR "[####] " . `date`;

##
my ( %ltr_hash );
my ( $ambiguous ) = ( 0 );

## check system for required programs
check_progs();

make_blastdbs ( $db_string );

open (my $SAM, "samtools view $sam_infile |") or die $!;
while (my $line = <$SAM>) {
  chomp $line;
  my @F = split ("\t", $line);
  ## parse reads and build hash
  if ( $F[5] =~ m/^(\d+)S/ ) {
    ## read overhangs left-side of LTR element
    if ( $1 >= $overhang_threshold ) {
      ## read has 'genome' tag >= $overhang_threshold
      push @{ $ltr_hash{$F[2]}{left_readnames} }, $F[0];
      push @{ $ltr_hash{$F[2]}{left_cigars} }, $F[5];
      push @{ $ltr_hash{$F[2]}{left_readseqs} }, $F[9];

    }
  } elsif ( $F[5] =~ m/(\d+)S$/ ) {
    ## read overhangs right-side of LTR element
    if ( $1 >= $overhang_threshold ) {
      ## read has 'genome' tag >= $overhang_threshold
      push @{ $ltr_hash{$F[2]}{right_readnames} }, $F[0];
      push @{ $ltr_hash{$F[2]}{right_cigars} }, $F[5];
      push @{ $ltr_hash{$F[2]}{right_readseqs} }, $F[9];

    }

  } else {
    print STDERR "[INFO] Skipping read $F[0] due to ambiguous cigar $F[5]\n" if ( $verbose );
    $ambiguous++;
  }
}
close $SAM;

## gather some information
if ( $verbose ) {
  foreach (nsort keys %ltr_hash) {
    my $lefties = @{$ltr_hash{$_}{left_readnames}} ? scalar(@{$ltr_hash{$_}{left_readnames}}) : "0";
    my $righties = @{$ltr_hash{$_}{right_readnames}} ? scalar(@{$ltr_hash{$_}{right_readnames}}) : "0";
    print STDERR "[INFO] $_ has $lefties left reads and $righties right reads\n";
  }
}

################### SUBS

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
  if (!( $parallel_path )) {
    print STDERR "[INFO] Didn't find parallel in \$PATH, OK\n";
  } else {
    print STDERR "[INFO] Found parallel at $parallel_path\n"; ##
    $result = 0;
  }
  return $result;
}

sub make_blastdbs {
  my $string = $_[0];
  my @split = split ",", $string;
  print STDERR "[INFO] Creating blastdb's from: @split\n";

  if ( &parallel == 0 ) {
    if ( system("parallel --dry-run -j $threads 'makeblastdb -in {} -dbtype nucl' ::: @split") != 0 ) {
      die "[ERROR] Something wrong with makeblastdb command\n";
    }
  } else {
    foreach (@split) {
      if ( system("makeblastdb -in $_ -dbtype nucl") != 0 ) {
        die "[ERROR] Something wrong with makeblastdb command: $_\n";
      }
    }
  }
}

sub run_blast {


}
