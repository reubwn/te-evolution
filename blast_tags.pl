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

my @databases = split( m/\,/, $db_string );
print STDERR "[INFO] Check/make blastdb's for: @databases\n";

@databases = @{ check_blastdbs(\@databases) };
make_blastdbs( \@databases );

open (my $SAM, "samtools view $sam_infile |") or die $!;
while (my $line = <$SAM>) {
  chomp $line;
  my @F = split ("\t", $line);
  ## parse reads and build hash
  if ( $F[5] =~ m/^(\d+)S/ ) {
    ## read overhangs left-side of LTR element
    if ( $1 >= $overhang_threshold ) {
      ## read has 'genome' tag >= $overhang_threshold
      $F[0] =~ s/\s.+// if ($F[0] =~ m/\s+/); ## remove trailing text with whitespaces, common in fastq headers
      push @{ $ltr_hash{$F[2]}{left_names} }, $F[0];
      push @{ $ltr_hash{$F[2]}{left_cigars} }, $F[5];
      push @{ $ltr_hash{$F[2]}{left_seqs} }, $F[9];

    }
  } elsif ( $F[5] =~ m/(\d+)S$/ ) {
    ## read overhangs right-side of LTR element
    if ( $1 >= $overhang_threshold ) {
      ## read has 'genome' tag >= $overhang_threshold
      $F[0] =~ s/\s.+// if ($F[0] =~ m/\s+/);
      push @{ $ltr_hash{$F[2]}{right_names} }, $F[0];
      push @{ $ltr_hash{$F[2]}{right_cigars} }, $F[5];
      push @{ $ltr_hash{$F[2]}{right_seqs} }, $F[9];

    }

  } else {
    print STDERR "[INFO] Skipping read $F[0] due to ambiguous cigar $F[5]\n" if ( $verbose );
    $ambiguous++;
  }
}
close $SAM;

print Dumper \%ltr_hash;

## print some information
if ( $verbose ) {
  foreach (nsort keys %ltr_hash) {
    my $lefties = ($ltr_hash{$_}{left_names}) ? scalar(@{$ltr_hash{$_}{left_names}}) : "0";
    my $righties = ($ltr_hash{$_}{right_names}) ? scalar(@{$ltr_hash{$_}{right_names}}) : "0";
    print STDERR "[INFO] $_ has $lefties left reads and $righties right reads\n";
  }
}

## print lefties and righties to file for BLASTing
foreach my $query (nsort keys %ltr_hash) {
  # my ( @left_names, @left_seqs, @right_names, @right_seqs );
  ## process lefties
  if ( ($ltr_hash{$query}{left_names}) && ($ltr_hash{$query}{left_seqs}) ) {
    my @left_names = @{$ltr_hash{$query}{left_names}};
    my @left_seqs = @{$ltr_hash{$query}{left_seqs}};
    my @left_cigars = @{$ltr_hash{$query}{left_cigars}};

    ## write to file
    open (my $F, ">$query.lefties.fa") or die $!;
    for my $i ( 0 .. $#left_names ) {
      print $F ">$left_names[$i]\n$left_seqs[$i]\n";
    }
    close $F;

    ## BLAST vs each database in turn
    foreach my $database (@databases) {
      open (my $BLAST, "blastn -num_threads $threads -task megablast -evalue 1e-20 -query $query -db -outfmt '6 std qcovus' |") or die $!;
      while (my $line = <$BLAST>) {

      }
    }

  }

  ## process righties
  # if ( ($ltr_hash{$query}{right_names}) && ($ltr_hash{$query}{right_seqs}) ) {
  #   my @right_names = @{$ltr_hash{$query}{right_names}};
  #   my @right_seqs = @{$ltr_hash{$query}{right_seqs}};
  #
  #   open (my $F, ">$query.righties.fa") or die $!;
  #   for my $i ( 0 .. $#right_names ) {
  #     print $F ">$right_names[$i]\n$right_seqs[$i]\n";
  #   }
  #   close $F;
  # }


}


# ################### SUBS

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
  my @databases = @{ $_[0] };
  for my $i ( 0 .. $#databases ) {
    if ( (-f "$databases[$i].nhr") && (-f "$databases[$i].nin") && (-f "$databases[$i].nsq") ) {
      splice( @databases, $i, 1);
    }
  }
  return \@databases;
}

sub make_blastdbs {
  my @databases = @{ $_[0] };

  return if ( scalar(@databases) == 0 );

  if ( &parallel == 0 ) {
    if ( system("parallel -j $threads 'makeblastdb -in {} -dbtype nucl' ::: @databases") != 0 ) {
      die "[ERROR] Something wrong with makeblastdb command\n";
    }
  } else {
    foreach (@databases) {
      if ( system("makeblastdb -in $_ -dbtype nucl") != 0 ) {
        die "[ERROR] Something wrong with makeblastdb command: $_\n";
      }
    }
  }
}

sub run_blast {
  my $seq = $_[0];
  open (my $BLAST, "blastn -query $seq")

}
