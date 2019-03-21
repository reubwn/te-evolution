#!/usr/bin/env perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "
SYNOPSIS
  Removes near-identical sequences from input.

OPTIONS [*] = required
  -i|--fasta   [FILE] : fasta file of sequences [*]
  -p|--percid  [INT]  : identity threshold for removing near-duplicates [98%]
  -t|--threads [INT]  : number of clustalo threads [4]
  -h|--help           : this message
  -v|--verbose        : verbose mode
\n";

## input
my (
  $fasta_infile,
  $help, $verbose
);
## defaults
my $identity_threshold = 98;
my $threads = 4;

GetOptions (
  'i|fasta=s' => \$fasta_infile,
  'p|percid:i' => \$identity_threshold,
  't|threads:i' => \$threads,
  'h|help' => \$help,
  'v|verbose' => \$verbose
);
## help and usage
die $usage if $help;
die $usage unless ($fasta_infile);

print STDERR "[####] TE-EVOLUTION create_tags.pl\n";
print STDERR "[####] " . `date`;

my %seq_hash;

my $seqio = Bio::SeqIO->new( -file => $fasta_infile, -format => "fasta" );
while ( my $seq_obj = $seqio->next_seq() ) {
  $seq_hash{$seq_obj->display_id()} = $seq_obj->seq();
}
print "[INFO] Number of raw sequences: ".scalar(keys %seq_hash)."\n";
print "[INFO] Identity threshold set to: $identity_threshold\%\n";
print "[INFO] Running clustalo with $threads CPUs... ";

my %filtered = %{ clustalo_pairwise( \%seq_hash ) };

print "done\n";
print "[INFO] Number of filtered sequences: ".scalar(keys %filtered)."\n";

foreach (sort {length($filtered{$b})<=>length($filtered{$a})} keys %filtered) {
  print ">$_\n$filtered{$_}\n";
}

############## SUBS

sub clustalo_pairwise {
  if ( system( "clustalo --help &>/dev/null" ) != 0 ) {
    die "[ERROR] Couldn't find clustalo in \$PATH!\n\n";
  }

  my %contigs = %{ $_[0] };
  my %contigs_new = %{ $_[0] };
  foreach my $seq1 (keys %contigs) {
    foreach my $seq2 (keys %contigs) {
      next if $seq1 eq $seq2;
      open (my $TMP, ">tmp.fa") or die $!;
      print $TMP ">$seq1\n$contigs{$seq1}\n>$seq2\n$contigs{$seq2}\n";
      close $TMP;
      if ( system ("clustalo -i tmp.fa -t DNA -o tmp.aln --force --threads=$threads") !=0 ) {
        die "[ERROR] Clustalo didn't work, is it operational and in \$PATH?\n";
      }
      my $in = Bio::AlignIO->new( -file => "tmp.aln", -format => "fasta" );
      my $aln_obj = $in->next_aln();
      if ( $aln_obj->percentage_identity() > $identity_threshold ) {
        my $delete = length( $contigs{$seq1} ) < length( $contigs{$seq2} ) ? $seq1 : $seq2;
        delete $contigs_new{$delete};
        if ($verbose) {
          if ($delete eq $seq1) { ## *** denotes sequence that is retained
            print STDERR join ("\t", $seq1, "$seq2***", $aln_obj->percentage_identity()) . "\n";
          } else {
            print STDERR join ("\t", "$seq1***", $seq2, $aln_obj->percentage_identity()) . "\n";
          }
        }
      }
    }
  }
  return \%contigs_new;
}

if ( system("rm tmp.fa tmp.aln") !=0 ) {
  die "[WARN] Could do clean up\n";
}
