#!/usr/bin/env perl

use strict;
use warnings;

use Bio::SeqIO;
use Data::Dumper;
use Sort::Naturally;
use Getopt::Long qw(:config no_ignore_case);

my $usage = "
SYNOPSIS
  Creates TE-tags from LTR_retriever GFF3 input.

OPTIONS [*] = required
  -g|--gff     [FILE] : LTR_retriever *.pass.list.gff3 file [*]
  -f|--fasta   [FILE] : reference fasta file [*]
  -o|--out   [STRING] : output prefix for outfiles ['results_createTEags']
  -g|--overhang [INT] : require at least this number of bases as genome tag [50]
  -h|--help           : this message
\n";

## input
my ( $in_gfffile, $in_fastafile, $help, $debug );
## defaults
my $out_prefix = "results";
my $overhang_threshold = 50;


GetOptions (
  'g|gff=s' => \$in_gfffile,
  'f|fasta=s' => \$in_fastafile,
  'o|out:s' => \$out_prefix,
  'h|help' => \$help,
  'd|debug' => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ($in_gfffile && $in_fastafile);

print STDERR "[####] TE-EVOLUTION create_tags.pl\n";
print STDERR "[####] " . `date`;

## make dir structure for sequences
my $base_dir = $out_prefix . "_createTEags";
system rm    => '-rf' => $base_dir if (-d $base_dir); ## delete it if exists, it's a brutal world
system mkdir => '-p'  => "$base_dir/sequence/";
# system mkdir => '-p'  => "$base_dir/sequence/righties";

## parse fasta sequence file
my %scaffolds_hash;
my $in = Bio::SeqIO -> new ( -file => $in_fastafile, -format => "fasta" );
while (my $seq_obj = $in -> next_seq()) {
  $scaffolds_hash{$seq_obj->display_id()} = $seq_obj;
}
print STDERR "[INFO] Retrieved ".commify(scalar(keys %scaffolds_hash))." sequences from $in_fastafile\n";

## parse GFF3 file for LTR positions
my $count = 1;
my (%seen_already, %repeat_regions);
my ($repeat_id, $ltr_id);
open (my $GFF, $in_gfffile) or die $!;
while (<$GFF>) {
  chomp;
  next if (m/^\#/);
  my @F = split (/\s+/, $_);
  ## get repeat ID
  if ( $F[2] eq "repeat_region" ) {
    if ($F[8] =~ m/ID=(\S+)/) {
      $repeat_id = $1;
      $repeat_regions{$repeat_id} = $F[0];
    }
  }
  ## get TE tag sequences
  if ( $F[2] eq "long_terminal_repeat" ) { ## possible it might change in future versions
    open (my $FA, ">>$base_dir/sequence/$repeat_id.fa");
    ## get fasta header construct
    if ($F[8] =~ m/Parent=(\S+)/) {
      $ltr_id = $1;
    }
    if ( $seen_already{$repeat_id} ) {
      ## LTR is a rightie
      my $left_coord = ($F[4]-$overhang_threshold+1) < 0 ? 0 : ($F[4]-$overhang_threshold+1); ## left coord cannot be < 0
      my $right_coord = $scaffolds_hash{$F[0]}->length() < ($F[4]+$overhang_threshold) ? $scaffolds_hash{$F[0]}->length() : ($F[4]+$overhang_threshold); ## right coord cannot be > seq length
      ## print TE tag
      print $FA ">$ltr_id:R:$left_coord..$right_coord\n";
      print $FA $scaffolds_hash{$F[0]} -> subseq( $left_coord,$right_coord ) . "\n";
      ## above should return this:
      ## -----[----LTR2----]---
      ##                 ===== tag2
    } else {
      ## LTR is a leftie
      my $left_coord = ($F[3]-$overhang_threshold) < 0 ? 0 : ($F[3]-$overhang_threshold); ## left coord cannot be < 0
      my $right_coord = $scaffolds_hash{$F[0]}->length() < ($F[3]+$overhang_threshold-1) ? $scaffolds_hash{$F[0]}->length() : ($F[3]+$overhang_threshold-1); ## right coord cannot be > seq length
      ## print TE tag
      print $FA ">$ltr_id:L:$left_coord..$right_coord\n";
      print $FA $scaffolds_hash{$F[0]} -> subseq( $left_coord,$right_coord ) . "\n";
      ## above should return this:
      ## -----[----LTR1----]---
      ##    ===== tag1
    }
    $seen_already{$repeat_id}++;
    # print STDERR "$seen_already{$F[0]}\n";
    close $FA if $seen_already{$repeat_id} == 2;
  }
}
close $GFF;
print STDERR "[INFO] Found ".commify(scalar(keys %repeat_regions))." LTR regions in $in_gfffile\n";
print STDERR "[INFO] Finished\n";


#################### SUBS

sub commify {
    my $text = reverse $_[0];
    $text =~ s/(\d\d\d)(?=\d)(?!\d*\.)/$1,/g;
    return scalar reverse $text;
}

sub percentage {
    my $numerator = $_[0];
    my $denominator = $_[1];
    my $places = "\%.2f"; ## default is two decimal places
    if (exists $_[2]){$places = "\%.".$_[2]."f";};
    my $float = (($numerator / $denominator)*100);
    my $rounded = sprintf("$places",$float);
    return "$rounded\%";
}

sub check_progs {
  chomp( my $samtools_path = `which samtools` );
  if (!( $samtools_path )) {
    die "[ERROR] Cannot find samtools in \$PATH\n";
  } else {
    print STDERR "[INFO] Found samtools at $samtools_path\n"; ##
  }
  chomp( my $bedtools_path = `which bedtools` );
  if (!( $bedtools_path )) {
    die "[ERROR] Cannot find BEDTools in \$PATH\n";
  } else {
    print STDERR "[INFO] Found BEDTools at $bedtools_path\n"; ##
  }
}

__END__
