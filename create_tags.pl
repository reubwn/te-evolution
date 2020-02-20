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

## check for samtools
check_progs();

## make dir structure for sequences
my $base_dir = $out_prefix . "_createTEags";
system rm    => '-rf' => $base_dir if (-d $base_dir); ## delete it if exists, it's a brutal world
system mkdir => '-p'  => "$base_dir/sequence/lefties";
system mkdir => '-p'  => "$base_dir/sequence/righties";

## parse fasta sequence file
my %scaffolds_hash;
my $in = Bio::SeqIO -> new ( -file => $in_fastafile, -format => "fasta" );
while (my $seq_obj = $in -> next_seq()) {
  $scaffolds_hash{$seq_obj->display_id()} = $seq_obj;
}
print STDERR "[INFO] Retrieved ".commify(scalar(keys %scaffolds_hash))." sequences from $in_fastafile\n";

## parse GFF3 file for LTR positions
my $count = 1;
my %seen_already;
open (my $GFF, $in_gfffile) or die $!;
while (<$GFF>) {
  chomp;
  next if (m/^\#/);
  my @F = split (/\s+/, $_);
  if ( $F[2] eq "long_terminal_repeat" ) { ## possible it might change in future versions
    if ( $seen_already{$F[0]} ) {
      ## LTR is a rightie
      my $filename = $F[0]."_R_".$seen_already{$F[0]}++.".fasta";
      open (my $FA, ">$base_dir/sequence/righties/$filename")
      print $FA ">$filename\n";
      print $FA $scaffolds_hash{$F[0]} -> subseq(($F[4]-$overhang_threshold),($F[4]+$overhang_threshold)) . "\n";
      close $FA;
      ## above should return this:
      ## -----[----LTR2----]---
      ##                 ===== tag2
    } else {
      ## LTR is a leftie
      my $filename = $F[0]."_L_".$seen_already{$F[0]}++.".fasta";
      open (my $FA, ">$base_dir/sequence/lefties/$filename")
      print $FA ">$filename\n";
      print $FA $scaffolds_hash{$F[0]} -> subseq(($F[3]-$overhang_threshold),($F[3]+$overhang_threshold)) . "\n";
      close $FA;
      ## above should return this:
      ## -----[----LTR1----]---
      ##    ===== tag1
    }
    $seen_already{$F[0]}++;
    $count++;
  }
}
close $GFF;

# my $usage = "
# SYNOPSIS
#   Creates TE-tags from SAM/BAM input.
#
# OPTIONS [*] = required
#   -i|--in      [FILE] : Reads mapped to LTRs (BAM/SAM) [*]
#   -o|--out   [STRING] : output prefix for outfiles ['results_createTEags']
#   -g|--overhang [INT] : require at least this number of bases as genome tag [50]
#   -h|--help           : this message
# \n";
#
# ## input
# my ( $in_samfile, $help, $debug );
# ## defaults
# my $out_prefix = "results";
# my $overhang_threshold = 50;
#
#
# GetOptions (
#   'i|in=s' => \$in_samfile,
#   'o|out:s' => \$out_prefix,
#   'h|help' => \$help,
#   'd|debug' => \$debug
# );
# ## help and usage
# die $usage if $help;
# die $usage unless ($in_samfile);
#
# print STDERR "[####] TE-EVOLUTION create_tags.pl\n";
# print STDERR "[####] " . `date`;
#
# ## check for samtools
# check_progs();
#
# ## make dir structure for sequences
# my $base_dir = $out_prefix . "_createTEags";
# system rm    => '-rf' => $base_dir if (-d $base_dir); ## delete it if exists, it's a brutal world
# system mkdir => '-p'  => "$base_dir/sequence/lefties";
# system mkdir => '-p'  => "$base_dir/sequence/righties";
#
# my ( %ltr_coords, %ltr_hash, %results );
# my ( $lefties, $righties, $ambiguous, $l_i, $r_i ) = (0,0,0,0,0);
#
# ## parse SAM input file -H
# ## to get ref lengths
# open (my $INSAM1, "samtools view -H $in_samfile |") or die $!;
# while (my $line = <$INSAM1>) {
#   chomp $line;
#   if ( $line =~ m/^\@SQ/ ) { ## probably unneccesary
#     if ( $line =~ m/SN:(\w+)\tLN:(\d+)/ ) {
#       $ltr_coords{$1} = $2; ## key= LTR name; val= LTR length
#     }
#   }
# }
# close $INSAM1;
#
# ## open BED file of overhangs
# open (my $BED, ">$base_dir/overhangs.bed") or die $!;
# foreach (nsort keys %ltr_coords) {
#   ## bedtools - select ONLY reads that overhang into TE-part by 50 bases or more
#   print $BED "$_\t0\t$overhang_threshold\tleft_overhang\n$_\t" . ($ltr_coords{$_} - $overhang_threshold) . "\t$ltr_coords{$_}\tright_overhang\n";
# }
# close $BED;
#
# ## open BAM file for output
# open (my $OUTSAM, ">$base_dir/accepted_reads.sam") or die $!;
# open (my $OUTFASTA, ">$base_dir/accepted_reads.fasta") or die $!;
# open (my $OUTMAP, ">$base_dir/names_map.txt") or die $!;
#
# ## parse SAM input file
# ## to generate te-tags with metadata - select ONLY reads that overhang into genome-part by 50 bases or more
# open (my $INSAM2, "bedtools intersect -a $in_samfile -b $base_dir/overhangs.bed | samtools view -h |") or die $!;
# while (my $line = <$INSAM2>) {
#   ## pipe header lines straight to outbam
#   if ( $line =~ m/^\@/ ) {
#     print $OUTSAM $line;
#     next;
#   }
#
#   my @F = split ("\t", $line);
#   ## parse reads and build hash
#   if ( $F[5] =~ m/^(\d+)S/ ) {
#     ## read overhangs left-side of LTR element
#     if ( $1 >= $overhang_threshold ) {
#       ## read has te-tag >= $overhang_threshold
#       $F[0] =~ s/\s.+// if ($F[0] =~ m/\s+/); ## remove trailing text with whitespaces, common in fastq headers
#       push @{ $ltr_hash{$F[2]}{left_names} }, $F[0];
#       push @{ $ltr_hash{$F[2]}{left_simplified} }, "left_$l_i"; ## simplified name, make things easier
#       push @{ $ltr_hash{$F[2]}{left_cigars} }, $F[5];
#       push @{ $ltr_hash{$F[2]}{left_seqs} }, $F[9];
#       print $OUTSAM $line; ## read is accepted
#       print $OUTFASTA ">$F[0]\n$F[9]\n";
#       print $OUTMAP "left_$l_i\t$F[0]\n";
#       $lefties++;
#       $l_i++;
#
#     }
#   } elsif ( $F[5] =~ m/(\d+)S$/ ) {
#     ## read overhangs right-side of LTR element
#     if ( $1 >= $overhang_threshold ) {
#       ## read has te-tag >= $overhang_threshold
#       $F[0] =~ s/\s.+// if ($F[0] =~ m/\s+/);
#       push @{ $ltr_hash{$F[2]}{right_names} }, $F[0];
#       push @{ $ltr_hash{$F[2]}{right_simplified} }, "right_$r_i";
#       push @{ $ltr_hash{$F[2]}{right_cigars} }, $F[5];
#       push @{ $ltr_hash{$F[2]}{right_seqs} }, $F[9];
#       print $OUTSAM $line; ## read is accepted
#       print $OUTFASTA ">$F[0]\n$F[9]\n";
#       print $OUTMAP "right_$r_i\t$F[0]\n";
#       $righties++;
#       $r_i++;
#
#     }
#
#   } else {
#     $ambiguous++;
#   }
# }
# close $INSAM2;
# close $OUTSAM;
# close $OUTFASTA;
# print Dumper \%ltr_hash if ($debug);
#
# print STDERR "[INFO] Number of LTR sequences: " . scalar(keys %ltr_hash) . "\n";
# print STDERR "[INFO] Total number of left-overhanging reads: $lefties\n";
# print STDERR "[INFO] Total number of right-overhanging reads: $righties\n";
# print STDERR "[INFO] Total number of reads skipped due to ambiguous CIGAR: $ambiguous\n" if ( $ambiguous > 0 );
#
# open (my $OUTTAB, ">$base_dir/te_tags.txt") or die $!;
#
# ## process LTR reads
# foreach my $query (nsort keys %ltr_hash) {
#
#   ## ~~~~~~~~~~~~~~~
#   ## process lefties
#   ## ~~~~~~~~~~~~~~~
#   if ( ($ltr_hash{$query}{left_names}) && ($ltr_hash{$query}{left_seqs}) ) {
#     my @names_arr = @{$ltr_hash{$query}{left_names}};
#     my @simplified_arr = @{$ltr_hash{$query}{left_simplified}};
#     my @seqs_arr = @{$ltr_hash{$query}{left_seqs}};
#     my @cigars_arr = @{$ltr_hash{$query}{left_cigars}};
#     my @overhangs_arr = map { m/^(\d+)S/; $1 } @cigars_arr;
#     my %cigars; @cigars{@names_arr} = @cigars_arr; ##key= seqid; val= cigar
#     my %boundaries; @boundaries{@names_arr} = @overhangs_arr; ##key= seqid; val=overhang value from cigar
#     print Dumper \%boundaries if ( $debug );
#
#
# }
#
#
#
#
#
#
#
#
#
#
# #####
# foreach (nsort keys %ltr_hash) {
#   if ( $ltr_hash{$_}{left_names} ) {
#     my @left_names = @{ $ltr_hash{$_}{left_names} };
#     my @left_names = @{ $ltr_hash{$_}{left_names} };
#     open (my $L, ">$base_dir/sequence/lefties/$_.fasta") or die $!;
#     for my $i ( 0 .. $#left_names ) {
#       print $L ">$left_names[$i]\n${$ltr_hash{$_}{left_seqs}}[$i]\n"; ## !
#     }
#     close $L;
#   }
#
#   if ( $ltr_hash{$_}{right_names} ) {
#     my @right_names = @{ $ltr_hash{$_}{right_names} };
#     open (my $R, ">$base_dir/sequence/righties/$_.fasta") or die $!;
#     for my $i ( 0 .. $#right_names ) {
#       print $R ">$right_names[$i]\n${$ltr_hash{$_}{right_seqs}}[$i]\n"; ## !
#     }
#     close $R;
#   }
# }

####################
#################### SUBS
####################

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


#
# my %seq_hash;
#
# my $seqio = Bio::SeqIO->new( -file => $fasta_infile, -format => "fasta" );
# while ( my $seq_obj = $seqio->next_seq() ) {
#   $seq_hash{$seq_obj->display_id()} = $seq_obj->seq();
# }
# print "[INFO] Number of raw sequences: ".scalar(keys %seq_hash)."\n";
# print "[INFO] Identity threshold set to: $identity_threshold\%\n";
# print "[INFO] Running clustalo with $threads CPUs... ";
#
# my %filtered = %{ clustalo_pairwise( \%seq_hash ) };
#
# print "done\n";
# print "[INFO] Number of filtered sequences: ".scalar(keys %filtered)."\n";
#
# foreach (sort {length($filtered{$b})<=>length($filtered{$a})} keys %filtered) {
#   print ">$_\n$filtered{$_}\n";
# }
#
# ############## SUBS
#
# sub clustalo_pairwise {
#   if ( system( "clustalo --help &>/dev/null" ) != 0 ) {
#     die "[ERROR] Couldn't find clustalo in \$PATH!\n\n";
#   }
#
#   my %contigs = %{ $_[0] };
#   my %contigs_new = %{ $_[0] };
#   foreach my $seq1 (keys %contigs) {
#     foreach my $seq2 (keys %contigs) {
#       next if $seq1 eq $seq2;
#       open (my $TMP, ">tmp.fa") or die $!;
#       print $TMP ">$seq1\n$contigs{$seq1}\n>$seq2\n$contigs{$seq2}\n";
#       close $TMP;
#       if ( system ("clustalo -i tmp.fa -t DNA -o tmp.aln --force --threads=$threads") !=0 ) {
#         die "[ERROR] Clustalo didn't work, is it operational and in \$PATH?\n";
#       }
#       my $in = Bio::AlignIO->new( -file => "tmp.aln", -format => "fasta" );
#       my $aln_obj = $in->next_aln();
#       if ( $aln_obj->percentage_identity() > $identity_threshold ) {
#         my $delete = length( $contigs{$seq1} ) < length( $contigs{$seq2} ) ? $seq1 : $seq2;
#         delete $contigs_new{$delete};
#         if ($verbose) {
#           if ($delete eq $seq1) { ## *** denotes sequence that is retained
#             print STDERR join ("\t", $seq1, "$seq2***", $aln_obj->percentage_identity()) . "\n";
#           } else {
#             print STDERR join ("\t", "$seq1***", $seq2, $aln_obj->percentage_identity()) . "\n";
#           }
#         }
#       }
#     }
#   }
#   return \%contigs_new;
# }
#
# if ( system("rm tmp.fa tmp.aln") !=0 ) {
#   die "[WARN] Could do clean up\n";
# }
