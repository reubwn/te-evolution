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
  Assembles input fasta sequences using super simple de Bruijn graph assembly approach
OPTIONS [*] = required
  -i|--fasta [FILE]   : fasta file of sequences [*]
  -k|--kmer  [INT]    : Kmer size [51]
  -h|--help           : this message
  -d|--debug          : debug mode
\n";

## input
my (
  $fasta_infile,
  $help, $print_kmers, $debug
);
## defaults
my $kmer_size = 51;

GetOptions (
  'i|fasta=s' => \$fasta_infile,
  'k|kmer:i' => \$kmer_size,
  'h|help' => \$help,
  'p|print' => \$print_kmers,
  'd|debug' => \$debug
);
## help and usage
die $usage if $help;
die $usage unless ($fasta_infile);

print STDERR "[####] TE-EVOLUTION assemble_reads.pl\n";
print STDERR "[####] " . `date`;

my $seqin = Bio::SeqIO->new( -file => $fasta_infile, -format => "fasta" );

my $DEAD;
if ($print_kmers) {
  open ($DEAD, ">trimmed_kmers.txt") or die $!;
}

## read in sequences and get unique kmers
my %kmer_population;
while ( my $seq_obj = $seqin->next_seq() ) {
  ## get forward kmers
  foreach (@{ get_kmers($seq_obj->seq()) }) {
    $kmer_population{$_}++; ##key= kmer; val= count
  }
  ## get revcomp kmers and add to population; pretty sure this is NOT the most efficient way to do it but...
  foreach (@{ get_kmers($seq_obj->revcom()->seq()) }) {
    $kmer_population{$_}++;
  }
}
my @kmers = keys %kmer_population; ## put kmers into array
print Dumper(\%kmer_population) if $debug;

# my @kmers_trimmed = @{ trim_kmers( \@kmers, 5 ) };
# print Dumper(\%graph) if $debug;

## find all paths through de Bruijn graph
my @paths = contig_generation( make_graph(\@kmers, $kmer_size) );
print Dumper \@paths if $debug;

## generate contig sequences from paths
my $index = 1;
my %contigs;
foreach my $path (@paths) {
  my $path_coverage = scalar @{ $path };
  my $sequence = substr $path->[0], 0, -1;
  $sequence .= join q{}, map { substr $_, -1 } @{$path};
  # push @contigs, $sequence;
  $contigs{"path".$index."|len".(length($sequence))."|cov".$path_coverage} = $sequence;
  $index++;
}

## remove identical revcomp $contigs
%contigs = %{ remove_revcomps( \%contigs ) };

# %contigs = %{ remove_similar (\%contigs, 5) };

%contigs = %{ clustalo_pairwise( \%contigs ) };

# my $test = "ATCTGGGGGGG";
# my $rc = revcomp($test);
# print "$rc\n";
## print @contigs
print STDERR "[INFO] Number of seqs in \%contigs: ".scalar(keys %contigs)."\n";
foreach (sort {length($contigs{$b})<=>length($contigs{$a})} keys %contigs) {
  print ">$_\n$contigs{$_}\n";
}
# printf "%s\n", join "\n", sort @contigs;

############## SUBS

sub get_kmers {
  my $seq = $_[0];
  my $len = length $seq;

  my @kmers;
  for (my $i = 0; $i + $kmer_size <= $len; $i++) {
    my $kmer = substr($seq, $i, $kmer_size);
    push ( @kmers, $kmer) unless $kmer =~ m/N/; ##discard kmers with Ns in them
  }

  return \@kmers;
 }

sub hamming {
  return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

sub trim_kmers {
  my @kmers = @{ $_[0] };
  my $threshold = $_[1];
  my %kmers_trimmed = map { $_ => 1 } @kmers;
  for my $i (0 .. $#kmers) {
    for my $j (0 .. $#kmers) {
      my $hamming = hamming ($kmers[$i], $kmers[$j]);
      next if $hamming == 0; ## skip identical kmers
      if ($hamming < $threshold) {
        my $delete = $kmer_population{$kmers[$i]} < $kmer_population{$kmers[$j]} ? $kmers[$i] : $kmers[$j]; ## get the kmer with LOWER count (more likely error)
        print $DEAD join ("\t", $delete, $kmer_population{$delete}, $hamming) . "\n" if $print_kmers;
        delete $kmers_trimmed{$delete};
      }
    }
  }
  print STDERR "[INFO] Trimmed ".(scalar(@kmers) - scalar(keys %kmers_trimmed))." from kmer population\n";
  my @kmers_trimmed = keys %kmers_trimmed;
  return \@kmers_trimmed;
}

 sub make_graph {
   ## make debruijn graph where nodes are k-1-mers and edges between nodes are kmers
   my @kmers = @{ $_[0] };

   my %graph;
   foreach my $kmer (@kmers) {
     my $head = substr ($kmer, 0, -1);
     my $tail = substr ($kmer, 1);
     $graph{$head}{$tail}++;
   }

   return \%graph;
 }

sub contig_generation {
  my ($graph) = @_;

  my %in;
  my %out;
  foreach my $node1 ( keys %{$graph} ) {
    foreach my $node2 ( keys %{ $graph->{$node1} } ) {
      $out{$node1} += $graph->{$node1}{$node2};
      $in{$node2}  += $graph->{$node1}{$node2};
      if ( !exists $in{$node1} ) {
        $in{$node1} = 0;
      }
      if ( !exists $out{$node2} ) {
        $out{$node2} = 0;
      }
    }
  }

  my @paths;

  foreach my $node ( keys %{$graph} ) {
    next if $in{$node} == 1 && $out{$node} == 1;
    foreach ( 1 .. $out{$node} ) {
      while ( exists $graph->{$node} ) {
        my $next_node = ( keys %{ $graph->{$node} } )[0];
        my @path = ( $node, $next_node );
        while ( $in{$next_node} == 1 && $out{$next_node} == 1 ) {
          $next_node = ( keys %{ $graph->{$next_node} } )[0];
          push @path, $next_node;
        }
        push @paths, \@path;
        $graph = remove_path( $graph, \@path );
      }
    }
  }

  return @paths;
}

sub remove_path {
  my ( $graph, $path ) = @_;

  foreach my $i ( 1 .. scalar @{$path} - 1 ) {
    my $node1 = $path->[ $i - 1 ];
    my $node2 = $path->[$i];
    $graph->{$node1}{$node2}--;
    if ( $graph->{$node1}{$node2} == 0 ) {
      delete $graph->{$node1}{$node2};
      if ( scalar keys %{ $graph->{$node1} } == 0 ) {
        delete $graph->{$node1};
      }
    }
  }

  return $graph;
}

sub revcomp {
  (my $revcomp = reverse $_[0]) =~ tr/ATCGatcg/TAGCtagc/;
  return $revcomp;
}

sub remove_revcomps {
  ## remove contigs that are direct revcomps
  my %contigs = %{ $_[0] };
  my %revcomp = map { revcomp $_ => 1} values %contigs;
  foreach (keys %contigs) {
    my $rc = revcomp $contigs{$_};
    if (exists $revcomp{$rc}) {
      delete $revcomp{$_};
      delete $revcomp{$rc};
      delete $contigs{$_};
      next;
    }
  }
  return \%contigs;
}

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
      if ( system ("clustalo -i tmp.fa -t DNA -o tmp.aln --force --threads=4") !=0 ) {
        die "[ERROR] Clustalo didn't work\n";
      }
      my $in = Bio::AlignIO->new( -file => "tmp.aln", -format => "fasta" );
      my $aln_obj = $in->next_aln();
      if ( $aln_obj->percentage_identity() > 98 ) {
        my $delete = length( $contigs{$seq1} ) < length( $contigs{$seq2} ) ? $seq1 : $seq2;
        delete $contigs_new{$delete};
        print join ("\t", $seq1, $seq2, hamming($contigs{$seq1},$contigs{$seq2}), $aln_obj->percentage_identity()) . "\n";
      }
    }
  }
  return \%contigs_new;
}


# sub remove_similar {
#   ## remove contigs that are very similar to other contigs
#   my %contigs = %{ $_[0] };
#   my $threshold = $_[1];
#   my %contigs_new = %contigs;
#   my %contigs_revcomp = map { $_ => revcomp $contigs{$_}} keys %contigs;
#   # print STDERR "[INFO] Number of seqs in \%contigs: ".scalar(keys %contigs)."\n";
#   foreach my $x (keys %contigs) {
#     foreach my $y1 (keys %contigs) {
#       next if $x eq $y1; ## skip comparisons to self
#       my $hamming = hamming( $contigs{$x}, $contigs{$y1} );
#       if ($hamming < $threshold) {
#         my $delete = length $contigs{$x} < length $contigs{$y1} ? $x : $y1; ## seqid of smaller
#         delete $contigs_new{$delete};
#       }
#     }
#     foreach my $y2 (keys %contigs_revcomp) {
#       next if $x eq $y2; ## skip comparisons to self
#       my $hamming = hamming( $contigs{$x}, $contigs{$y2} );
#       if ($hamming < $threshold) {
#         my $delete = length $contigs{$x} < length $contigs{$y2} ? $x : $y2; ## seqid of smaller
#         delete $contigs_new{$delete};
#       }
#     }
#   }


  # foreach my $x (keys %contigs) {
  #   foreach my $y1 (@forward) {
  #     my $hamming = hamming( $contigs{$x}, $y1 );
  #     next if $hamming == 0; ## skip comparisons to self
  #     print "$x\n$contigs{$x}\n$y1\n$hamming\n" if $hamming <= $threshold;
  #   }
  #
  #   foreach my $y2 (@revcomp) {
  #     my $hamming = hamming( $contigs{$x}, $y2 );
  #     next if $hamming == 0;
  #     print "$x\n$contigs{$x}\n$y2 (rev)\n$hamming\n" if $hamming <= $threshold;
  #   }
  # }
#   return \%contigs_new;
# }

__END__

TGCAA
TTGCA
