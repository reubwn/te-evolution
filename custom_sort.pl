#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;

#my @sort_order = qw ( B C A D );
my @sort_order = qw (
 Unknown Other DNA DNA/Academ DNA/CMC DNA/Crypton DNA/Dada DNA/Ginger DNA/IS3EU
 DNA/IS5 DNA/Kolobok DNA/MULE DNA/Maverick DNA/Merlin DNA/P DNA/PIF DNA/PiggyBac
 DNA/Sola DNA/TcMar DNA/Zator DNA/hAT RC RC/Helitron LTR LTR/Copia LTR/DIRS LTR/Gypsy
 LTR/Juno LTR/Mag LTR/Ngaro LTR/Pao LTR/TelKA LINE LINE/CR1 LINE/CRE LINE/Dong LINE/I
 LINE/L1 LINE/L2 LINE/Penelope LINE/Proto2 LINE/R1 LINE/R2 LINE/RTE LINE/Rex
 LINE/Soliton SINE SINE/5S SINE/I SINE/R2 SINE/tRNA
 tRNA rRNA Satellite Simple_repeat Low_complexity ARTEFACT RNA Retroposon repeat
);

my @array_to_sort = qw ( DNA DNA/IS5 DNA/P Unknown LTR DNA LINE SINE LTR/Pao RC Unknown );

my $count = 0;
my %position_of;
$position_of{$_} = $count++ for @sort_order;

print Dumper \%position_of;

sub by_order {
  #my @a = split m/\//, $a;
  #my @b = split m/\//, $b;
  $position_of{$a} <=> $position_of{$b};
}

my @new = sort by_order @array_to_sort;

print Dumper \@new;
