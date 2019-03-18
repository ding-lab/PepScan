#!/usr/local/bin/perl
#
# Step#1: ensembl_gtf_to_bed.pl
# Step#2: merge_junctions.pl
# Step#3: filter_known_transcripts.pl
# Step#4: filter_A.pl
# Step#5: protein_seq_from_genome_using_bed.pl
# Step#6: protein_seq_from_genome_using_junctions.pl

# This script takes as input 2 *.bed files: 
# 1. Junctions from RNA-Seq data aligned to the genome (*-merged-junctions.bed) and
# 2. A collection of gene models (e.g. Ensembl) where each tuple for a particular transcript has the structure described here: http://genome.ucsc.edu/FAQ/FAQformat.html#format1

# Two output files: 
# 1. Junctions from RNA-Seq data aligned to the genome that have not matched the gene model (*-junctions.filter.bed);
# 2. Statistics for junctions from RNA-Seq data alighed to the genome that have matched the gene model based on "C", "B", and "E" critiria(*-junctions.filter.stat.txt);

use strict;
use Data::Dumper;

my $error=0;
my $filename1="";
my $filename2="";
if ($ARGV[0]=~/\w/) { $filename1=$ARGV[0];} else { $error=1; } # RNAseq
if ($ARGV[1]=~/\w/) { $filename2=$ARGV[1];} else { $error=1; } # gene models

$filename1=~s/\\/\//g;
$filename2=~s/\\/\//g;
my $filename1_=$filename1;
$filename1_=~s/\.bed//i;
my $dir=$filename1;
if ($dir!~/\/([\/]+)$/) { $dir="."; }
my $line = "";
my $count_model_exons=0;
my $count_model_begin=0;
my $count_model_end=0;
my $chromosome = "";
my %stat = ();
my $num_observ_max=0;
my %junctions_model=();
my $count_junctions=0;
my %count_junctions=();

if ($error==0)
{
	# Goes over the DNA database file(gene models) and identifies introns (begin_intron and end_intron) for each gene;
	# Creates a hash "junctions_model" where key contains chromosome, begin_intron or end_intron, or both begin_intron and end_intron, and value is either "C", "B", or "E":
	# "C" - for beginning and end of each intron;
	# "B" - for the begging of the 1st intron in the transcript, if transcript is "+" strand, or for the end of the last intron, if transcript is "-" strand;
	# "E" - for the end of the last intron in the transcript, if transcript is "+" strand, or for the beginning of the first intron if transcript is "-" strand; 
	if(open (IN_DNA, qq!$filename2!))
	{
		my $line_number = 0;
		while($line = <IN_DNA>)
		{
			chomp($line);
			if ($line=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
			{
				$chromosome = $1;
				my $begin_exon = $2;
				my $protein = $4;
				my $strand = $6;
				my $block_count = $10;
				my $block_sizes = $11;
				my $block_starts = $12; 
				my @sizes = split(',',$block_sizes);
				my @starts = split(',',$block_starts);
				my $start = $begin_exon;
				 
				for (my $i = 0; $i < $block_count-1; $i++)
				{
					my $begin_intron = $start + $sizes[$i] + 1;
					my $end_intron = $begin_exon + $starts[$i+1] - 1; 
					$junctions_model{qq!$chromosome#$begin_intron#$end_intron!}="C";
					if ($i==0) 
					{ 
						if ($strand=~/^\+$/) { $junctions_model{qq!$chromosome#$begin_intron!}="B"; $count_model_begin++; }
						else { $junctions_model{qq!$chromosome#$begin_intron!}="E"; $count_model_end++; } 
					}
					if ($i==($block_count-1-1)) 
					{ 
						if ($strand=~/^\-$/) { $junctions_model{qq!$chromosome#$end_intron!}="B"; $count_model_begin++; }
						else { $junctions_model{qq!$chromosome#$end_intron!}="E"; $count_model_end++; } 
					}
					$start = $begin_exon + $starts[$i+1];
					$count_model_exons++;
				}
			}
			$line_number++;
		}
		close(IN_DNA);
	}
	
	# Goes over RNAseq data (in *-merged-junctions.bed), and assignes status ("C", "B", or "E") for each intron (junction) if it has been matched to the database intron;
	# If no match was found for RNAseq intron, the intron tuple is saved in the *-junctions.filter.bed file;
	open (OUT, qq!>$filename1_.filter.bed!);
	if(open (IN_RNA, qq!$filename1!))
	{
		while($line = <IN_RNA>)
		{
			chomp($line);
			if ($line=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
			{
				$chromosome= $1;
				my $begin_exon_1= $2;
				my $end_exon_2= $3;
				my $junc_number = $4;
				my $num_observ = $5; 
				my $length_exons = $11;
				my $begin_exons = $12;
				my @length_exon = split(',', $length_exons,2);
				my @begin_exon = split(',', $begin_exons,2);
				my $size = ();
				my $non_matches = 0;
				my $junction_status="";
				# variable $num_observ_max finds the number of maximum observations for a particular intron in RNAseq data;
				if ($num_observ_max<$num_observ) { $num_observ_max=$num_observ; }
				if ($num_observ >= 1) 
				{	
					my $begin_intron = $begin_exon_1 + $length_exon[0] + 1;
					my $end_intron = $begin_exon_1 + $begin_exon[-1] - 1;
					my $length_intron = $end_intron - $begin_intron;
					if ($junctions_model{qq!$chromosome#$begin_intron#$end_intron!}=~/\w/)
					{
						$junction_status="C";
					}
					else
					{
						if ($junctions_model{qq!$chromosome#$begin_intron!}=~/\w/)
						{
							if ($junction_status!~/\w/) { $junction_status=$junctions_model{qq!$chromosome#$begin_intron!}; }
						}
						else
						{
							if ($junctions_model{qq!$chromosome#$end_intron!}=~/\w/)
							{
								if ($junction_status!~/\w/) { $junction_status=$junctions_model{qq!$chromosome#$end_intron!}; }
							}
						}
					}
					if ($junction_status!~/\w/) 
					{
						print OUT qq!$line\n!;
					}
					else
					{
						$stat{"$junction_status#$num_observ"}++;
						$count_junctions{"$junction_status"}++;
					}
					$count_junctions++;
				}
			}
		}
		close(IN_RNA);
	}
	close(OUT);
	
	# This file contains statistics for matched introns from RNAseq data with status "C", "B", and "E", where
	#	num_observ is the number of times intron was observed during RNA sequencing,
	#	C, B, and E columns state the number of introns for "C", "B" and "E" matches;
	if (open(OUT, qq!>$filename1_.filter.stat.txt!))
	{
		print OUT qq!num_observ\tC_matches\tB_matches\tE_matches\n!;
		for(my $num_observ=1;$num_observ<=$num_observ_max;$num_observ++)
		{
			print OUT qq!$num_observ!;
			foreach my $type ("C","B","E")
			{
				if ($stat{"$type#$num_observ"}!~/\w/) { $stat{"$type#$num_observ"}=0;}
				print OUT qq!\t$stat{"$type#$num_observ"}!;
			}
			print OUT qq!\n!;
		}
		close(OUT);
	}

	print qq!Model exons: $count_model_exons\n!;
	print qq!Model transcript begin: $count_model_begin\n!;
	print qq!Model transcript end: $count_model_end\n!;
	print qq!RNA-Seq introns: $count_junctions\n!;
	foreach my $type ("C","B","E",)
	{
		print qq!RNA-Seq introns $type: $count_junctions{$type}\n!;
	}
}

sub evaluate {
    (my $val, my @array) = @_;
    foreach (@array)
    {
        if ($val == $_)
        {
            return $_;
        }
    }
}
