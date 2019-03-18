#!/usr/local/bin/perl
#
# Step#1: ensembl_gtf_to_bed.pl
# Step#2: merge_junctions.pl
# Step#3: filter_known_transcripts.pl
# Step#4: filter_A.pl
# Step#5: protein_seq_from_genome_using_bed.pl
# Step#6: protein_seq_from_genome_using_junctions.pl

# This script merges 2 *-junctions.bed files (*_L005-junctions.bed with *_L007-junctions.bed, or *_L006-junctions.bed with *_L008-junctions.bed);
# The description of the structure for the tuple (a horizontal row of data items) in input files could be found here: http://genome.ucsc.edu/FAQ/FAQformat.html#format1
# The output files (4):
# 1) *-num-observed.txt that contains count for each junction from both input files;
# 2) *-merged-junctions.bed that has tuples for junctions that are in both input files as well as in each input file separetely;
# 3) *inputfile#1-only-junctions.bed that has tuples for junctions that are only in the 1st input file;
# 4) *inputfile#2-only-junctions.bed that has tuples for junctions that are only in the 2nd input files

use strict;
use Data::Dumper;

my $error=0;
my $dir="";
my $dir2="";
my $min_sum_length="";
my $min_length="";

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $error=1; } 
if ($ARGV[1]=~/\w/) { $dir2=$ARGV[1];} else { $dir2=$dir; } 
if ($ARGV[2]=~/\w/) { $min_sum_length=$ARGV[2];} else { $min_sum_length=45; } 
if ($ARGV[3]=~/\w/) { $min_length=$ARGV[3];} else { $min_length=8; } 

$dir=~s/\\/\//g;
my $line = "";
my %junctions=();
my %num_observ=();

if ($error==0)
{		
	if (opendir(dir,"$dir"))
	{
		my @allfiles=readdir dir;
		closedir dir;
		foreach my $filename (@allfiles)
		{
			if (($filename=~/\-junctions\.bed$/ or $filename=~/^junctions\.bed$/) and $filename!~/merged-junctions.bed/)
			{
				if(open (IN_RNA, qq!$dir/$filename!))
				{
					my $label=$filename;
					$label=~s/^.*[_\-]([^_\-]+)\-junctions\.bed$/$1/i;
					while($line = <IN_RNA>)
					{
						chomp($line);
						# Takes one tuple at a time from 1st input file and creates 3 hashes:
						# "num_observ" where key contains chromosome id, beggining of intron, and end of intron, and value is number of times this intron was observed in 1st input file;
						# "num_observ1" where key contains chromosome id, beggining of intron, and end of intron, and value is number of times this intron was observed in 1st input file;
						# "junction1" where key contains chromosome id, beggining of intron, and end of intron, and value is the same as tuple for a particular junction except 4th column
						#	(junction name has an extention specifing the file where tuple was taken from); 
						if ($line=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
						{
							my $chromosome= $1;
							my $begin_exon_1= $2;
							my $end_exon_2= $3;
							my $junc_number = $4;
							my $num_observ = $5; 
							my $length_exons = $11;
							my $begin_exons = $12;
							my @length_exon = split(',', $length_exons,2);
							my @begin_exon = split(',', $begin_exons,2);
							my $begin_intron = $begin_exon_1 + $length_exon[0] + 1;
							my $end_intron = $begin_exon_1 + $begin_exon[-1] - 1;
							if ($length_exon[0]>=$min_length and $length_exon[1]>=$min_length and $length_exon[0]+$length_exon[1]>=$min_sum_length)
							{
								$num_observ{qq!$chromosome#$begin_intron#$end_intron!}=+$num_observ;
								$junctions{qq!$chromosome#$begin_intron#$end_intron!}.=qq!$junc_number-$label,!;
							}
						}
					}
					close(IN_RNA);
				}
			}
		}
	
		if (open (OUT, qq!>$dir2/merged-junctions.bed!))
		{
			my $count=0;
			foreach my $key (sort keys %num_observ)
			{
				# Prints out the number of observations for particular junction (into *-num-observ.txt=STAT)if junction is in both input files;
				# Prints out the tuple for a particular junction (into *-merged-junctions.bed=OUT)if junction is in both input files or in any of the two input files,
				# 	and changes junctions' labels (4th column);
				# Variable $both counts how many junctions were observed in both files;
				if ($key=~ /^([^#]+)#([^#]+)#([^#]+)$/)
				{
					my $chromosome = $1;
					my $begin_intron = $2;
					my $end_intron = $3;
					my $begin_exon_1=$begin_intron-50-1;
					my $end_exon_2=$end_intron+50+1;
					my $start_exon_2=$end_intron+1-$begin_exon_1;
					print OUT qq!$chromosome\t$begin_exon_1\t$end_exon_2\tj$count\t$num_observ{$key}\t+\t$begin_exon_1\t$end_exon_2\t0\t2\t50,50\t0,$start_exon_2\n!;
					$count++;
				}
			}
			close(OUT);
		}
	}
}
