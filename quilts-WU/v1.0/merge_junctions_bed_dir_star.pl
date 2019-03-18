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
my $min_sum_length="";
my $min_length=""; 
if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $error=1; } 
if ($ARGV[1]=~/\w/) { $min_sum_length=$ARGV[1];} else { $min_sum_length=40; } 
if ($ARGV[2]=~/\w/) { $min_length=$ARGV[2];} else { $min_length=8; } 

#$dir=~s/\\/\//g;
my $line = "";
my %junctions=();
my %num_observ=();

if ($error==0)
{		
	if (opendir(dir2,"$dir"))
	{
		my @allfolders=readdir dir2;
		closedir dir2;
		foreach my $folder (@allfolders)
		{
			if ($folder=~/^UNCID/)
			{
				if (opendir(dir2,"$dir/$folder"))
				{
					my @allfiles=readdir dir2;
					close dir2; 
					foreach my $filename (@allfiles)
					{
						if ($filename=~/^SJ/ && $filename=~/\.tab$/)
						{
							if(open (IN_RNA, qq!<$dir/$folder/$filename!))
							{
								my $label=$folder;
								while($line = <IN_RNA>)
								{
									chomp($line);
									# Takes one tuple at a time from 1st input file and creates 3 hashes:
									# "num_observ" where key contains chromosome id, beggining of intron, and end of intron, and value is number of times this intron was observed in 1st input file;
									# "num_observ1" where key contains chromosome id, beggining of intron, and end of intron, and value is number of times this intron was observed in 1st input file;
									# "junction1" where key contains chromosome id, beggining of intron, and end of intron, and value is the same as tuple for a particular junction except 4th column
									#	(junction name has an extention specifing the file where tuple was taken from); 
									if ($line=~ /^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
									{
										my $chromosome= $1;
										my $begin_intron=$2-1;
										my $end_intron=$3-1;
										my $num_observ=$7+$8;
										my $overhang=$9; 
										my $begin_exon_1= $begin_intron-$overhang-1;
										my $begin_exon_2=$end_intron+1; 
										my $end_exon_2= $end_intron+$overhang;
										my $length_exons = "$overhang, $overhang";
										my $begin_exons = "$begin_exon_1,$begin_exon_2";
										my @length_exon = split(',', $length_exons,2);
										my @begin_exon = split(',', $begin_exons,2);
										if ($length_exon[0]>=$min_length and $length_exon[1]>=$min_length and $length_exon[0]+$length_exon[1]>=$min_sum_length)
										{
											$num_observ{qq!$chromosome#$begin_intron#$end_intron!}=+$num_observ;
											$junctions{qq!$chromosome#$begin_intron#$end_intron!}.=qq!$label,!;
										}
									}
								}
								close(IN_RNA);
							}
							else {print "error opening STAR junctions file input";}
						}
					}
				
					if (open (OUT, qq!>$dir/merged-junctions.bed!))
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
								my $begin_exon_1=$begin_intron-25;
								my $end_exon_2=$end_intron+25+1;
								my $start_exon_2=$end_intron+1-$begin_exon_1;
								print OUT qq!chr$chromosome\t$begin_exon_1\t$end_exon_2\tj$count\t$num_observ{$key}\t+\t$begin_exon_1\t$end_exon_2\t0\t2\t25,25\t0,$start_exon_2\n!;
								$count++;
							}
						}
						close(OUT);
					}
					else
					{
						print "error opening merged-junctions file";
					}
				}
			}
		}
	}
}
