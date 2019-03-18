#!c:/perl/bin/perl.exe
#
# Step#1: ensembl_gtf_to_bed.pl
# Step#2: merge_junctions.pl
# Step#3: filter_known_transcripts.pl
# Step#4: filter_A.pl
# Step#5: protein_seq_from_genome_using_bed.pl
# Step#6: protein_seq_from_genome_using_junctions.pl

# Script takes Homo_sapiens.GRCh37.68(69).gtf file (or .gff Version 2 file ('Gene-Finding Format' or 'General Feature Format')) which contains features (tuple) in either
# 	complete gene or RNA transcript or protein structures;
# Structure of GTF/GFF file: <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
#	more information on file format (http://useast.ensembl.org/info/website/upload/gff.html);
# Script makes 3 output files:
# 1) Homo_sapiens.GRCh37.68.gtf.rna.bed (transcripts with features of type exon)
# 2) Homo_sapiens.GRCh37.68.gtf.prot.bed (proteins with features of type CDS (coding regions of the seq))
# 3) Homo_sapiens.GRCh37.68.gtf.ensembl_gtf_to_bed.log (usually empty -> no errors)

use strict;
my $filename="";

my $error=0;

if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }

if ($error==0)
{
	open (OUT,">$filename.prot.bed");
	open (OUT_RNA,">$filename.rna.bed");
	open (LOG,">$filename.ensembl_gtf_to_bed.log");
	my $line="";
	my %loc_prot=();
	my %loc_rna=();
	my %gene=();
	if (open (IN,"$filename"))
	{
		# script goes over parses *.gtf file and creates two hashes:
		# 1) "gene" with type1 key = protein name and type2 key = transcript name, value contains gene name;
		# 2) "loc" with type1 key = all protein names and type2 key = all transcript names, value contains chromosome, strand direction,
		# 	beginning and end of coding region of a sequence;
		while ($line=<IN>)
		{
			chomp($line);
			
			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
			{
				my $chr=$1;
				my $source=$2;
				my $feature=$3;
				my $start=$4;
				my $end=$5;
				my $score=$6;
				my $strand=$7;
				my $frame=$8;
				my $attribute=$9;
				if ($chr=~/^\s*[0-9XY]+$\s*/ and $feature=~/^\s*CDS\s*/)
				{
					if ($attribute=~/protein_id\s*\"([^\"]+)\"\s*;/)
					{
						my $protein=$1;
						if ($attribute=~/gene_name\s*\"([^\"]+)\"\s*;/)
						{
							my $gene=$1;
							$gene{$protein}=$gene;
						} 
						else
						{
							$gene{$protein}=""; print LOG qq!Gene not found: $line\n!;
						}
						$loc_prot{$protein}.="#$chr$strand:$start-$end#";
					} else
					{
						print LOG qq!Protein not found: $line\n!;
					}
				}
				if ($chr=~/^\s*[0-9XY]+$\s*/ and $feature=~/^\s*exon\s*/)
				{
					if ($attribute=~/transcript_id\s*\"([^\"]+)\"\s*;/)
					{
						my $rna=$1;
						if ($attribute=~/gene_name\s*\"([^\"]+)\"\s*;/)
						{
							my $gene=$1;
							$gene{$rna}=$gene;
						} else
						{
							$gene{$rna}=""; print LOG qq!Gene not found: $line\n!;
						}
						$loc_rna{$rna}.="#$chr$strand:$start-$end#";
					} else
					{
						print LOG qq!transcript not found: $line\n!;
					}
				}
			}
		}
		close(IN);
		
		my $count=0;
		foreach my $protein (sort keys %loc_prot)
		{
			my $chr="";
			my $strand="";
			my $chr_same=1;
			my $strand_same=1;
			my $exon_num=0;
			my $min=100000000000;
			my $max=0;
			my @exon_pos=();
			my @exon_length=();
			my $temp=$loc_prot{$protein};
			while($temp=~s/^#([^#]+)#//)
			{
				my $temp_=$1;
				if ($temp_=~/^([0-9XY]+)([\+\-])\:([0-9]+)\-([0-9]+)$/)
				{
					my $chr_=$1;
					my $strand_=$2;
					my $start=$3;
					my $end=$4;
					if ($min>$start) { $min=$start; }
					if ($max<$start) { $max=$start; }
					if ($min>$end) { $min=$end; }
					if ($max<$end) { $max=$end; }
					if ($chr!~/\w/)
					{
						$chr=$chr_;
					} else
					{
						if ($chr!~/^$chr_$/)
						{
							$chr_same=0;
						}
					}
					if ($strand!~/\w/)
					{
						$strand=$strand_;
					}  else
					{
						if ($strand!~/^$strand_$/)
						{
							$strand_same=0;
						}
					}
				} else
				{
					print LOG qq!Error parsing location: $protein: $temp_\n!;
				}
			}
			if ($chr_same==0 or $strand_same==0)
			{
				print LOG qq!Error: $loc_prot{$protein} $chr_same $strand_same\n!;
			}
			else
			{
				my $temp=$loc_prot{$protein};
				if ($strand=~/\+/)
				{
					while($temp=~s/^#([^#]+)#//)
					{
						my $temp_=$1;
						if ($temp_=~/^([0-9XY]+)([\+\-])\:([0-9]+)\-([0-9]+)$/)
						{
							my $start=$3;
							my $end=$4;
							if ($temp!~/\w/) { $end+=3; }
							$exon_pos[$exon_num]=$start-$min;
							$exon_length[$exon_num]=$end-$start+1;
							$exon_num++;
						}
					}
					$min--;
					$max+=3;
				}
				else
				{
					$min-=4;
					while($temp=~s/#([^#]+)#$//)
					{
						my $temp_=$1;
						if ($temp_=~/^([0-9XY]+)([\+\-])\:([0-9]+)\-([0-9]+)$/)
						{
							my $start=$3;
							my $end=$4;
							if ($exon_num==0) { $start-=4; } else { $start--; }
							$end--;
							$exon_pos[$exon_num]=$start-$min;
							$exon_length[$exon_num]=$end-$start+1;
							$exon_num++;
						}
					}
				}
				print OUT qq!chr$chr\t$min\t$max\t$protein-$gene{$protein}\t1000\t$strand\t$min\t$max\t0\t$exon_num!;
				print OUT qq!\t!;
				for(my $i=0; $i<$exon_num;$i++) { if ($i>0) { print OUT ","; } print OUT "$exon_length[$i]"; }
				print OUT qq!\t!;
				for(my $i=0; $i<$exon_num;$i++) { if ($i>0) { print OUT ","; } print OUT "$exon_pos[$i]"; }
				print OUT qq!\n!;
				$count++;
			}
		}
		print LOG "$count proteins\n";
		
		$count=0;
		foreach my $rna (sort keys %loc_rna)
		{
			my $chr="";
			my $strand="";
			my $chr_same=1;
			my $strand_same=1;
			my $exon_num=0;
			my $min=100000000000;
			my $max=0;
			my @exon_pos=();
			my @exon_length=();
			my $temp=$loc_rna{$rna};
			while($temp=~s/^#([^#]+)#//)
			{
				my $temp_=$1;
				if ($temp_=~/^([0-9XY]+)([\+\-])\:([0-9]+)\-([0-9]+)$/)
				{
					my $chr_=$1;
					my $strand_=$2;
					my $start=$3;
					my $end=$4;
					if ($min>$start) { $min=$start; }
					if ($max<$start) { $max=$start; }
					if ($min>$end) { $min=$end; }
					if ($max<$end) { $max=$end; }
					if ($chr!~/\w/)
					{
						$chr=$chr_;
					} else
					{
						if ($chr!~/^$chr_$/)
						{
							$chr_same=0;
						}
					}
					if ($strand!~/\w/)
					{
						$strand=$strand_;
					}  else
					{
						if ($strand!~/^$strand_$/)
						{
							$strand_same=0;
						}
					}
				} 
				else
				{
					print LOG qq!Error parsing location: $rna: $temp_\n!;
				}
			}
			if ($chr_same==0 or $strand_same==0)
			{
				print LOG qq!Error: $loc_rna{$rna} $chr_same $strand_same\n!;
			}
			else
			{
				my $temp=$loc_rna{$rna};
				if ($strand=~/\+/)
				{
					while($temp=~s/^#([^#]+)#//)
					{
						my $temp_=$1;
						if ($temp_=~/^([0-9XY]+)([\+\-])\:([0-9]+)\-([0-9]+)$/)
						{
							my $start=$3;
							my $end=$4;
							if ($temp!~/\w/) { $end+=3; }
							$exon_pos[$exon_num]=$start-$min;
							$exon_length[$exon_num]=$end-$start+1;
							$exon_num++;
						}
					}
					$min--;
					$max+=3;
				}
				else
				{
					$min-=4;
					while($temp=~s/#([^#]+)#$//)
					{
						my $temp_=$1;
						if ($temp_=~/^([0-9XY]+)([\+\-])\:([0-9]+)\-([0-9]+)$/)
						{
							my $start=$3;
							my $end=$4;
							if ($exon_num==0) { $start-=4; } else { $start--; }
							$end--;
							$exon_pos[$exon_num]=$start-$min;
							$exon_length[$exon_num]=$end-$start+1;
							$exon_num++;
						}
					}
				}
				print OUT_RNA qq!chr$chr\t$min\t$max\t$rna-$gene{$rna}\t1000\t$strand\t$min\t$max\t0\t$exon_num!;
				print OUT_RNA qq!\t!;
				for(my $i=0; $i<$exon_num;$i++) { if ($i>0) { print OUT_RNA ","; } print OUT_RNA "$exon_length[$i]"; }
				print OUT_RNA qq!\t!;
				for(my $i=0; $i<$exon_num;$i++) { if ($i>0) { print OUT_RNA ","; } print OUT_RNA "$exon_pos[$i]"; }
				print OUT_RNA qq!\n!;
				$count++;
			}
		}
		print LOG "$count transcripts\n";
	}
}