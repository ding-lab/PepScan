#!c:/perl/bin/perl.exe
#
use strict;

#my $root="/ifs/data/proteomics/tcga/";
# (rjm)
my $root="/gscuser/rmashl/gsc/working/cptac/test_quilts";


my %title=('colon'=>'Appendix B. Colon','breast'=>'Appendix C. Breast','ovarian'=>'Appendix D. Ovarian',);
my @test_title=("DNA (germline)","DNA (tumor)", "RNA (tumor)", "Protein (tumor)");
my @test_dir=("dna-germline","dna", "rna", "protein");
my @test_title1=("fastq dwnld,bam dwnld, vcf dwnld", "fastq dwnld,bam dwnld,junctions", "fasta,peptides,bed");
my @test_dir1=("fastq,bam,vcf","fastq,bam,vcf", "fastq,bam,junctions", "fasta,mzident,bed");
my @test_ext=(".fastq,.bam,.vcf",".fastq,.bam,.vcf", ".fastq,.bam,-junctions.bed", ".xml,.mzident,.bed");
my @test_count=(3,3,3,3);

if (opendir(dir,"$root/samples"))
{
	my @allfiles=readdir dir;
	closedir dir;
	foreach my $type (@allfiles)
	{
		if ($type!~/^template$/i)
		{
			if ($type=~/\w/)
			{
				if ($type=~/xenograft/)
				{
					@test_title=("DNA (germline)","DNA (tumor)", "RNA (tumor)", "Protein (tumor)");
					@test_dir=("dna-germline","dna", "rna", "protein");
					@test_title1=("fastq dwnld,bam dwnld, vcf dwnld", "fastq dwnld,bam dwnld,junctions", "fasta,peptides,bed");
					@test_dir1=("fastq,bam,vcf","fastq,bam,vcf", "fastq,bam,junctions", "fasta,peptides,bed");
					@test_ext=(".fastq,.bam,.vcf",".fastq,.bam,.vcf", ".fastq,.bam,-junctions.bed", ".xml,.psm,.bed");
					@test_count=(3,3,3,3);
				}
				else
				{
					@test_title=("DNA (GL)", "DNA", "RNA", "Protein");
					@test_dir=("dna-germline","dna", "rna", "protein");
					@test_title1=("bam dwnld","bam dwnld,vcf dwnld", "bam dwnld,tophat,star", "fasta,peptides,bed");
					@test_dir1=("bam", "bam,vcf/TCGA-20130502", "bam,tophat/junction_tophat208bowtie2_gMx,STAR", "fasta,peptides,bed");
					@test_ext=(".gto", ".gto,.vcf", ".gto,.bed,.bed", ".xxml,.psm,.bed");
					@test_count=(1,2,3,3);
				}
				if (opendir(dir,"$root/samples/$type"))
				{
					my @allfiles1=readdir dir;
					closedir dir;
					if(open(OUT,">$root/samples/$type/$type.html"))
					{
						my $temp=$type;
						if ($title{$type}=~/\w/) { $temp=$title{$type}; }
						#open(OUT_SCRIPT,">$root/samples/$type/$type.sh");
						print OUT qq!<html><head></head><body><h2>$temp</h2><table border="1">\n!;
						print OUT qq!<tr>!;
						print OUT qq!<td><b>Sample</b></td>!;
						for(my $i=0;$test_title[$i]=~/\w/;$i++)
						{
							print OUT qq!<td colspan="$test_count[$i]"><b>$test_title[$i]</b></td>!;
						}
						print OUT qq!</tr>!;
						print OUT qq!<tr>!;
						print OUT qq!<td><b></b></td>!;
						for(my $i=0;$test_title[$i]=~/\w/;$i++)
						{
							my $temp="$test_title1[$i],";
							while($temp=~s/^([^\,]+)\,//)
							{
								print OUT qq!<td><b>$1</b></td>!;
							}
						}
						print OUT qq!</tr>!;
						foreach my $sample (sort @allfiles1)
						{
							if ($sample=~/^whim/i or $sample=~/^TCGA/i)
							{
								if (opendir(dir,"$root/samples/$type/$sample"))
								{
									closedir dir;
									print OUT qq!<tr>!;
									print OUT qq!<td><b>$sample</b></td>!;
									my $first=0;
									for(my $i=0;$test_title[$i]=~/\w/;$i++)
									{
										my $temp="$test_dir1[$i],";
										my $temp_ext="$test_ext[$i],";
										my @done=();
										my $j=0;
										while($temp=~s/^([^\,]+)\,//)
										{
											my $test_dir1=$1;
											if ($temp_ext=~s/^([^\,]+)\,//)
											{
												my $test_ext=$1;
												my $count=0;
												if (opendir(dir,"$root/samples/$type/$sample/$test_dir[$i]/$test_dir1"))
												{
													my @allfiles=readdir dir;
													closedir dir;
													foreach my $filename (@allfiles)
													{
														if ($filename=~/$test_ext$/)
														{
															$count++;
														}
													}
												}
												if ($count==0) { $done[$j]=0; }
												else { $done[$j]=1; }
											}
											$j++;
										}
										for(my $j=0;$j<$test_count[$i];$j++)
										{
											if ($done[$j]==0) 
											{
												my $passed=0;
												for(my $k=$j+1;$k<$test_count[$i];$k++)
												{
													if ($done[$k]!=0) 
													{
														$passed=1;
													}
												}
												if ($passed==0)
												{
													if ($first==0)
													{
														$first=1;
														my $found=0;
														if ($test_title[$i]=~/^Protein/i and $j==0)
														{
															if ($type=~/xenograft/ or $type=~/test/)
															{
																if (open(TEST,"$root/samples/$type/$sample/$test_dir[$i]/fasta/create.sh"))
																{
																	close(TEST);
																	$found=1;
																	print OUT qq!<td bgcolor="#FFA500">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>!; 
																}
																else
																{
																	#if (open(OUT_SCRIPT_,">$root/samples/$type/$sample/$test_dir[$i]/fasta/create.sh"))
																	#{
																		#print OUT_SCRIPT qq!qsub $root/samples/$type/$sample/$test_dir[$i]/fasta/create.sh\n!;
																		#print OUT_SCRIPT_ qq!perl $root/scripts/quilts/v0.8/quilts_xenograft.pl $type/$sample\n!;
																		#close(OUT_SCRIPT_);
																	#}
																}
															}
														}
														if ($found==0)
														{
															print OUT qq!<td bgcolor="#000000">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>!; 
															#print OUT qq!<td bgcolor="#FF0000">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>!; 
														}
													}
													else
													{
														print OUT qq!<td bgcolor="#000000">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>!; 
													}
												}
												else
												{
													print OUT qq!<td bgcolor="#000000">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>!; 
													#print OUT qq!<td bgcolor="#999999">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>!; 
												}
											}
											else 
											{
												print OUT qq!<td bgcolor="#00FF00">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</td>!; 
											}
										}
									}
									print OUT qq!</tr>!;
								}
							}
						}
						print OUT qq!</table></body>\n!;
						close(OUT);
						#close(OUT_SCRIPT);
						#system(qq!chmod u+x $root/samples/$type/$type.sh!);
						#system(qq!$root/samples/$type/$type.sh!);
					}
				}
			}
		}
	}
}
	
	
