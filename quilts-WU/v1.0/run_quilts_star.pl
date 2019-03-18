#!c:/perl/bin/perl.exe
#
use strict;
use File::Spec;
my $script_dir = File::Spec->rel2abs(__FILE__);
my $version=$script_dir;
$version=~s/^.*\/([^\/]+)\/[^\/]+$/$1/;
print "version=$version\n";

#my $root="/ifs/data/proteomics/tcga";
# (rjm)
my $root="/gscuser/rmashl/gsc/working/cptac/test_quilts";


#my $scripts="$root/scripts/quilts/$version";
# (rjm)
my $scripts="$root/scripts/quilts-WU/$version";


my $type_selected="breast";


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
				if ($type=~/^$type_selected$/)
				{
					if (opendir(dir,"$root/samples/$type"))
					{
						my @allfiles1=readdir dir;
						closedir dir;
						foreach my $sample (sort @allfiles1)
						{
							my $vcf_check=0;
							my $star_check=0; 
							if ($sample=~/^whim/i or $sample=~/^TCGA/i)
							{
								if (opendir(dir,"$root/samples/$type/$sample"))
								{
									closedir dir;
									open(OUT_SH,">$root/samples/$type/$sample/quilts.$type.$sample.sh");
									#print OUT_SH "#\!/bin/bash\n#\$ -cwd\n";
									if (opendir(dir,"$root/samples/$type/$sample/dna/vcf/TCGA-20130502")) #check vcf file 
									{
										my @allfiles=readdir dir;
										closedir dir;
										foreach my $filename (@allfiles)
										{
											if ($filename=~/.vcf$/ and $filename!~/\.germline\.vcf/ and $filename!~/\.somatic_only\.vcf/) {$vcf_check=1;}
										}
									}
									if (opendir(dir,"$root/samples/$type/$sample/rna/STAR")) #check junction file
									{
										my @allfolders=readdir dir;
										closedir dir;
										foreach my $folder (@allfolders)
										{
											if (opendir(dir,"$root/samples/$type/$sample/rna/STAR/$folder")) #check junction file
											{
												my @allfiles=readdir dir;
												close dir;
												foreach my $filename (@allfiles)
												{
													if ($filename=~/^SJ/) {$star_check=1;}
												}
											}
										}
									}
									if (open(IN,"$root/samples/$type/$sample/protein/fasta/version.txt"))
									{
										print "$sample DONE vcf:$vcf_check star:$star_check\n";
										close(IN);
									}
									else
									#if (($star_check==1) && ($vcf_check==1))
									{
										print OUT_SH qq!perl $scripts/quilts_tcga_star.pl $type/$sample\n!;
										#system(qq!qsub -pe perl $scripts/quilts_tcga_star.pl $root/samples/$type/$sample!);
										print "$sample submitted vcf:$vcf_check star:$star_check\n";
										close(OUT_SH);
										system(qq!chmod u+x $root/samples/$type/$sample/quilts.$type.$sample.sh!);
										system(qq!qsub -cwd -S /bin/bash $root/samples/$type/$sample/quilts.$type.$sample.sh!);
									}
								}
							}
						}
					}
				}
			}	
		}
	}
}
	
	
