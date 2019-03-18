#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $type="";

#my $root="/ifs/data/proteomics/tcga";
# (rjm)
my $root="/gscuser/rmashl/gsc/working/cptac/test_quilts";


my $dir1="";
my $dir2="";
my $ext="";

if ($ARGV[0]=~/\w/) { $type=$ARGV[0];} else { $type="breast"; }
if ($ARGV[1]=~/\w/) { $dir1=$ARGV[1];} else { $dir1="dna/vcf"; }
if ($ARGV[2]=~/\w/) { $dir2=$ARGV[2];} else { $dir2="dna/vcf/TCGA-20130502"; }
if ($ARGV[3]=~/\w/) { $ext=$ARGV[3];} else { $ext="vcf"; }

if ($error==0)
{
	opendir(DIR,"$root/samples/$type");
	my @samples = readdir DIR;
	closedir DIR;
	foreach my $sample (sort @samples)
	{
		if ($sample=~/\w/)
		{
			if ($sample=~/^whim/i or $sample=~/^TCGA/i)
			{
				print(qq!mv $root/samples/$type/$sample/$dir1/*.$ext $root/samples/$type/$sample/$dir2\n!);
				system(qq!mv $root/samples/$type/$sample/$dir1/*.$ext $root/samples/$type/$sample/$dir2!);
			}
		}
	}
}
