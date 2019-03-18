#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $type="";

#my $root="/ifs/data/proteomics/tcga";
# (rjm)
my $root="/gscuser/rmashl/gsc/working/cptac/test_quilts";

my $level1="";

if ($ARGV[0]=~/\w/) { $type=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $level1=$ARGV[1];} else { $error=1; }

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
				print qq!$root/samples/$type/$sample/$level1\n!;
				system(qq!rm -R $root/samples/$type/$sample/$level1!);
			}
		}
	}
}
