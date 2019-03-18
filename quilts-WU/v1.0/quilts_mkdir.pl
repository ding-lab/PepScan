
#!c:/perl/bin/perl.exe
#
use strict;

# (rjm)
use File::Path(qw(make_path));

my $error=0;
my $type="";

#my $root="/ifs/data/proteomics/tcga";
# (rjm)
my $root="/gscuser/rmashl/gsc/working/cptac/test_quilts";


my $level1="";
my $level2="";
my $level3="";

if ($ARGV[0]=~/\w/) { $type=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $level1=$ARGV[1];} else { $error=1; }
if ($ARGV[2]=~/\w/) { $level2=$ARGV[2];} else { $level2=""; }
if ($ARGV[3]=~/\w/) { $level3=$ARGV[3];} else { $level3=""; }

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
			    #mkdir(qq!$root/samples/$type/$sample/$level1!);
                            # (rjm)
			    &make_path("$root/samples/$type/$sample/$level1");
				if ($level2=~/\w/)
				{
				    #mkdir(qq!$root/samples/$type/$sample/$level1/$level2!);
				# (rjm)
				    &make_path("$root/samples/$type/$sample/$level1/$level2");
				}
				if ($level3=~/\w/)
				{
				    #mkdir(qq!$root/samples/$type/$sample/$level1/$level2/$level3!);
				# (rjm)
				    &make_path("$root/samples/$type/$sample/$level1/$level2/$level3");
				}
			}
		}
	}
}
