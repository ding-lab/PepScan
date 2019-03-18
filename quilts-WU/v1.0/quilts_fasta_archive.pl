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


if ($ARGV[0]=~/\w/) { $type=$ARGV[0];} else { $error=1; }

if ($error==0)
{
    #mkdir("$root/samples/$type/archive");
    # (rjm)
    &make_path("$root/samples/$type/archive");

	my $datetime=GetDateTime();
#	mkdir("$root/samples/$type/archive/fasta-$datetime");
# (rjm)
    &make_path("$root/samples/$type/archive/fasta-$datetime");

	opendir(DIR,"$root/samples/$type");
	my @samples = readdir DIR;
	closedir DIR;
	my @samples_ = ();
	foreach my $sample (sort @samples)
	{
		if ($sample=~/\w/)
		{
			if (opendir(DIR,"$root/samples/$type/$sample/protein/fasta"))
			{
				closedir DIR;
				system(qq!mv $root/samples/$type/$sample/protein/fasta $root/samples/$type/archive/fasta-$datetime/$sample!);
				#mkdir(qq!$root/samples/$type/$sample/protein/fasta!);
# (rjm)
				&make_path("$root/samples/$type/$sample/protein/fasta");
			}
		}
	}
}

sub GetDateTime
{
	my $sec="";
	my $min="";
	my $hour="";
	my $mday="";
	my $mon="";
	my $year="";

	($sec,$min,$hour,$mday,$mon,$year) = localtime();

	if ($sec<10) { $sec="0$sec"; }
	if ($min<10) { $min="0$min"; }
	if ($hour<10) { $hour="0$hour"; }
	if ($mday<10) { $mday="0$mday"; }
	$mon++;
	if ($mon<10) { $mon="0$mon"; }
	$year+=1900;
	my $date="$year$mon$mday-$hour$min$sec";
	
	return $date;
}
