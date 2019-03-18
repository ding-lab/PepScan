#!c:/perl/bin/perl.exe
#
use strict;
use File::Spec;

# (rjm)
use File::Path(qw(make_path));

my $script_path = File::Spec->rel2abs(__FILE__);
my $script_dir=$script_path;
my $version=$script_path;
$script_dir=~s/^(.*)\/[^\/]+$/$1/;
$version=~s/^.*\/([^\/]+)\/[^\/]+$/$1/;
print "version=$version\n";

my $error=0;
my $type="";
my @dirs=("altsplice_transcripts","other","variants","frameshift");
my %dirs=("altsplice_transcripts"=>"alternative","other"=>"other","variants"=>"variants","frameshift"=>"frameshift",);

#my $root="/ifs/data/proteomics/tcga";
# (rjm)
my $root="/gscuser/rmashl/gsc/working/cptac/test_quilts";


my %stat=();
my $max_length=0;

if ($ARGV[0]=~/\w/) { $type=$ARGV[0];} else { $type="breast"; }

if ($error==0)
{
	opendir(DIR,"$root/samples/$type");
	my $datetime=GetDateTime();
	my @samples = readdir DIR;
	closedir DIR;
	my @samples_ = ();
	my %dirs_ = ();
	foreach my $sample (sort @samples)
	{
		if ($sample=~/\w/)
		{
			if (opendir(DIR,"$root/samples/$type/$sample/protein/fasta"))
			{
				my @versions = readdir DIR;
				closedir DIR;
				@samples_=(@samples_,$sample);
				print qq!Publishing $sample\n!;
				foreach my $version (@versions)
				{
					if ($version=~/\w/)
					{
					    #mkdir("$root/samples/$type/archive/$version-$datetime");
					#	mkdir("$root/samples/$type/archive/$version-$datetime-$dirs{$dir}");
# (rjm)
					    &make_path("$root/samples/$type/archive/$version-$datetime");
					    &make_path("$root/samples/$type/archive/$version-$datetime-$dirs{$dir}");

						foreach my $dir (@dirs)
						{
							system("cp $root/samples/$type/$sample/protein/fasta/$version/$dir/proteome.fasta $root/samples/$type/archive/$version-$datetime-$dirs{$dir}/$sample-$dirs{$dir}.fasta");
							$dirs_{"$version-$datetime"}=1;
						}
					}
				}
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
	
