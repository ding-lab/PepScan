#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my @dirs=(
	"Protected_Mutations_COAD/UCSC__IlluminaGA_DNASeq_Cont/Level_2",
	"Protected_Mutations_COAD/UCSC__SOLiD_DNASeq_Cont/Level_2",
	
	
	"Protected_Mutations_READ/BCM__SOLiD_DNASeq_Cont/Level_2",
	"Protected_Mutations_READ/BI__IlluminaGA_DNASeq_Cont/Level_2",
	"Protected_Mutations_READ/BI__SOLiD_DNASeq_Cont/Level_2",
	"Protected_Mutations_READ/UCSC__IlluminaGA_DNASeq_Cont/Level_2",
	"Protected_Mutations_READ/UCSC__SOLiD_DNASeq_Cont/Level_2",
	);

if ($error==0)
{
	foreach my $dir (@dirs)
	{
		print qq!$dir\n!;
		opendir(DIR,"$dir");
		my @files = readdir DIR;
		closedir DIR;
		foreach my $filename (@files)
		{
			if ($filename=~/(TCGA-[^-_]+-[^-_]+)[-_].*\.vcf$/i)
			{
				my $sample=$1;
				system(qq!cp $dir/$filename $sample/dna/vcf!);
				print qq!$dir/$filename\n!;
			}
		}
	}
}
