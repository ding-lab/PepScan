#!/usr/bin/perl
use strict;
use File::Spec;
my $script_dir = File::Spec->rel2abs(__FILE__);
my $version=$script_dir;
$version=~s/^.*\/([^\/]+)\/[^\/]+$/$1/;
my $usage = '
This script will run quilts for a given sample directory
';
die $usage unless scalar @ARGV == 1;
my ($sample_dir)=@ARGV;
my $sample= (split(/\//,$sample_dir))[-1];
my $vcf_dir=$sample_dir."/vcf"; 
my $junction_dir=$sample_dir."/junction"; 
my $protein_dir=$sample_dir."/protein"; 
print "version=$version\n";
my $root="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/Quilts/Song_adapted";
my $scripts="$root/scripts/quilts-WU/$version";

if (opendir(dir,$sample_dir))
{
	my @allfiles=readdir dir;
	closedir dir;
	my $vcf_check=0; 
        my $tophat_check=0; 
	foreach my $file (@allfiles)
	{								
	open(OUT_SH,">$sample_dir/quilts.$sample.sh");
	opendir(dir,$vcf_dir); 
	my @allfiles=readdir dir;
	closedir dir;
	foreach my $filename (@allfiles)
	 {
	 if ($filename=~/.vcf$/ and $filename!~/\.germline\.vcf/ and $filename!~/\.somatic_only\.vcf/) {$vcf_check=1;}
         }
        }
        if (open(IN,"$junction_dir/junctions.bed")) #check junction file
	{
	  close(IN);
	  $tophat_check=1;
        }
        if (open(IN,"$protein_dir/fasta/version.txt"))
	{
	  print "$sample DONE vcf:$vcf_check tophat:$tophat_check\n";
	  close(IN);
	}
	else
	{
	 print OUT_SH qq!perl $scripts/quilts_tcga_tophat.pl $sample_dir\n!;
         print "$sample submitted vcf:$vcf_check tophat:$tophat_check\n";
	 close(OUT_SH);
	 system(qq!chmod u+x $sample_dir/quilts.$sample.sh!);
	 system(qq!bash $sample_dir/quilts.$sample.sh!);
	}
   }		
