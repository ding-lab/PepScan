# creat the db as the input of msgf+ 
# Song Cao, Ding Lab
# scao@genome.wustl.edu 

#!/usr/bin/perl
use strict;
use warnings;

my $usage = '
The input is the directory of mzml file and will output the db from Quilts. 
';

die $usage unless scalar @ARGV == 1;
my ($input_dir_mzml) = @ARGV;
my $mzml_name=(split(/\//,$input_dir_mzml))[-1]; 
my @sample_name=split(/\_/,$mzml_name);  
my $sample1="TCGA-".$sample_name[1]; 
my $sample2="TCGA-".$sample_name[2]; 
my $sample3="TCGA-".$sample_name[3]; 
my @working_name= split(/\//,$input_dir_mzml);
my $quilts_dir=$working_name[0]; 
my $indel_dir=$working_name[0];


my $out_file=$input_dir_mzml."/".$mzml_name.".fasta"; 
open(OUT,">$out_file"); 

for(my $i=1;$i<=scalar @working_name-3;$i++)
 {
   $quilts_dir=join("/",$quilts_dir,$working_name[$i]); 	
 }

for(my $i=1;$i<=scalar @working_name-3;$i++)
 {
   $indel_dir=join("/",$indel_dir,$working_name[$i]);
 }

$quilts_dir=join("/",$quilts_dir,"quilts");
$indel_dir=join("/",$indel_dir,"indel");
 
print $quilts_dir, "\n";
print $indel_dir,"\n";

my $sample1_variant=$quilts_dir."/".$sample1."/protein/fasta/QUILTSv1.0-ensembl_human_37.70-var-vcf-rna-junction/variants/proteome.fasta"; 
my $sample2_variant=$quilts_dir."/".$sample2."/protein/fasta/QUILTSv1.0-ensembl_human_37.70-var-vcf-rna-junction/variants/proteome.fasta";
my $sample3_variant=$quilts_dir."/".$sample3."/protein/fasta/QUILTSv1.0-ensembl_human_37.70-var-vcf-rna-junction/variants/proteome.fasta";

my $sample1_indel=$indel_dir."/".$sample1."/result/proteome-indel-mod.fasta";
my $sample2_indel=$indel_dir."/".$sample2."/result/proteome-indel-mod.fasta";
my $sample3_indel=$indel_dir."/".$sample3."/result/proteome-indel-mod.fasta";

my $db_ensembl="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/Quilts/databases/ensembl_human_37.70/proteome.fasta"; 
my $db_crap="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/Quilts/databases/crap/proteome.fasta"; 

open(IN_crap,"<$db_crap"); 
open(IN_ensembl,"<$db_ensembl"); 
open(IN_S1,"<$sample1_variant"); 
open(IN_S2,"<$sample2_variant"); 
open(IN_S3,"<$sample3_variant"); 

#print $sample1_indel,"\n"; 

open(IN_S4,"<$sample1_indel");
open(IN_S5,"<$sample2_indel");
open(IN_S6,"<$sample3_indel");

while(<IN_ensembl>)
{
  my $line=$_; 
  print OUT $line; 
}

while(<IN_crap>)
{
  my $line=$_; 
  print OUT $line; 
}

while(<IN_S1>)
{
  my $line=$_; 
  if($line=~/^\>/) { $line=~s/variant/variant-$sample1/g; }
  print OUT $line; 
}

while(<IN_S2>)
{
  my $line=$_;
  if($line=~/^\>/) { $line=~s/variant/variant-$sample2/g; }
  print OUT $line;
}

while(<IN_S3>)
{
  my $line=$_;
  if($line=~/^\>/) { $line=~s/variant/variant-$sample3/g; }
  print OUT $line;
}

while(<IN_S4>)
{
  my $line=$_;
  if($line=~/^\>/) { $line=~s/indel/indel-$sample1/g; }
  print OUT $line;
}

while(<IN_S5>)
{
  my $line=$_;
  if($line=~/^\>/) { $line=~s/indel/indel-$sample2/g; }
  print OUT $line;
}

while(<IN_S6>)
{
  my $line=$_;
  if($line=~/^\>/) { $line=~s/indel/indel-$sample3/g; }
  print OUT $line;
}

close IN_ensembl; 
close IN_crap; 
close IN_S1; 
close IN_S2; 
close IN_S3;
close IN_S4; 
close IN_S5; 
close IN_S6;
close OUT;  
