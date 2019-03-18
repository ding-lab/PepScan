#!c:/perl/bin/perl.exe
#
use strict;
use File::Spec;
my $script_path = File::Spec->rel2abs(__FILE__);
my $script_dir=$script_path;
my $version=$script_path;
$script_dir=~s/^(.*)\/[^\/]+$/$1/;
$version=~s/^.*\/([^\/]+)\/[^\/]+$/$1/;
print "version=$version\n";

my $error=0;
my $db="";
my $genome="";
my $root="";

if ($ARGV[0]=~/\w/) { $db=$ARGV[0];} else { $db="ensembl_human_37.70"; }
if ($ARGV[1]=~/\w/) { $genome=$ARGV[1];} else { $genome="genome_human"; }

$root="/ifs/data/proteomics/tcga";
if ($db!~/\//) { $db="$root/databases/$db"; }
if ($genome!~/\//) { $genome="$root/databases/$genome"; }
my $scripts="$script_dir";
#my $scripts_pgx="$script_dir/../../pgx/v0.7";
# (rjm)
my $scripts_pgx="$script_dir/../../pgx-WU/v0.7";
print qq!$db\n$genome\n$scripts\n!;

my $db_original=$db;
$db_original=~s/(ensembl_[^\_]+)_(.*)$/$1\_original\_$2/;
if (opendir(dir,"$db"))
{
	my $original_fasta="";
	my $original_fasta_count=0;
	my $gtf="";
	my $gtf_count=0;
	my @allfiles=readdir dir;
	closedir dir;
	foreach my $filename (@allfiles)
	{
		print qq!$filename\n!;
		if ($filename=~/\.pep\.all\.fa\.gz$/)
		{
			my $filename_=$filename;
			$filename_=~s/\.gz$//;
			$original_fasta=$filename_;
			$original_fasta.="sta";
			$original_fasta_count++;
			system(qq!gunzip $db/$filename!);
			system(qq!mv $db/$filename_ $db/$original_fasta!);
		}
		if ($filename=~/\.pep\.all\.fa$/)
		{
			$original_fasta=$filename;
			$original_fasta.="sta";
			$original_fasta_count++;
			system(qq!mv $db/$filename $db/$original_fasta!);
		}
		if ($filename=~/\.pep\.all\.fasta$/)
		{
			$original_fasta=$filename;
			$original_fasta_count++;
		}
		if ($filename=~/\.gtf\.gz$/)
		{
			$gtf=$filename;
			$gtf=~s/\.gz$//;
			$gtf_count++;
			system(qq!gunzip $db/$filename!);
		}
		if ($filename=~/\.gtf$/)
		{
			$gtf=$filename;
			$gtf_count++;
		}
	}
	print qq!original_fasta: $original_fasta ($original_fasta_count)\n!;
	print qq!gtf: $gtf ($gtf_count)\n!;
	if ($original_fasta_count!=1) { print qq!Error: more than one fasta file\n!; $error=1; }
	if ($gtf_count!=1) { print qq!Error: more than one gtf file\n!; $error=1; }

	if ($error==0)
	{
		if (open(OUT,">$db/prepare_ensembl.sh"))
		{
			print OUT qq!
cd $db
dos2unix *.gtf
dos2unix *.fasta
perl $scripts/get_ensembl_description.pl proteome
perl $scripts/check_and_clean_fasta.pl .
perl $scripts/ensembl_gtf_to_bed.pl $gtf
mv $gtf.prot.bed proteome-first.bed
mv $gtf.rna.bed transcriptome.bed
perl $scripts/protein_seq_from_genome_using_bed.pl proteome-first.bed $genome
rm -f proteome-first.bed-mod.fasta
rm -f proteome-first.bed-mod.log
rm -f proteome-first.bed-mod.stat
rm -f proteome-first.bed-frame-shift.fasta
perl $scripts/compare_seq_fasta.pl proteome-first.fasta $original_fasta.cleaned.fasta
mv proteome-first-corrected.bed proteome.bed
perl $scripts/protein_seq_from_genome_using_bed.pl proteome.bed $genome
rm -f proteome.bed-mod.fasta
rm -f proteome.bed-mod.log
rm -f proteome.bed-mod.stat
rm -f proteome.bed-frame-shift.fasta
perl $scripts/compare_seq_fasta.pl proteome.fasta $original_fasta.cleaned.fasta
python $scripts_pgx/pgx_index.py $db

mkdir $db_original
cd $db_original
cp $db/$original_fasta.cleaned.fasta proteome.fasta
dos2unix proteome.fasta
python $scripts_pgx/pgx_index.py $db_original
cd $db
			!;
			close(OUT);
			system(qq!qsub $db/prepare_ensembl.sh!);
		}
	}
}
