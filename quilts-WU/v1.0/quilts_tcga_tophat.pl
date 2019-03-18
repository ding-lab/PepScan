#!/usr/bin/perl
use strict;
use File::Spec;
use File::Path(qw(make_path));
my $script_path = File::Spec->rel2abs(__FILE__);
my $script_dir=$script_path;
my $version=$script_path;
$script_dir=~s/^(.*)\/[^\/]+$/$1/;
$version=~s/^.*\/([^\/]+)\/[^\/]+$/$1/;

my $usage = '
This script will run quilts for a given sample directory
';
die $usage unless scalar @ARGV >= 1;

print "version=$version\n";
#my $root="/ifs/data/proteomics/tcga";
my $root="/gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/Quilts";
my $scripts=$script_dir;

my $error=0;
my $sample="";
my $db="";
my $genome="";
my $vcf_source="";
my $junction_source="";

if ($ARGV[0]=~/\w/) { $sample=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $db=$ARGV[1];} else { $db="refseq_human_20130727_proteome"; }
if ($ARGV[2]=~/\w/) { $genome=$ARGV[2];} else { $genome="genome_human"; }
if ($ARGV[3]=~/\w/) { $vcf_source=$ARGV[3];} else { $vcf_source="vcf"; }

if ($ARGV[4]=~/\w/) { $junction_source=$ARGV[4];} else { $junction_source="junction"; }

my $junction_source_=$junction_source;
$junction_source_=~s/^tophat\/junction_//;
$junction_source_=~s/_gMx$//;

print $junction_source,"\n"; 

my $dir_version="QUILTS$version-$db-var-$vcf_source-rna-$junction_source_";
my $vcf_dir="$sample/$vcf_source";
my $junctions_dir="$sample/$junction_source";
my $db_dir="$root/databases/$db";
my $genome_dir="$root/databases/$genome";
my $result_dir="$sample/protein";
my $pgx_index="$root/Song_adapted/scripts/pgx-WU/v0.7/pgx_index.py";

&make_path("$result_dir/fasta/$dir_version");
&make_path("$result_dir/fasta/$dir_version/log");

open(LOG,">$result_dir/fasta/$dir_version/version.txt");
print LOG qq!version=$version\n!;
print LOG qq!root=$root\n!;
print LOG qq!scripts=$scripts\n!;

my $vcf_count=0;
my $vcf_gl="";
my $vcf="";
my $i=1; 
my %somatic = ();
my %germline = ();
my %header=(); 
if (opendir(dir,"$vcf_dir"))
{
	my @allfiles=readdir dir;
	closedir dir;
	foreach my $filename (@allfiles)
	{
		if ($filename=~/\.vcf$/ and $filename!~/\.germline\.vcf/ and $filename!~/\-somatic_only\.vcf/)
		{
			print "$filename\n"; 
			$vcf=$filename;
			$vcf_count++; 
		}
	}
	if ($vcf_count==1)
	{
		if (open (IN, "<$vcf_dir/$vcf"))
		{
			while (my $line=<IN>)
			{
				chomp($line);
				if ($line=~/^#/)
				{
					$header{$i}="$line\n";
					$i++; 
				}
				elsif ($line=~/^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)/)
				{
					my $GS=$8;
					if ($GS=~/SOMATIC/)
					{
						$somatic{$line}=0; 
					}
					else
					{
						$germline{$line}=0; 
					}
				}
			}
			close IN;
			if(open (OUTG, ">$vcf_dir/$vcf.germline.vcf") or die "can't open germline output file!")
			{
				foreach my $key (sort  {$a <=> $b} keys %header)
				{
					print OUTG "$header{$key}"; 
				}
				foreach my $key (keys %germline)
				{
					print OUTG "$key\n"; 
				}
				close OUTG; 
			}
		}
		$vcf_gl="$vcf.germline.vcf";
	}
}

print "vcf=",$vcf_count,"\n";

if ($vcf_count==0) { print qq!Warning: No tumor VCF file found\n!; $error=0; }
if ($vcf_count>1) { print qq!Error: More than one tumor VCF file found\n!; $error=1; }

if ($error==0)
{
	print "run merge_junctions_bed_dir.pl","\n"; 
	system(qq!perl $scripts/merge_junctions_bed_dir.pl $junctions_dir $result_dir/fasta/$dir_version/log/!);
	print "run filter_known_transcripts.pl","\n";
	system(qq!perl $scripts/filter_known_transcripts.pl $result_dir/fasta/$dir_version/log/merged-junctions.bed $db_dir/transcriptome.bed!);
	print "run filter_A.pl","\n"; 

	system(qq!perl $scripts/filter_A.pl $result_dir/fasta/$dir_version/log/merged-junctions.filter.bed $db_dir/proteome.bed!);
	#<STDIN>;
	print "cp proteme.bed","\n"; 

	system(qq!cp $db_dir/proteome.bed $result_dir/fasta/$dir_version/log!);

# (rjm)
	if ($vcf   =~/\w/) {print qq!Have somatic  vcf $vcf\n!}
	if ($vcf_gl=~/\w/) {print qq!Have germline vcf $vcf_gl\n!}

	print qq!Running protein_seq_from_genome_using_bed.pl...\n!;

	if ($vcf=~/\w/ and $vcf_gl=~/\w/)
	{
		system(qq!perl $scripts/protein_seq_from_genome_using_bed.pl $result_dir/fasta/$dir_version/log/proteome.bed $genome_dir $vcf_dir/$vcf_gl $vcf_dir/$vcf!);
	}
	else
	{
		system(qq!perl $scripts/protein_seq_from_genome_using_bed.pl $result_dir/fasta/$dir_version/log/proteome.bed $genome_dir!);
	}
       
        print qq!Finishing protein_seq_from_genome_using_bed.pl...\n!; 

	system(qq!mkdir -p $result_dir/fasta/$dir_version/variants!);

	system(qq!cp $result_dir/fasta/$dir_version/log/proteome.bed-mod.fasta $result_dir/fasta/$dir_version/variants/proteome.fasta!);

# (rjm) -- added conditional provisionally
	#if ( -e "$result_dir/fasta/$dir_version/log/proteome.bed-mod.bed" ) {
	system(qq!cp $result_dir/fasta/$dir_version/log/proteome.bed-mod.bed $result_dir/fasta/$dir_version/variants/proteome.bed!);
	#}
	
	print qq!Running quilts_fasta_compare_variants.pl...\n!; 

	system(qq!perl $scripts/quilts_fasta_compare_variants.pl $result_dir/fasta/$dir_version/variants!);
	print qq!Running check_fasta.pl...\n!;

	system(qq!perl $scripts/check_fasta.pl $result_dir/fasta/$dir_version/variants!);
	system(qq!python $pgx_index $result_dir/fasta/$dir_version/variants!);

	#system(qq!rm -f $result_dir/fasta/$dir_version/log/proteome.bed!);
 	print qq!Running protein_seq_from_genome_using_bed.pl\n!; 

	system(qq!perl $scripts/protein_seq_from_genome_using_bed.pl $result_dir/fasta/$dir_version/log/merged-junctions.filter.bed.A.bed $genome_dir $vcf_dir/$vcf_gl $vcf_dir/$vcf!);
	system(qq!mkdir -p $result_dir/fasta/$dir_version/altsplice_transcripts!);
	system(qq!cp $result_dir/fasta/$dir_version/log/merged-junctions.filter.bed.A.fasta $result_dir/fasta/$dir_version/altsplice_transcripts/proteome.fasta!);
	system(qq!cp $result_dir/fasta/$dir_version/log/merged-junctions.filter.bed.A.bed $result_dir/fasta/$dir_version/altsplice_transcripts/proteome.bed!);
	system(qq!perl $scripts/quilts_fasta_compare_alternative.pl $result_dir/fasta/$dir_version/altsplice_transcripts!);
	system(qq!perl $scripts/check_fasta.pl $result_dir/fasta/$dir_version/altsplice_transcripts!);
	system(qq!python $pgx_index $result_dir/fasta/$dir_version/altsplice_transcripts!);
	system(qq!mkdir -p $result_dir/fasta/$dir_version/altsplice!);
# (rjm) -- added conditional provisionally
	#if ( -e "$result_dir/fasta/$dir_version/log/merged-junctions.filter.bed.A.bed.alt.pep.fasta") {
	system(qq!cp $result_dir/fasta/$dir_version/log/merged-junctions.filter.bed.A.bed.alt.pep.fasta $result_dir/fasta/$dir_version/altsplice/proteome.fasta!);
	#}

	system(qq!perl $scripts/check_fasta.pl $result_dir/fasta/$dir_version/altsplice!);
	system(qq!python $pgx_index $result_dir/fasta/$dir_version/altsplice!);

	system(qq!perl $scripts/protein_seq_from_genome_using_junction.pl $result_dir/fasta/$dir_version/log/merged-junctions.filter.bed.notA.bed $genome_dir $vcf_dir/$vcf_gl $vcf_dir/$vcf!);
	system(qq!mkdir -p $result_dir/fasta/$dir_version/other!);
	system(qq!cp $result_dir/fasta/$dir_version/log/merged-junctions.filter.bed.notA.bed.fasta $result_dir/fasta/$dir_version/other/proteome.fasta!);
	system(qq!perl $scripts/check_fasta.pl $result_dir/fasta/$dir_version/other!);
	system(qq!python $pgx_index $result_dir/fasta/$dir_version/other!);
			
	if (open(OUT,">$result_dir/fasta/$dir_version/taxonomy.xml"))
	{
		print OUT qq!<?xml version="1.0"?>
<bioml label="x\! taxon-to-file matching list">
	<taxon label="sample">
	
		<file format="peptide" URL="$result_dir/fasta/$dir_version/variants/proteome.fasta" />
		<file format="peptide" URL="$result_dir/fasta/$dir_version/altsplice/proteome.fasta" />
		<file format="peptide" URL="$result_dir/fasta/$dir_version/other/proteome.fasta" />
		
		<file format="peptide" URL="$root/databases/$db/proteome.fasta" />
		<file format="mod" URL="$root/databases/$db/mod.xml" />
		<file format="saps" URL="$root/databases/$db/saps.xml" />
	
		<file format="peptide" URL="$root/databases/crap/proteome.fasta" />
		<file format="mod" URL="$root/databases/crap/crap_mod.xml" />
		<file format="saps" URL="$root/databases/crap/crap_saps.xml" />
		
	</taxon>
<bioml>
	!;
		close(OUT);
	}

	close(LOG);
}
