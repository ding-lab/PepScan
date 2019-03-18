#!c:/perl/bin/perl.exe

# Script a file with protein descriptions;
# Input is a folder containing .dat files and

use strict;
use Data::Dumper;

my $error=0;
my $dir="";
my $line="";

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $error=1; }

if ($error==0)
{
	opendir ( DIR, "$dir" ) || die "Error in opening dir $dir\n";
	my @allDat=readdir DIR;
	my %descriptions=();
	my %description_gene=();
	my %gene=();
	foreach my $file (sort @allDat)
	{
		if ($file =~/\.dat$/i)
		{
			print qq!$file\n!;
			#goes through each .dat file and makes a hash "gene_protein" (key is gene and value all proteins in that gene)
			open (IN, qq!$dir/$file!) || die "Could not open $dir/$file\n";
			my $status="";
			my $gene_name="";
			my $gene="";
			my $description_started=0;
			my $protein_gene="";
			my $protein="";
			while($line = <IN>)
			{	
				chomp($line);
				if ($line=~/^FT   (\S+)/)
				{
					$status = $1;
					$gene_name="";
					$gene="";
					$description_started=0;
					$protein_gene="";
					$protein="";
				}
				if($status=~/^CDS$/)
				{
					if ($line=~/FT\s+\/gene=\"([^\"]+)\"/)
					{
						$protein_gene=$1;
						$protein="";
					}
					if ($line=~/FT\s+\/protein_id=\"([^\"]+)\"/)
					{
						$protein=$1;
						$descriptions{"$protein-$gene{$protein_gene}"}=$description_gene{$protein_gene};
					}
				}
				if($status=~/^gene$/)
				{
					if ($description_started==1)
					{
						if ($line=~/FT\s+([^\"]+)\"?/)
						{
							$description_gene{$gene_name}.=" $1";
						}
						if ($line=~/FT\s+([^\"]+)\"/)
						{
							$description_started=0;
						}
					}
					if ($line=~/FT\s+\/gene=\"?([^\"]+)\"?/)
					{
						$gene_name=$1;
						$gene{$gene_name}="";
						$description_gene{$gene_name}="";
					}
					if ($line=~/FT\s+\/locus_tag=\"([^\"]+)\"/)
					{
						$gene{$gene_name}=$1;
					}
					if ($line=~/FT\s+\/note=\"([^\"]+)\"?/)
					{
						$description_gene{$gene_name}=$1;
						$description_started=1;
					}
				}
			}
			close (IN);
		}
	}

	open(OUT, qq!>$dir-descriptions.txt!);
	foreach my $key (sort keys %descriptions)
	{
		print OUT "$key\t$descriptions{$key}\n";
	}
	closedir DIR;

}