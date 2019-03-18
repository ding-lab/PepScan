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
my $type="";
my @dirs=("altsplice","altsplice_transcripts","other","variants",);

#my $root="/ifs/data/proteomics/tcga";
# (rjm)
my $root="/gscuser/rmashl/gsc/working/cptac/test_quilts";


my %stat=();
my $max_length=0;

if ($ARGV[0]=~/\w/) { $type=$ARGV[0];} else { $type="breast"; }

if ($error==0)
{
	opendir(DIR,"$root/samples/$type");
	my @samples = readdir DIR;
	closedir DIR;
	my @samples_ = ();
	my %versions_ = ();
	foreach my $sample (sort @samples)
	{
		if ($sample=~/\w/)
		{
			if (opendir(DIR,"$root/samples/$type/$sample/protein/fasta"))
			{
				my @versions = readdir DIR;
				closedir DIR;
				@samples_=(@samples_,$sample);
				print qq!$sample\n!;
				foreach my $version (@versions)
				{
					foreach my $dir (@dirs)
					{
						if (open (IN,"$root/samples/$type/$sample/protein/fasta/$version/$dir/proteome.fasta"))
						{
							if (open (LOG,">$root/samples/$type/$sample/protein/fasta/$version/$dir/proteome.log"))
							{
								$versions_{$version}=1;
								my $name="";
								my $description="";
								my $sequence="";
								my $line="";
								my %names=();
								my %sequences=();
								while ($line=<IN>)
								{
									chomp($line);
									if ($line=~/^>\s*(\S+)\s?(.*)$/)
									{
										my $name_=$1;
										my $description_=$2;
										$sequence=~s/\*$//;
										if ($name=~/\w/ and $sequence=~/\w/)
										{
											$names{$name}++;
											$sequences{$name}.="#$sequence#";
											my $length=length($sequence);
											if ($max_length<$length) { $max_length=$length; }
											$stat{"$sample#$version#$dir#$length"}++;
											if ($dir=~/^variants$/)
											{
												my $temp=$description;
												while($temp=~s/chr[0-9XYM]+\-[0-9\(\)]+\-([SG]):([A-Z\*])([0-9]+)([A-Z\-a-z\*]+),\s*$//)
												{
													$stat{"$sample#$version#$dir#$1"}++;
													$stat{"$sample#$version#$dir#$2#$4#$1"}++;
												}
											}
											my %characters=();
											while ($sequence=~s/^(.)//)
											{
												$characters{$1}++;
											}
											foreach my $char (keys %characters)
											{
												if ($char!~/[ACDEFGHIKLMNPQRSTUVWXY]/)
												{
													print LOG qq!Strange character: $characters{$char} $char\t$name\n!;
													$stat{"$sample#$version#$dir#strange"}++;
												}
											}
										}
										else
										{
											if ($name=~/\w/)
											{
												if ($sequence!~/\w/)
												{
													print LOG qq!Sequence missing\t$name\n!;
													$stat{"$sample#$version#$dir#0"}++;
													$stat{"$sample#$version#$dir#missing"}++;
												}
											}
										}
										$name=$name_;
										$description=$description_;
										$sequence="";
									}
									else
									{
										$sequence.="\U$line";
									}
								}
										$sequence=~s/\*$//;
										if ($name=~/\w/ and $sequence=~/\w/)
										{
											$names{$name}++;
											$sequences{$name}.="#$sequence#";
											my $length=length($sequence);
											if ($max_length<$length) { $max_length=$length; }
											$stat{"$sample#$version#$dir#$length"}++;
											if ($dir=~/^variants$/)
											{
												my $temp=$description;
												while($temp=~s/chr[0-9XYM]+\-[0-9\(\)]+\-([SG]):([A-Z\*])([0-9]+)([A-Z\-a-z\*]+),\s*$//)
												{
													$stat{"$sample#$version#$dir#$1"}++;
													$stat{"$sample#$version#$dir#$2#$4#$1"}++;
												}
											}
											my %characters=();
											while ($sequence=~s/^(.)//)
											{
												$characters{$1}++;
											}
											foreach my $char (keys %characters)
											{
												if ($char!~/[ACDEFGHIKLMNPQRSTUVWXY]/)
												{
													print LOG qq!Strange character: $characters{$char} $char\t$name\n!;
													$stat{"$sample#$version#$dir#strange"}++;
												}
											}
										}
										else
										{
											if ($name=~/\w/)
											{
												if ($sequence!~/\w/)
												{
													print LOG qq!Sequence missing\t$name\n!;
													$stat{"$sample#$version#$dir#0"}++;
													$stat{"$sample#$version#$dir#missing"}++;
												}
											}
										}
										
								foreach my $name (keys %names)
								{
									if ($names{$name}>1)
									{
										my $temp=$sequences{$name};
										my $same="Seqences identical";
										my $previous="";
										while($temp=~s/^#([^#]+)#//)
										{
											my $seq=$1;
											if ($previous=~/\w/)
											{
												if ($previous!~/^$seq$/) 
												{ 
													$same="Seqences different"; 
													$stat{"$sample#$version#$dir#duplicate_different"}++;
												}
												else
												{
													$stat{"$sample#$version#$dir#duplicate_same"}++;
												}
											}
											$previous=$seq;
										}
										print LOG qq!ID repeated $names{$name} times\t$name\t$same\n!;
									}
								}
								close(LOG);
							}
							close(IN);
						}
					}
				}
			}
		}
	}
	
	foreach my $version (keys %versions_)
	{
		if (open (OUT,">$root/samples/$type/stat/$version-fasta-length.stat"))
		{
			print OUT qq!protein-length!;
			foreach my $sample (@samples_)
			{
				foreach my $dir (@dirs)
				{
					print OUT qq!\t$sample-$dir!;
				}
			}
			print OUT qq!\n!;
			for(my $i=0;$i<=$max_length;$i++)
			{
				print OUT qq!$i!;
				foreach my $sample (@samples_)
				{
					foreach my $dir (@dirs)
					{
						if ($stat{"$sample#$version#$dir#$i"}!~/\w/) { $stat{"$sample#$version#$dir#$i"}=0; }
						print OUT qq!\t$stat{"$sample#$version#$dir#$i"}!;
					}
				}
				print OUT qq!\n!;
			}
			close(OUT);
		}	
		if (open (OUT,">$root/samples/$type/stat/$version-fasta-errors.stat"))
		{
			print OUT qq!error-type!;
			foreach my $sample (@samples_)
			{
				foreach my $dir (@dirs)
				{
					print OUT qq!\t$sample-$dir!;
				}
			}
			print OUT qq!\n!;
			foreach my $error_type ("missing","strange","duplicate_different","duplicate_same",)
			{
				print OUT qq!$error_type!;
				foreach my $sample (@samples_)
				{
					foreach my $dir (@dirs)
					{
						if ($stat{"$sample#$version#$dir#$error_type"}!~/\w/) { $stat{"$sample#$version#$dir#$error_type"}=0; }
						print OUT qq!\t$stat{"$sample$version#$dir#$error_type"}!;
					}
				}
				print OUT qq!\n!;
			}
			close(OUT);
		}
		
		
		if (open (OUT,">$root/samples/$type/stat/$version-fasta-germline-variants.stat"))
		{
			my $dir="variants";
			print OUT qq!error-type!;
			foreach my $sample (@samples_)
			{
				print OUT qq!\t$sample-$dir!;
			}
			print OUT qq!\n!;
			foreach my $aa1 ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*",)
			{
				foreach my $aa2 ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*",)
				{
					print OUT qq!$aa1->$aa2!;
					foreach my $sample (@samples_)
					{
						if ($stat{"$sample#$version#$dir#$aa1#$aa2#G"}!~/\w/) { $stat{"$sample#$version#$dir#$aa1#$aa2#G"}=0; }
						print OUT qq!\t$stat{"$sample#$version#$dir#$aa1#$aa2#G"}!;
					}
					print OUT qq!\n!;
				}
			}
			close(OUT);
		}
		if (open (OUT,">$root/samples/$type/stat/$version-fasta-somatic-variants.stat"))
		{
			my $dir="variants";
			print OUT qq!error-type!;
			foreach my $sample (@samples_)
			{
				print OUT qq!\t$sample-$dir!;
			}
			print OUT qq!\n!;
			foreach my $aa1 ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*",)
			{
				foreach my $aa2 ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","*",)
				{
					print OUT qq!$aa1->$aa2!;
					foreach my $sample (@samples_)
					{
						if ($stat{"$sample#$version#$dir#$aa1#$aa2#S"}!~/\w/) { $stat{"$sample#$version#$dir#$aa1#$aa2#S"}=0; }
						print OUT qq!\t$stat{"$sample#$version#$dir#$aa1#$aa2#S"}!;
					}
					print OUT qq!\n!;
				}
			}
			close(OUT);
		}
	}
	system(qq!chmod u+x $script_dir/plot_stat.sh!);
	system(qq!$script_dir/plot_stat.sh!);
}
