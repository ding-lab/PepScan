#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $dir="";
my $db="";


#my $root="/ifs/data/proteomics/tcga/databases";
# (rjm)
my $root="/gscuser/rmashl/gsc/working/cptac/test_quilts/databases";

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $db=$ARGV[1];} else { $db="ensembl_human_37.70"; }

if ($error==0)
{
	if (open (IN,"$root/$db/proteome.fasta"))
	{
		my $name="";
		my $description="";
		my $sequence="";
		my $line="";
		my %names=();
		my %sequences=();
		my %descriptions=();
		while ($line=<IN>)
		{
			chomp($line);
			$line=~s/[\n\r]+$//;
			if ($line=~/^>\s*(\S+)\s?(.*)$/)
			{
				my $name_=$1;
				my $description_=$2;
				if ($name=~/\w/ and $sequence=~/\w/)
				{
					$names{$name}++;
					$sequences{$name}=$sequence;
					$descriptions{$name}=$description;
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
				if ($name=~/\w/ and $sequence=~/\w/)
				{
					$names{$name}++;
					$sequences{$name}=$sequence;
					$descriptions{$name}=$description;
				}
		opendir(DIR,"$dir");
		my @files = readdir DIR;
		closedir DIR;
		foreach my $filename (sort @files)
		{
			if ($filename=~/\.fasta$/)
			{
				print "$dir/$filename\n";
				if (open (IN,"$dir/$filename"))
				{
					if (open (LOG,">$dir/$filename.compare.log"))
					{
						my $name="";
						my $description="";
						my $sequence="";
						my $line="";
						my $status="";
						while ($line=<IN>)
						{
							chomp($line);
							if ($line=~/^>\s*(\S+)\s?(.*)$/)
							{
								my $name_=$1;
								my $description_=$2;
								$name_=~s/^([^\-]+\-[^\-]+).*/$1/;
								if ($name=~/\w/ and $sequence=~/\w/)
								{
									if ($sequences{$name}=~/^$sequence$/) { $status="Equal"; } else { $status="Different"; }
									print LOG qq!\n>($names{$name}) $status $name $description\n!;
									print LOG qq!$sequence\n!;
									print LOG qq!$sequences{$name}\n!;
									if ($status=~/Different/)
									{
										my $string1=$sequence;
										my $string2=$sequences{$name};
										my $pre = length( ( ( $string1 ^ $string2 ) =~ /^(\0*)/ )[ 0 ] );
										my $post = length( ( ( reverse( $string1 ) ^ reverse( $string2 ) ) =~ /^(\0*)/ )[ 0 ] );    
										my $string1_ = substr( $string1, $pre, -$post ); 
										my $string2_ = substr( $string2, $pre, -$post );
										my $end=length($string2_)+$pre;
										print LOG "$pre-$end: $string2_->$string1_\n";
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
								if ($name=~/\w/ and $sequence=~/\w/)
								{
									if ($sequences{$name}=~/^$sequence$/) { $status="Equal"; } else { $status="Different"; }
									print LOG qq!\n>($names{$name}) $status $name $description\n!;
									print LOG qq!$sequence\n!;
									print LOG qq!$sequences{$name}\n!;
									if ($status=~/Different/)
									{
										my $string1=$sequence;
										my $string2=$sequences{$name};
										my $pre = length( ( ( $string1 ^ $string2 ) =~ /^(\0*)/ )[ 0 ] );
										my $post = length( ( ( reverse( $string1 ) ^ reverse( $string2 ) ) =~ /^(\0*)/ )[ 0 ] );    
										my $string1_ = substr( $string1, $pre, -$post ); 
										my $string2_ = substr( $string2, $pre, -$post );
										my $end=length($string2_)+$pre;
										print LOG "$pre-$end: $string2_->$string1_\n";
									}
								}
					}
				}
			}
		}
	}
}
