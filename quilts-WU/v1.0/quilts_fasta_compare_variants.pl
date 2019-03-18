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
									print LOG qq!($names{$name}) $status $name $description   (!;
									if ($status=~/Different/)
									{
										my $string1=$sequence;
										my $string2=$sequences{$name};
										my $pre_sum=0;
										my $diff_old="";
										my $diff_new="";
										my $diff_start="";
										my $diff_end="";
										while ($string1=~/\w/)
										{
											my $pre = length( ( ( $string1 ^ $string2 ) =~ /^(\0*)/ )[ 0 ] );
											if ($pre>0 and length($diff_new)>0)
											{
												if ($diff_start!=$diff_end) { print LOG qq!$diff_start-$diff_end: !; }
												else { print LOG qq!$diff_start: !; }
												print LOG "$diff_old->";
												print LOG "$diff_new, ";
												$diff_old="";
												$diff_new="";
												$diff_start="";
												$diff_end="";
											}
											$pre_sum+=$pre;
											if (length($string1)!=$pre)
											{
												my $new=substr $string1,$pre,1;
												my $old=substr $string2,$pre,1;
												if ($diff_new!~/\w/) { $diff_start=$pre_sum; }
												$diff_end=$pre_sum;
												$diff_old.=$old;
												$diff_new.=$new;
												substr $string1,0,$pre+1,"";
												substr $string2,0,$pre+1,"";
												$pre_sum++
											}
											else
											{
												substr $string1,0,$pre,"";
												substr $string2,0,$pre,"";
											}
										}
										if (length($diff_new)>0)
										{
											if ($diff_start!=$diff_end) { print LOG qq!$diff_start-$diff_end: !; }
											else { print LOG qq!$diff_start: !; }
											print LOG "$diff_old->";
											print LOG "$diff_new, ";
											$diff_old="";
											$diff_new="";
											$diff_start="";
											$diff_end="";
										}
										if ($string2=~/\w/)
										{
											print LOG qq!$pre_sum: $string2-> , !;
										}
									}
									print LOG ")\n";
									#print LOG qq!$sequence\n!;
									#print LOG qq!$sequences{$name}\n!;
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
									print LOG qq!($names{$name}) $status $name $description   (!;
									if ($status=~/Different/)
									{
										my $string1=$sequence;
										my $string2=$sequences{$name};
										my $pre_sum=0;
										my $diff_old="";
										my $diff_new="";
										my $diff_start="";
										my $diff_end="";
										while ($string1=~/\w/)
										{
											my $pre = length( ( ( $string1 ^ $string2 ) =~ /^(\0*)/ )[ 0 ] );
											if ($pre>0 and length($diff_new)>0)
											{
												if ($diff_start!=$diff_end) { print LOG qq!$diff_start-$diff_end: !; }
												else { print LOG qq!$diff_start: !; }
												print LOG "$diff_old->";
												print LOG "$diff_new, ";
												$diff_old="";
												$diff_new="";
												$diff_start="";
												$diff_end="";
											}
											$pre_sum+=$pre;
											if (length($string1)!=$pre)
											{
												my $new=substr $string1,$pre,1;
												my $old=substr $string2,$pre,1;
												if ($diff_new!~/\w/) { $diff_start=$pre_sum; }
												$diff_end=$pre_sum;
												$diff_old.=$old;
												$diff_new.=$new;
												substr $string1,0,$pre+1,"";
												substr $string2,0,$pre+1,"";
												$pre_sum++
											}
											else
											{
												substr $string1,0,$pre,"";
												substr $string2,0,$pre,"";
											}
										}
										if (length($diff_new)>0)
										{
											if ($diff_start!=$diff_end) { print LOG qq!$diff_start-$diff_end: !; }
											else { print LOG qq!$diff_start: !; }
											print LOG "$diff_old->";
											print LOG "$diff_new, ";
											$diff_old="";
											$diff_new="";
											$diff_start="";
											$diff_end="";
										}
										if ($string2=~/\w/)
										{
											print LOG qq!$pre_sum: $string2-> , !;
										}
									}
									print LOG ")\n";
									#print LOG qq!$sequence\n!;
									#print LOG qq!$sequences{$name}\n!;
								}
					}
				}
			}
		}
	}
}
