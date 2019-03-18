#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $dir="";

if ($ARGV[0]=~/\w/) { $dir=$ARGV[0];} else { $dir="."; }

if ($error==0)
{
	opendir(DIR,"$dir");
	my @files = readdir DIR;
	closedir DIR;
	foreach my $filename (@files)
	{
		if ($filename=~/\.fasta$/i)
		{
			if (open (IN,"$dir/$filename"))
			{
				open (OUT,">$dir/$filename.cleaned.fasta");
				if (open (LOG,">$dir/$filename.log"))
				{
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
							$sequence=~s/\*+$//;
							if ($name=~/\w/ and $sequence=~/\w/)
							{
								my $this_error=0;
								$names{$name}++;
								$sequences{$name}.="#$sequence#";
								my $sequence_=$sequence;
								my %characters=();
								while ($sequence_=~s/^(.)//)
								{
									$characters{$1}++;
								}
								foreach my $char (keys %characters)
								{
									if ($char!~/[ACDEFGHIKLMNPQRSTUVWXY]/)
									{
										print LOG qq!Strange character: $characters{$char} $char\t$name\n!;
										$error=1;
										$this_error=1;
									}
								}
								if($this_error==0)
								{
									print OUT qq!>$name!;
									if ($description=~/\w/) { print OUT qq! $description!; }
									print OUT qq!\n$sequence\n!;
								}
							}
							else
							{
								if ($name=~/\w/)
								{
									if ($sequence!~/\w/)
									{
										print LOG qq!Sequence missing\t$name\n!;
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
							$sequence=~s/\*+$//;
							if ($name=~/\w/ and $sequence=~/\w/)
							{
								my $this_error=0;
								$names{$name}++;
								$sequences{$name}.="#$sequence#";
								my $sequence_=$sequence;
								my %characters=();
								while ($sequence_=~s/^(.)//)
								{
									$characters{$1}++;
								}
								foreach my $char (keys %characters)
								{
									if ($char!~/[ACDEFGHIKLMNPQRSTUVWXY]/)
									{
										print LOG qq!Strange character: $characters{$char} $char\t$name\n!;
										$error=1;
										$this_error=1;
									}
								}
								if($this_error==0)
								{
									print OUT qq!>$name!;
									if ($description=~/\w/) { print OUT qq! $description!; }
									print OUT qq!\n$sequence\n!;
								}
							}
							else
							{
								if ($name=~/\w/)
								{
									if ($sequence!~/\w/)
									{
										print LOG qq!Sequence missing\t$name\n!;
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
									if ($previous!~/^$seq$/) { $same="Seqences different"; }
								}
								$previous=$seq;
							}
							print LOG qq!ID repeated $names{$name} times\t$name\t$same\n!;
						}
					}
					close(LOG);
				}
				close(IN);
				close(OUT);
			}
		}
	}
}
