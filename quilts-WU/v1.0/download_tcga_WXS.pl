#!c:/perl/bin/perl.exe
#
use strict;

# (rjm)
use File::Path(qw(make_path));

my $error=0;
my $filename="";

if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }

#my $stratgies_to_use="#WXS#WGS#RNA-Seq#";
my $stratgies_to_use="#WXS#";
if ($error==0)
{
	if (open (IN,"$filename"))
	{
		if (open (LOG,">$filename.log"))
		{
			my $line="";
			while ($line=<IN>)
			{
				chomp($line);
				#print qq!$line\n!;
				if ($line=~/^(TCGA\-[^\-]+\-[^\-]+)/)
				{
					my $participant=$1;
					if (opendir(dir,"$participant"))
					{
						closedir(dir);
					}
					else
					{
						system(qq!cp -R ../template $participant!);
					}
					my $count=0;
					my $count_=0;
					my $file_count=0;
					my $file_name="";
					my %checksum=();
					my %data=();
					print qq!$participant\n!;
					system(qq!cgquery "legacy_sample_id=$participant*" > $participant/$participant.cgquery!);
					if (open (IN_,"$participant/$participant.cgquery"))
					{
						while ($line=<IN_>)
						{
							chomp($line);
							if ($line=~/^\s*Analysis\s+([0-9]+)\s*$/)
							{
								$count=$1;
								$count_++;
								if ($count==$count_)
								{
									#print LOG qq!Started $count\n!;
								}
								else
								{
									print LOG qq!Error: $count\!=$count_\n!;
								}
								%data=();
								$file_count=0;
								$file_name="";
								%checksum=();
							}
							if ($line=~/^\s*file\s+([0-9]+)\s*$/)
							{
								$file_count=$1;
								$file_name="";
							}
							if ($line=~/^\s*filename\s*:\s*(.+)\s*$/)
							{
								$file_name=$1;
							}
							if ($line=~/^\s*checksum\s*:\s*(.+)\s*$/)
							{
								my $checksum=$1;
								if ($file_count>0 and $file_name=~/\w/)
								{
									$checksum{$file_name}=$checksum;
								}
							}
							if ($line=~/^\s*(\w+)\s*:\s*(.+)\s*$/)
							{
								$data{$1}=$2;
							}
							if ($line!~/\w/)
							{
								if($count>0 and $data{'analysis_id'}=~/\w/)
								{
									if ($stratgies_to_use=~/#$data{'library_strategy'}#/i)
									{
										my $ok=0;
										my $target="rna";
										if ($data{'library_strategy'}=~/^W.S$/) 
										{
											$target="dna"; 
											if ($data{'sample_type'}=~/^10$/) { $target="dna-germline"; }
											if ($data{'sample_type'}=~/^01$/ or $data{'sample_type'}=~/^10$/) { $ok=1; }
										}
										else
										{
											if ($data{'sample_type'}=~/^01$/) { $ok=1; }
										}
										if ($ok==1)
										{
										    # mkdir("$participant/$target/bam/$data{'analysis_id'}");
										  # (rjm)
										    &make_path("$participant/$target/bam/$data{'analysis_id'}");
											if (open (OUT_,">$participant/$target/bam/$data{'analysis_id'}/tcga.md5"))
											{
												foreach my $file_name (keys %checksum)
												{
													my $file_name_=$file_name;
													$file_name_=~s/\.bam//;
													print OUT_ qq!$checksum{$file_name}  $participant/$target/bam/$data{'analysis_id'}/$file_name\n!;
												}
												close(OUT);
											}
											print LOG qq!*\t!;
											print qq!$participant $data{'sample_type'} $data{'analysis_id'} $data{'library_strategy'} $data{'refassem_short_name'} $data{'disease_abbr'}!;
											my $ok_=0;
											my $ok__=1;
											if (open (IN__,"$participant/$target/bam/$data{'analysis_id'}/md5sumcheck.log"))
											{
												my $line="";
												while ($line=<IN__>)
												{
													if ($line=~/\w/)
													{
														if ($line!~/\: OK/)
														{
															$ok__=0;
														}
														$ok_++;
													}
												}
												close(IN__);
											}
											if ($ok_==0 or $ok__==0)
											{
												system(qq!md5sum --check $participant/$target/bam/$data{'analysis_id'}/tcga.md5 > $participant/$target/bam/$data{'analysis_id'}/md5sumcheck.log!);
												$ok_=0;
												$ok__=1;
												if (open (IN__,"$participant/$target/bam/$data{'analysis_id'}/md5sumcheck.log"))
												{
													my $line="";
													while ($line=<IN__>)
													{
														if ($line=~/\w/)
														{
															if ($line!~/\: OK/)
															{
																$ok__=0;
															}
															$ok_++;
														}
													}
													close(IN__);
												}
												if ($ok_==0 or $ok__==0)
												{
													system(qq!gtdownload -d $data{'analysis_id'} -c /ifs/data/proteomics/tcga/cghub.key -p $participant/$target/bam!);
													system(qq!md5sum --check $participant/$target/bam/$data{'analysis_id'}/tcga.md5 > $participant/$target/bam/$data{'analysis_id'}/md5sumcheck.log!);
													print qq! Done\n!;
													print LOG qq!Done\t!;
												}
												else
												{
													print qq! Skipped\n!;
													print LOG qq!Skipped\t!;
												}
											}
											else
											{
												print qq! Skipped\n!;
												print LOG qq!Skipped\t!;
											}
										}
										else
										{
											print LOG qq!\t\t!;
										}
									}
									else
									{
										print LOG qq!\t\t!;
									}
									print LOG qq!$participant\t$data{'participant_id'}\t$data{'sample_type'}\t$data{'analysis_id'}\t$data{'library_strategy'}\t$data{'platform'}\t$data{'refassem_short_name'}\t$data{'disease_abbr'}\n!;
								}
								%data=();
							}
						}
					}
				}
			}
			close(LOG);
		}
		close(IN);
	}
}
