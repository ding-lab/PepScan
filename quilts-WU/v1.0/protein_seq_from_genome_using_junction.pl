#!c:/perl/bin/perl.exe
#
use strict;

my $error=0;
my $filename="";
my $filename_vcf_germline="";
my $filename_vcf_somatic="";
my $dir=".";

if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }
if ($ARGV[1]=~/\w/) { $dir=$ARGV[1];} else { $dir="."; }
if ($ARGV[2]=~/\w/) { $filename_vcf_germline=$ARGV[2];} else { $filename_vcf_germline=""; }
if ($ARGV[3]=~/\w/) { $filename_vcf_somatic=$ARGV[3];} else { $filename_vcf_somatic=""; }

my %mapping = (	"TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L",
				"CTT"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L",
				"ATT"=>"I","ATC"=>"I","ATA"=>"I","ATG"=>"M",
				"GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
				
				"TCT"=>"S","TCC"=>"S","TCA"=>"S","TCG"=>"S",
				"CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P",
				"ACT"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
				"GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
				
				"TAT"=>"Y","TAC"=>"Y","TAA"=>"*","TAG"=>"*",
				"CAT"=>"H","CAC"=>"H","CAA"=>"Q","CAG"=>"Q",
				"AAT"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K",
				"GAT"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
				
				"TGT"=>"C","TGC"=>"C","TGA"=>"*","TGG"=>"W",
				"CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
				"AGT"=>"S","AGC"=>"S","AGA"=>"R","AGG"=>"R",
				"GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G");
				
if ($error==0)
{
	open (OUT,">$filename.fasta");
	open (OUT_MOD,">$filename-mod.fasta");
	open (BED,">$filename-mod.bed");
	open (LOG,">$filename-mod.log");
	open (STAT,">$filename-mod.stat");
	my $proteins_count=0;
	my $proteins_modified=0;
	my $count_stop_removed=0;
	my $count_stop_introduced=0;
	my $protein_modifications_count=0;
	my @protein_modifications_distr=();
	my $count_variant_in_exon=0;
	my $count_variant_in_exon_old_error=0;
	my $count_variant_in_exon_syn=0;
	my $count_variant_in_exon_nonsyn=0;
	my $both_vcf=0;
	my $line="";
	my %chr=();
	my %bed=();
	my %seq=();
	if (open (IN,"$filename"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
			{
				my $chr=$1;
				my $name=$4;
				$bed{$name}=$line;
				$chr{$chr}.="#$name#";
			} else { print LOG qq!Error parsing: $line!; }
		}
		close(IN);
	}
	my %vcf_old=();
	my %vcf_new=();
	my %vcf_type=();
	if (open (IN,"$filename_vcf_germline"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
			{
				my $chr="chr$1";
				my $pos=$2;
				my $id=$3;
				my $old=$4;
				my $new=$5;
				my $qaul=$6;
				$pos--;
				$new=~s/\,.*$//; #####
				$vcf_old{"$chr#$pos"}=$old;
				$vcf_new{"$chr#$pos"}=$new;
				$vcf_type{"$chr#$pos"}="G";
				#print qq!$chr#$pos#$new\n!;
			} else { print LOG qq!Error parsing: $line!; }
		}
		close(IN);
	}
	if (open (IN,"$filename_vcf_somatic"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
			{
				my $chr="chr$1";
				my $pos=$2;
				my $id=$3;
				my $old=$4;
				my $new=$5;
				my $qaul=$6;
				$pos--;
				$new=~s/\,.*$//; #####
				if ($vcf_old{"$chr#$pos"}=~/\w/ and $vcf_new{"$chr#$pos"}!~/^$new$/)
				{
					print LOG qq!$chr $pos: germline:$vcf_old{"$chr#$pos"}->$vcf_new{"$chr#$pos"} somatic:$old->$new\n!;
					$both_vcf++;
				}
				$vcf_old{"$chr#$pos"}=$old;
				$vcf_new{"$chr#$pos"}=$new;
				$vcf_type{"$chr#$pos"}="S";
				#print qq!$chr#$pos#$new\n!;
			} else { print LOG qq!Error parsing: $line!; }
		}
		close(IN);
	}
	foreach my $chr (sort keys %chr)
	{
		print qq!$chr\n!;
		if (open (IN,"$dir/$chr.fa"))
		{
			print qq!opened $chr\n!;
			my $sequence="";
			$line=<IN>;
			chomp($line);
			if ($line=~/^>$chr\s*$/)
			{
				while ($line=<IN>)
				{
					chomp($line);
					if ($line=~/^>/)
					{
						print qq!Error: > not expected: $line\n!;
					}
					else
					{
						$line=~s/\s+//g;
						if ($line!~/^[atcgATCGnN]+$/)
						{
							print qq!Error: unexpected character: $line\n!;
						}
						else
						{
							$sequence .= "\U$line";
						}
					}
				}
				my $temp=$chr{$chr};
				while ($temp=~s/^#([^#]+)#//)
				{
					my $name=$1;
					my $modified=0;
					print LOG qq!\n$name: $bed{$name}\n!;
					if ($bed{$name}=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
					{
						my $chr=$1;
						my $start=$2;
						my $end=$3;
						my $num_observ=$5;
						my $strand=$6;
						my $num=$10;
						my $length_exons = $11;
						my $begin_exons = $12;
						my @length_exon = split(',', $length_exons,2);
						my @begin_exon = split(',', $begin_exons,2);
						my $begin_intron = $start + $length_exon[0] + 1;
						my $end_intron = $start + $begin_exon[-1] - 1;
						my $length_intron = $end_intron - $begin_intron;
						$start=$begin_intron-300-1;
						$end=$end_intron+300+1;
						my $start2=300+$length_intron+2;
						my $segment_lengths="300,300,";
						my $segment_starts="0,$start2,";
						print LOG qq!$chr\t$start\t$end\t$name\t1000\t$strand\t$start\t$end\t200,0,0\t2\t$segment_lengths\t$segment_starts\n!;
						my $seq="";
						my $seq_original="";
						my $segment_starts_=$segment_starts;
						my $segment_lengths_=$segment_lengths;
						my %variants=();
						my %variants_=();
						my $segment_start_first=0;
						if ($segment_starts_=~/^([0-9\-]+)\,/) { $segment_start_first=$1; }
						while ($segment_starts_=~s/^([0-9\-]+)\,//)
						{
							my $segment_start=$1;
							if ($segment_lengths_=~s/^([0-9]+)\,//)
							{
								my $segment_length=$1;
								my $seq_=substr $sequence,$start+$segment_start,$segment_length;
								$seq_original.=$seq_; 
								print LOG qq!$start+$segment_start,$segment_length: $seq_\n!;
								for(my $i=0,my $j=$start+$segment_start;$i<$segment_length;$i++,$j++)
								{
									my $n=substr $seq_,$i,1;
									#print qq!$chr#$j: $n $vcf_new{"$chr#$j"}\n!;
									my $temp="";
									if ($vcf_new{"$chr#$j"}=~/\w/ and $vcf_new{"$chr#$j"}!~/^$vcf_old{"$chr#$j"}$/)
									{
										if ($n!~/^$vcf_old{"$chr#$j"}$/) 
										{ 
											$temp=qq! (Warning: $n\!\=$vcf_old{"$chr#$j"}) - Ignored!; 
											my $i_=$i+$segment_start_first+length($seq);
											print LOG qq!Variant $chr $j: $n,$vcf_old{"$chr#$j"}$temp->$vcf_new{"$chr#$j"} $i_\n!;
											if ($temp=~/\w/) { $count_variant_in_exon_old_error++; }
										}
										else
										{
											$modified++;
											substr $seq_,$i,1,$vcf_new{"$chr#$j"};
											my $i_=$i+$segment_start_first+length($seq);
											print LOG qq!Variant $chr $j: $n,$vcf_old{"$chr#$j"}$temp->$vcf_new{"$chr#$j"} $i_\n!;
											if ($temp=~/\w/) { $count_variant_in_exon_old_error++; }
											$count_variant_in_exon++;
											$variants{$i_}=qq!$vcf_old{"$chr#$j"}#$vcf_new{"$chr#$j"}!;
										}
									}
								}
								$seq.=$seq_; 
							} else { print qq!Error parsing $bed{$name}\n!; }
						}
						my $name_=$name; $name_=~s/\-[^\-]+$//;
						my $gene=$name; $gene=~s/^[^\-]+\-//;
						my $length=length($seq);
						my $seq__=$seq;
						my $seq_original__=$seq_original;
						foreach $strand ("+","-")
						{
							for(my $k=0;$k<3;$k++)
							{
								my $protein="";
								my $protein_original="";
								my $description_="";
								if ($strand=~/\-/)
								{
									my $seq_ = reverse $seq__;
									$seq=$seq_;
									$seq=~tr/ATCG/TAGC/;
									$seq_ = reverse $seq_original__;
									$seq_original=$seq_;
									$seq_original=~tr/ATCG/TAGC/;
									foreach my $key (keys %variants) 
									{ 
										my $key_=$length+$segment_start_first-$key-1; 
										$variants_{$key_}=$variants{$key}; 
										#print qq!$key->$key_\n!;
									}
								}
								else
								{
									$seq=$seq__;
									$seq_original = $seq_original__;
									foreach my $key (keys %variants) { $variants_{$key}=$variants{$key}; }
								}
								for(my $k_=0;$k_<=$k;$k_++)
								{
									$seq=~s/^[A-Z]//;
									$seq_original=~s/^[A-Z]//;
								}
								$modified=0;
								my $stop_found=0;
								for(my $n=0;$n<$length;$n=$n+3)
								{
									my $n_=$n+2;
									my $triplet = substr($seq_original, $n, 3);
									if (length($triplet)==3)
									{
										if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
										$protein_original.=$mapping{$triplet}; 
									}
									$triplet = substr($seq, $n, 3);
									if (length($triplet)==3)
									{
										if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
										$protein.=$mapping{$triplet}; 
										my $triplet_old=$triplet;
										for (my $i=$n, my $j=0; $i<=$n_; $i++,$j++) 
										{
											if ($variants_{$i}=~/^([^#]+)#([^#]+)$/)
											{
												my $old=$1;
												if ($strand=~/\-/) { $old=~tr/ATCG/TAGC/; } 
												substr $triplet_old,$j,1,$old;
											}
										}
										if ($mapping{$triplet_old}!~/[\w\*]/) { $mapping{$triplet_old}="X"; }
										#print qq!$name $n: $triplet_old->$triplet: $mapping{$triplet_old}->$mapping{$triplet}\n!;
										if ($mapping{$triplet}=~/\*/)
										{
											$stop_found=1;
										}
										if ($mapping{$triplet}=~/\*/ and $mapping{$triplet_old}!~/\*/)
										{
											print LOG qq!$n-$n_: $triplet_old->$triplet: $mapping{$triplet_old}->$mapping{$triplet}\n!;
											$description_.=qq!$n-$n_:$triplet_old->$triplet:$mapping{$triplet_old}->$mapping{$triplet},!;
											$modified++;
											$count_variant_in_exon_nonsyn++;
											$count_stop_introduced++;
										}
										else
										{
											if ($mapping{$triplet}!~/\*/ and $mapping{$triplet_old}=~/\*/)
											{
												print LOG qq!$n-$n_: $triplet_old->$triplet: $mapping{$triplet_old}->$mapping{$triplet}\n!;
												$description_.=qq!$n-$n_:$triplet_old->$triplet:$mapping{$triplet_old}->$mapping{$triplet},!;
												$modified++;
												$count_variant_in_exon_nonsyn++;
												$count_stop_removed++;
											}
											else
											{
												if ($mapping{$triplet}!~/^$mapping{$triplet_old}$/)
												{
													print LOG qq!$n-$n_: $triplet_old->$triplet: $mapping{$triplet_old}->$mapping{$triplet}\n!;
													$description_.=qq!$n-$n_:$triplet_old->$triplet:$mapping{$triplet_old}->$mapping{$triplet},!;
													$modified++;
													$count_variant_in_exon_nonsyn++;
												}
											}
										}
									}
								}
								$proteins_count++;
								$protein_modifications_count+=$modified;
								$protein_modifications_distr[$modified]++;
								#$protein=~s/\*$//;
								#$protein=~s/^([^\*]+)\*.*$/$1/;
								#$protein_original=~s/\*$//;
								#$protein_original=~s/^([^\*]+)\*.*$/$1/;
								my $temp1=substr $protein_original,0,100;
								$temp1=~s/^.*\*([^\*]*)$/$1/;
								my $cleavage=0;
								my $temp1_="";
								while($temp1=~s/(.)(.)$//)
								{
									my $aa1=$1;
									my $aa2=$2;
									if ($cleavage<2)
									{
										$temp1_="$aa2$temp1_";
										$temp1.="$aa1";
										if ($aa1=~/^[RK]$/ and $aa2=~/^[^P]$/)
										{
											$cleavage++;
										}
									}
								}
								if ($cleavage==0) { $temp1_=""; }
								my $count_RK = $temp1=~tr/RK/RK/;
								if ($temp1!~s/^[^RK]*[RK]//) { $temp1=""; }
								for(my $l=0;$l<$count_RK-2;$l++) { $temp1=~s/^[^RK]*[RK]/$1/; }
								my $temp2=substr $protein_original,100;
								$temp2=~s/([^\*]*)\*.*$/$1/;
								$cleavage=0;
								my $temp2_="";
								while($temp2=~s/^(.)(.)//)
								{
									my $aa1=$1;
									my $aa2=$2;
									if ($cleavage<2)
									{
										$temp2_.="$aa1";
										$temp2="$aa2$temp2";
										if ($aa1=~/^[RK]$/ and $aa2=~/^[^P]$/)
										{
											$cleavage++;
										}
									}
								}
								if ($cleavage<2) { $temp2_.=$temp2; }
								if ($cleavage==0) { $temp2_=""; }
								if (length($temp1_)>0 and length($temp2_)>0 and length($temp1_)+length($temp2_)>6)
								{
									print OUT qq!>$name-$num_observ-$k$strand-bridge \n$temp1_$temp2_\n!;
									#print OUT qq!>$name-$num_observ-$k$strand-bridge\n$temp1_-$temp2_\n$protein_original\n!;
								}
								if ($modified!=0)
								{
									$proteins_modified++;
									my $temp1=substr $protein,0,100;
									$temp1=~s/^.*\*([^\*]*)$/$1/;
									my $cleavage=0;
									my $temp1_="";
									while($temp1=~s/(.)(.)$//)
									{
										my $aa1=$1;
										my $aa2=$2;
										if ($cleavage<2)
										{
											$temp1_="$aa2$temp1_";
											$temp1.="$aa1";
											if ($aa1=~/^[RK]$/ and $aa2=~/^[^P]$/)
											{
												$cleavage++;
											}
										}
									}
									if ($cleavage==0) { $temp1_=""; }
									my $count_RK = $temp1=~tr/RK/RK/;
									if ($temp1!~s/^[^RK]*[RK]//) { $temp1=""; }
									for(my $l=0;$l<$count_RK-2;$l++) { $temp1=~s/^[^RK]*[RK]/$1/; }
									my $temp2=substr $protein,100;
									$temp2=~s/([^\*]*)\*.*$/$1/;
									$cleavage=0;
									my $temp2_="";
									while($temp2=~s/^(.)(.)//)
									{
										my $aa1=$1;
										my $aa2=$2;
										if ($cleavage<2)
										{
											$temp2_.="$aa1";
											$temp2="$aa2$temp2";
											if ($aa1=~/^[RK]$/ and $aa2=~/^[^P]$/)
											{
												$cleavage++;
											}
										}
									}
									if ($cleavage<2) { $temp2_.=$temp2; }
									if ($cleavage==0) { $temp2_=""; }
									if (length($temp1_)>0 and length($temp2_)>0 and length($temp1_)+length($temp2_)>6)
									{
										print OUT_MOD qq!>$name-$num_observ-$k$strand-bridge-variant $description_\n$temp1_$temp2_\n!;
										#print OUT_MOD qq!>$name-$num_observ-$k$strand-bridge-variant $description_\n$protein\n!;
									}
									my $protein_length=length($protein);
									my $segment_starts_=$segment_starts;
									my $segment_lengths_=$segment_lengths;
									my $segment_starts__="";
									my $segment_lengths__="";
									my $segment_count=0;
									my $length_sum=0;
									my $length_sum_=0;
									my $end_found=0;
									my $start_=$start;
									my $end_=$end;
									my $segment_start_first=0;
									if ($strand=~/\+/)
									{
										while ($segment_starts_=~s/^([0-9\-]+)\,//)
										{
											my $segment_start=$1;
											if ($segment_lengths_=~s/^([0-9]+)\,//)
											{
												my $segment_length=$1;
												$length_sum+=$segment_length;
												if ($end_found==0)
												{
													if ($length_sum>3*$protein_length)
													{
														$end_found=1;
														$segment_length-=$length_sum-3*$protein_length;
														$segment_length+=3;
													}
													$segment_starts__.="$segment_start,";
													$segment_lengths__.="$segment_length,";
													$end_=$start+$segment_start+$segment_length;
													$segment_count++;
												}
											}
										}
									}
									else
									{
										my $segment_start_first=0;
										if ($segment_starts_=~/^([0-9\-]+)\,/) { $segment_start_first=$1; }
										while ($segment_starts_=~s/([0-9\-]+)\,$//)
										{
											my $segment_start=$1;
											if ($segment_lengths_=~s/([0-9]+)\,$//)
											{
												my $segment_length=$1;
												$length_sum+=$segment_length;
												if ($end_found==0)
												{
													my $segment_end=$segment_start+$segment_length;
													if ($length_sum>3*$protein_length)
													{
														$end_found=1;
														$segment_length-=$length_sum-3*$protein_length;
														$segment_length+=3;
													}
													$segment_start-=$segment_start_first;
													$segment_starts__="$segment_start,$segment_starts__";
													$segment_lengths__="$segment_length,$segment_lengths__";
													$start_=$start+$segment_end-$segment_length;
													$segment_count++;
												}
											}
										}
									}
									$segment_lengths__=~s/\,$//;
									$segment_starts__=~s/\,$//;
									print BED qq!$chr\t$start_\t$end_\t$name-variant\t1000\t$strand\t$start_\t$end_\t0\t$segment_count\t$segment_lengths__\t$segment_starts__\n!;
									print LOG qq!---$name\t$chr\t$start_\t$end_\t$name-variant\t1000\t$strand\t$start_\t$end_\t0\t$segment_count\t$segment_lengths__\t$segment_starts__\n!;
								}
							}
						}
					} else { print qq!Error parsing $bed{$name}\n!; }
				}
			} else { print qq!Error in name $chr: $line\n!; }
			
			close(IN);
		}
	}
	print STAT qq!
		$proteins_count proteins
		$proteins_modified proteins modified
		$count_stop_removed stop codons removed
		$count_stop_introduced stop codons introduced
		$protein_modifications_count total modifications
		$count_variant_in_exon variants in exons
		$count_variant_in_exon_old_error variants in exons where old does not match genome
		$count_variant_in_exon_nonsyn non-synonymous variants
		$both_vcf disagreements between the two vcf files
	!;
	print STAT qq!\nnumber of modifications\tnumber of proteins\n!;
	for (my $k=0;$k<=100;$k++) { print STAT qq!$k\t$protein_modifications_distr[$k]\n!; }
	close(OUT);
	close(LOG);
	close(STAT);
	close(BED);
}