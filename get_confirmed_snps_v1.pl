
# get SNPs confirmed by MS spectroscopy 
# Song Cao, Ding Lab
# scao@genome.wustl.edu

#v1, select the longest peptide as the representive peptide, 1-17-2014

#!/usr/bin/perl
use strict;
use warnings;
my $usage='
perl get_confirmed_snps.pl $run_dir $db_peptide $out_dir
'; 
die $usage unless @ARGV == 3;
my($run_dir, $db_peptide, $out_dir) = @ARGV;
my %count_p=();
my @gene_list; 
my @sample_list; 
my %tot_gsnp=(); 
my %tot_ssnp=(); 
my %num_expr_gsnp=(); 
my %num_expr_ssnp=(); 
my %table_peptide=(); 
my %gene_expr_gsnp=(); 
my %gene_expr_ssnp=(); 
my %gene_list_gsnp=(); 
my %gene_list_ssnp=(); 
my %variant_list_ssnp=(); 
my %variant_list_gsnp=();
my %peptide_list_ssnp=(); 
my %peptide_list_gsnp=(); 
my %num_pep_gsnp=(); 
my %num_pep_ssnp=(); 
my @gene_list_germline=(); 
my @gene_list_somatic=(); 

if ($run_dir =~/(.+)\/$/) {
        $run_dir = $1;
}

if(!(-d $out_dir)) { `mkdir $out_dir`; }

my $run_dir_quilts;
$run_dir_quilts=$run_dir."/quilts";
opendir(DH, $run_dir_quilts) or die "Cannot open dir $run_dir"; 

my @sample_dir_list_quilts = readdir DH;
close DH;
open(IDB,"<$db_peptide"); 
while(<IDB>) 
 {
   chomp; 
   my $line=$_; 
   my @list=split("\t",$line);
   my $find_sample=0;
   if(scalar @sample_list==0) { $sample_list[0]=$list[0]; }
   else { for(my $i=0;$i<=scalar @sample_list-1;$i++) {
          if($list[0] eq $sample_list[$i]) {
          $find_sample=1;
          last;
           }
          }

  if($find_sample==0) {
        $sample_list[scalar @sample_list]=$list[0]; }
      }
 
  if(defined $table_peptide{$list[0]}{$list[1]}) { 
						   my $new_peptide=join(";",$table_peptide{$list[0]}{$list[1]},$list[2]);  
						   $table_peptide{$list[0]}{$list[1]}=$new_peptide;
                     						}
  						
  else { $table_peptide{$list[0]}{$list[1]}=$list[2]; }
 
}

close IDB; 


print scalar @sample_list,"\t",@sample_list,"\n";
#for(my $i=0;$i<=scalar @sample_list -1;$i++) 
 #{
  #  print 
 #}

for (my $i=0;$i<@sample_dir_list_quilts;$i++) {
my $sample_name = $sample_dir_list_quilts[$i];
my $sample_name_abr=$sample_name; 
$sample_name_abr=~s/TCGA-//g;
 
my $find_sample=0; 
for(my $i=0;$i<=scalar @sample_list-1;$i++) 
   {
	if($sample_name_abr eq $sample_list[$i]) 
	 {
          $find_sample=1;
          last;
         }
    }

if (!($sample_name =~ /\./) && $find_sample==1)
     {
       my $sample_full_path = $run_dir_quilts."/".$sample_name;
       print $sample_full_path,"\n";
       if (-d $sample_full_path)
       {
       my $file_vcf=$sample_full_path."/vcf/".$sample_name.".vcf";
       my $file_fasta=$sample_full_path."/protein/fasta/QUILTSv1.0-ensembl_human_37.70-var-vcf-rna-junction/variants/proteome.fasta";
       &count_tot($file_vcf,$file_fasta,$sample_name_abr);               
       }
      }
     }

my $out1=$out_dir."/germline_snp_confirmed_by_ms.txt"; 
my $out2=$out_dir."/somatic_snp_confirmed_by_ms.txt"; 
open(OUT1,">$out1");
open(OUT2,">$out2"); 
#print OUT1 "Sample","\t", "# of tot missense variants (germline)","\t","# of snp confirmed by ms", "\t", "Ratio", "\n"; 

#foreach my $sample(@sample_list)
 #{ 
  # if(defined $num_expr_gsnp{$sample})
   #{
   # print OUT1 $sample,"\t", $tot_gsnp{$sample},"\t", $num_expr_gsnp{$sample},"\t", $num_expr_gsnp{$sample}/$tot_gsnp{$sample}, "\n";}
   #else { print OUT1 $sample,"\t", $tot_gsnp{$sample},"\t", "0", "\t", "0", "\n"; } 
 #}

print OUT1 "Sample","\t", "# of tot missense variants (germline)","\t","# of variants confirmed by ms", "\t", "Ratio (\%)", "\t", "name of confirmed variants", "\t", "gene name", "\t", "sequence of one respective peptide","\t","number of spectra","\n";

foreach my $sample(@sample_list)
 {
   if(defined $num_expr_gsnp{$sample})
   {
   my $glist="";
   my $vlist="";
   my $plist="";
   my $nlist="";
   foreach my $variant (sort {$num_pep_gsnp{$sample}{$b} <=> $num_pep_gsnp{$sample}{$a}} keys %{$num_pep_gsnp{$sample}})
    {
    if($vlist eq "") { $vlist=$variant; }
    else { $vlist=join(";",$vlist,$variant); }
    if($glist eq "") { $glist=$gene_list_gsnp{$sample}{$variant}; }
    else { $glist=join(";",$glist,$gene_list_gsnp{$sample}{$variant}); }
    if($plist eq "") { $plist=$peptide_list_gsnp{$sample}{$variant}; }
    else { $plist=join(";",$plist,$peptide_list_gsnp{$sample}{$variant}); }
    if($nlist eq "") { $nlist=$num_pep_gsnp{$sample}{$variant}; }
    else { $nlist=join(";",$nlist,$num_pep_gsnp{$sample}{$variant}); }
    }
    print OUT1 $sample,"\t", $tot_gsnp{$sample},"\t", $num_expr_gsnp{$sample},"\t",$num_expr_gsnp{$sample}/$tot_gsnp{$sample}*100, "\t", $vlist, "\t", $glist,"\t",$plist,"\t",$nlist,"\n";
   }
   else { print OUT1 $sample,"\t", $tot_ssnp{$sample},"\t", "0", "\t", "0","\t", "None", "\t", "None", "\t", "None","\t","None","\n"; }
  }

print OUT2 "Sample","\t", "# of tot missense variants (somatic)","\t","# of variants confirmed by ms", "\t", "Ratio (\%)", "\t", "name of confirmed variants", "\t", "gene name", "\t", "sequence of one respective peptide","\t","number of spectra","\n";

foreach my $sample(@sample_list)
 {
   if(defined $num_expr_ssnp{$sample})
   {
   my $glist=""; 
   my $vlist=""; 
   my $plist=""; 
   my $nlist="";  
   foreach my $variant (sort {$num_pep_ssnp{$sample}{$b} <=> $num_pep_ssnp{$sample}{$a}} keys %{$num_pep_ssnp{$sample}})
    {
    if($vlist eq "") { $vlist=$variant; }
    else { $vlist=join(";",$vlist,$variant); }
    if($glist eq "") { $glist=$gene_list_ssnp{$sample}{$variant}; }
    else { $glist=join(";",$glist,$gene_list_ssnp{$sample}{$variant}); }
    if($plist eq "") { $plist=$peptide_list_ssnp{$sample}{$variant}; }
    else { $plist=join(";",$plist,$peptide_list_ssnp{$sample}{$variant}); }
    if($nlist eq "") { $nlist=$num_pep_ssnp{$sample}{$variant}; }
    else { $nlist=join(";",$nlist,$num_pep_ssnp{$sample}{$variant}); }
    }
    print OUT2 $sample,"\t", $tot_ssnp{$sample},"\t", $num_expr_ssnp{$sample},"\t",$num_expr_ssnp{$sample}/$tot_ssnp{$sample}*100, "\t", $vlist, "\t", $glist,"\t",$plist,"\t",$nlist,"\n";  
   }
   else { print OUT2 $sample,"\t", $tot_ssnp{$sample},"\t", "0", "\t", "0","\t", "None", "\t", "None", "\t", "None","\t","None","\n"; }
  }

open(OUT3,">somatic_genes_list_by_ms.txt");

print "tot number of somatic genes:", scalar @gene_list_somatic, "\n";

for(my $i=0;$i<=scalar @gene_list_somatic-1;$i++) {
print OUT3 $gene_list_somatic[$i],"\n";
print $gene_list_somatic[$i],"\t"; };
print "\n"; 

open(OUT4,">germline_genes_list_by_ms.txt"); 

print "tot number of germline genes:", scalar @gene_list_germline, "\n"; 
for(my $i=0;$i<=scalar @gene_list_germline-1;$i++) { 
print OUT4 $gene_list_germline[$i],"\n"; 
print $gene_list_germline[$i],"\t"; }; 

sub count_tot()
 {
   # $file1, the original vcf file
   # $file2, the variant file from quilts   
   my($file1,$file2,$sample)=@_; 
   my %variant=(); 
   #open(IN1,"<$file1"); 
   #close IN1; 
   open(IN2,"<$file2");
   my %snp_to_gid=(); 
   my %gid_protein=(); 
   my $gid;
   while(<IN2>) 
   {
    chomp; 
    my $line=$_; 
    if($line=~/\>/) { my $pos_variant=index($line,"variant");
		      my $pos_var=index($line,"VAR:"); 
		      $gid=substr($line,1,$pos_variant+6); 
		      my $s2=substr($line,$pos_var+4,length($line)-$pos_var-6);
		      my @s2a=split(/,/,$s2); 
		      for(my $i=0;$i<=scalar @s2a-1;$i++) { $snp_to_gid{$s2a[$i]}=$gid;}
	             }
    else 
	{
	 $gid_protein{$gid}=$line; 
	}
   }
  close IN2; 

  my %haveg_var1; 
  my %haves_var1; 
  my %haveg_var2;
  my %haves_var2;

  foreach my $var (keys %snp_to_gid)
  {

    my @t=split(":",$var);
    #print $t[0],"\t",$t[1],"\n";

    if($var=~/-G:/ && !defined $haveg_var1{$t[0]}) {  $tot_gsnp{$sample}++; $haveg_var1{$t[0]}=1; }
    if($var=~/-S:/ && !defined $haves_var1{$t[0]}) {  $tot_ssnp{$sample}++; $haves_var1{$t[0]}=1; } 

    #my @t=split(":",$var); 
    my $aa=$t[1];
    my $aa_pos;
    if($aa=~/(\d+)/) { $aa_pos=$1; }
    my $gid=$snp_to_gid{$var};
  #  print $sample,"\t",$var, "\t", $aa,"\t",$aa_pos,"\t",$gid,"\n"; 
   # <STDIN>;

    if($var=~/-G:/)
    {
      if(defined $table_peptide{$sample}{$gid})
       {
        my @peptide=split(/;/,$table_peptide{$sample}{$gid});
        my $count_gsnp=0;
        for(my $i=0;$i<=scalar @peptide-1;$i++)
          {
            my $pt=$peptide[$i];
            my $pos=index($gid_protein{$gid},$pt);
           # print "G", "\t", $pt, "\t", $pos,"\t", $aa_pos,"\n";  <STDIN>;
            my $len_pt=length($pt);
	    if($pos+1<=$aa_pos && $pos+$len_pt>=$aa_pos) {
                                                                 if($count_gsnp==0)
                                                                     {
                                                                        if(!defined $haveg_var2{$t[0]}) { $num_expr_gsnp{$sample}++; $haveg_var2{$t[0]}=1; }
                                                                        $gene_list_gsnp{$sample}{$var}=$gid;
                                                                        $peptide_list_gsnp{$sample}{$var}=$pt;
                                                                        my @genes=split(/\-/,$gid);
                                                                        my $gene_list_size=scalar @gene_list_germline;
                                                                        my $find_gene=0;
                                                                        if($gene_list_size==0) { $gene_list_germline[0]=$genes[1]; }
                                                                        else { for(my $i=0;$i<=$gene_list_size-1;$i++) {
                                                                                if($genes[1] eq $gene_list_germline[$i]) {
                                                                                        $find_gene=1;
                                                                                        last;
                                                                                        }
                                                                                }
                                                                                if($find_gene==0) {
                                                                                 $gene_list_germline[scalar @gene_list_germline]=$genes[1]; }
                                                                             }
                                                                         $count_gsnp=1;
                                                                        # $haveg_var2{$t[0]}=1; 
                                                                        }

								       if(defined $peptide_list_gsnp{$sample}{$var}) { my $len_curr=length($peptide_list_gsnp{$sample}{$var});

                                                                                               if($len_pt>$len_curr) { $peptide_list_gsnp{$sample}{$var}=$pt;} }
                                                                        $num_pep_gsnp{$sample}{$var}++;
							 }
        }
       }
    }

    if($var=~/-S:/)
    {
     # print $var,"\n";
      if(defined $table_peptide{$sample}{$gid})
       {
        my @peptide=split(/;/,$table_peptide{$sample}{$gid});
        my $count_ssnp=0;
        for(my $i=0;$i<=scalar @peptide-1;$i++)
          {
            my $pt=$peptide[$i];
            my $pos=index($gid_protein{$gid},$pt);
           # print $sample,"\t","S:", "\t", $pt, "\t", $pos,"\t", $aa_pos,"\n";  <STDIN>;
            my $len_pt=length($pt);

            if($pos+1<=$aa_pos && $pos+$len_pt>=$aa_pos) {
            #print $var,"\t",$sample,"\t","S:", "\t", $pt, "\t", $pos,"\t", $t[0],"\t",$aa_pos,"\n";  <STDIN>;
                                                                 if($count_ssnp==0)    
                                                                 {

                                                                        if(!defined $haves_var2{$t[0]}) { $num_expr_ssnp{$sample}++; $haves_var2{$t[0]}=1; }

                                                                        $gene_list_ssnp{$sample}{$var}=$gid;
                                                                        $peptide_list_ssnp{$sample}{$var}=$pt;
                                                                        my @genes=split(/\-/,$gid);
                                                                        my $gene_list_size=scalar @gene_list_somatic;
                                                                        my $find_gene=0;
                                                                        if($gene_list_size==0) { $gene_list_somatic[0]=$genes[1]; }
                                                                        else { for(my $i=0;$i<=$gene_list_size-1;$i++) {
                                                                                if($genes[1] eq $gene_list_somatic[$i]) {
                                                                                        $find_gene=1;
                                                                                        last;
                                                                                        }
                                                                                }
                                                                                if($find_gene==0) {
                                                                                 $gene_list_somatic[scalar @gene_list_somatic]=$genes[1]; }
                                                                             }
                                                                         $count_ssnp=1;
                                                                        }

                                                                       if(defined $peptide_list_ssnp{$sample}{$var}) { 
								       my $len_curr=length($peptide_list_ssnp{$sample}{$var});
								      # print $len_curr,"\n";
                                                                       if($len_pt>$len_curr) { $peptide_list_ssnp{$sample}{$var}=$pt;} }
                                                                       
                                                                       $num_pep_ssnp{$sample}{$var}++;
                                                         }
        }
       }
    }

   }
}
