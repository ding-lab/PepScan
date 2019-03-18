# Song Cao, Ding Lab
# scao@genome.wustl.edu
my $usage='perl /gscmnt/gc2108/info/medseq/proteogenomics/analyses/pipeline_step_1_peptides_detection/cptac/breast $out_file
'; 

die $usage unless @ARGV == 2;
my ( $run_dir, $out_dir) = @ARGV;

my $out_file; 

my %count_p=();
my @gene_list; 
my @sample_list; 
 
if ($run_dir =~/(.+)\/$/) {
        $run_dir = $1;
}

if(!(-d $out_dir)) {`mkdir $out_dir`; }

my $run_dir_msgf=$run_dir."/proteome"; 


if(-e $out_file)
{
#$print $out_file,"\n"; 
my $com="rm $out_file";
print $com;
system($com);
#<STDIN>;
 }

opendir(DH, $run_dir_msgf) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list_msgf = readdir DH;
close DH;

for (my $i=0;$i<@sample_dir_list_msgf;$i++) {
   $sample_name = $sample_dir_list_msgf[$i];
   if (!($sample_name =~ /\./)) {
   $sample_full_path = $run_dir_msgf."/".$sample_name;
   print $sample_full_path,"\n";
   if (-d $sample_full_path)
    {
    opendir(DH,$sample_full_path);
    my @sample_mzmlgz = readdir DH;
    {
     foreach my $tsv_file_name (@sample_mzmlgz)
     {
      if (($tsv_file_name =~ /\.tsv$/))
        {
        $sample_full_name=$run_dir_msgf."/".$sample_name."/".$tsv_file_name;
	print $sample_full_name,"\n"; 
   	&read_tsv_file($sample_full_name); 
	     }

      }
     }
     }
     }
     close DH;
    }

### reading variants ##

sub read_tsv_file() 
 {

   my ($stsv)=@_; 
   open(IN,"<$stsv"); 
   open(OUT,">>$out_file");
   my $tsv_name=(split(/\//,$stsv))[-1]; 
   my %scans=();
 
   while(<IN>)
   {
     chomp; 
     my $line=$_; 
     my $nn;
     if($line=~/scan=(\d+)/) {
     $nn=$1; 
     #print $nn,"\n"; <STDIN>;
     if(!defined $scans{$nn})
            {
     if($line=~/variant/ && !($line=~/XXX_/)) 
	  { 
 	   
	    my @lines=split(/\t/,$line);
	    my @proteins=split(/;/,$lines[9]);
	    my $prosize=scalar @proteins; 
	    for(my $i=0;$i<$prosize;$i++)
	    {
	    my @genes=split(/-/,$proteins[$i]);
            if($genes[2]=~/variant/)
	    {
            #print @genes,"\n"; 
	    my $gene=$genes[0]."-".$genes[1]."-".$genes[2];  
 	    my $sample=join("-",$genes[4],$genes[5]);  
	    $sample=~s/\(pre(\S+)//g; 
	    my $peptide=$lines[8];
	    #print $peptide,"\n"; 
	    $peptide=~s/\+//g;
	    $peptide=~s/\.//g;
	    $peptide =~s/\d//g;
	    my $evalue=$lines[12];
	    if($evalue<=3e-9 && $gene ne "") { 
		print OUT $sample,"\t", $gene, "\t", $peptide,"\t",$tsv_name,"\n"; } 
		last;    
		}
	    #last; 
	   } 
   	}
	$scans{$nn}=1; 
	}	
 	}
	}
} 
