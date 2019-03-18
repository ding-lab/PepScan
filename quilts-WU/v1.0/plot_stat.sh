module load python/2.7

#cd /ifs/data/proteomics/tcga/samples/breast/stat
# (rjm)
BASE="/gscuser/rmashl/gsc/working/cptac/test_quilts"
cd $BASE/samples/breast/stat


rm -f *.png
#python /ifs/data/proteomics/tcga/scripts/quilts/v1.0/plot_stat_variants.py
#python /ifs/data/proteomics/tcga/scripts/quilts/v1.0/plot_stat_seq_length.py
# (rjm)
python $BASE/scripts/quilts-WU/v1.0/plot_stat_variants.py
python $BASE/scripts/quilts-WU/v1.0/plot_stat_seq_length.py
