#!/bin/bash -l

### path to RADAR-seq repo
RADARSEQ_DIR="/PATH/TO/RADARSEQ_INSTALLATION_DIRECTORY"

export PATH=$RADARSEQ_DIR/bin:$PATH

###############################################################################
### Download example PacBio BAM subreads file                               ###
###############################################################################

wget https://sra-download.ncbi.nlm.nih.gov/traces/sra70/SRZ/019043/SRR19043774/lib26.bam

###############################################################################
### Define input data files and output directory                            ###
###############################################################################

# input PacBio BAM subreads file
subreads=$PWD/lib26.bam

# reference FASTA file
reference=$RADARSEQ_DIR/references/tko.fasta

# output directory
rundir=$PWD/radarseq_output

echo ""
echo "subreads   : $subreads"
echo "reference  : $reference"
echo "output dir : $rundir"

###############################################################################
### Demultiplexing (optional)                                               ###
###############################################################################

# echo ""
# echo "Demultiplexing reads"

# ### create output directory for demultiplexed reads
# jobdir=$rundir/00-demux
# mkdir -p $jobdir
# cd $jobdir

# ### demultiplex sequencing data
# lima $subreads $barcodes movie.bam \
#     --same \
#     --split-bam \
#     --num-threads 16 \
#     --log-level TRACE \
#     --log-file movie.log

###############################################################################
### Mapping                                                                 ###
###############################################################################

echo ""
echo "Mapping reads"

### create output directory for mapped reads
jobdir=$rundir/01-mapping
mkdir -p $jobdir
cd $jobdir

### align reads
time TMPDIR=$PWD pbmm2 align $reference $subreads aligned_reads.bam \
    --sort \
    --min-concordance-perc 75.0 \
    --min-length 500 \
    --sample \"\" \
    --num-threads 16 \
    --log-level TRACE \
    --log-file aligned_reads.log

### index reads
time pbindex aligned_reads.bam

### find non-overlapping reads
time radarseq-strata.pl --prefix aligned_reads aligned_reads.bam

###############################################################################
### Modification detection                                                  ###
###############################################################################

echo ""
echo "Modification detection"

### result files from the previous step
aligned_reads="$rundir/01-mapping/aligned_reads.bam"
strata_file="$rundir/01-mapping/aligned_reads.strata.csv"

num_strata=`cat $strata_file | cut -d, -f5 | sort -g | tail --lines 1`

### create output directory for modification detection files
jobdir=$rundir/02-mods
mkdir -p $jobdir
cd $jobdir

### by default, use parallel to execute multiple jobs
### '-j 16' defines the number of jobs running in parallel
### adjust according to the number of CPU cores in your computer divided by two (two CPU threads per job)
time seq 1 $num_strata | parallel -j 16 -k $RADARSEQ_DIR/scripts/mods-parallel.sh '{}' $strata_file $reference $aligned_reads

# ### if cluster envrionment is available, use the script below to submit array job on the cluster
# # qsub -t 1-$num_strata -v strata_file=$strata_file,reference=$reference,aligned_reads=$aligned_reads $RADARSEQ_DIR/scripts/mods-cluster.sh

###############################################################################
### Patch detection                                                         ###
###############################################################################

echo ""
echo "Patch detection"

### directory containing modification files from the previous step
moddir=$rundir/02-mods/modfiles

### create output directory for patches
jobdir=$rundir/03-patches
mkdir -p $jobdir
cd $jobdir

### find all modfiles
find $moddir -name "*.csv.gz" | sort > modfiles.txt

### detect patches
radarseq-detect.pl --msfile $reference.methylation $reference modfiles.txt

### filter & convert patches
radarseq-feat.pl patches.csv > patches_svm_features.txt

svm-scale -r $RADARSEQ_DIR/svm_model/model.scale patches_svm_features.txt > patches_svm_features_scaled.txt
svm-predict patches_svm_features_scaled.txt $RADARSEQ_DIR/svm_model/model.weights patches_svm_prediction.tmp

### merge patch features & SVM predictions
(echo "SVM_Prediction" && cat patches_svm_prediction.tmp) > patches_svm_prediction.txt
paste -d, patches.csv patches_svm_prediction.txt > patches_scored.csv

### generate final list
radarseq-table.pl patches_scored.csv > radarseq_patches.csv

### remove temporary files
rm -f patches*
rm -f features*

###############################################################################
### Results                                                                 ###
###############################################################################

cp $rundir/03-patches/radarseq_patches.csv $rundir/
