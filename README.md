# RADAR-seq

This is an official repository for the RADAR-seq computational workflow as adapted for the Current Protocols publication. Please see citing information below.

# Prerequisites

The RADAR-seq computational workflow requires a number of tools to be installed and available from the command line in your system.

## PacBio tools

* ```lima``` (optional) - Demultiplexing PacBio data. Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbbioconda).
* ```pbmm2``` - minimap2 SMRT wrapper for PacBio data. Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbbioconda).
* ```pbindex``` - Indexing PacBio sequencing data. Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbbioconda).
* ```zmwfilter``` - Simple utility for filtering PacBio BAM data on ZMW ID(s). Available as part of [pbbioconda](https://github.com/PacificBiosciences/pbbioconda).
* ```ipdSummary``` - Tool for detecting DNA base-modifications from kinetic signatures. Available as part of PacBio [SMRT Link](https://www.pacb.com/support/software-downloads/).

## Third-party tools

* ```samtools``` - Handling high-throughput sequencing data. Available through [conda](https://anaconda.org/bioconda/samtools) or [GitHub](https://github.com/samtools/samtools).
* ```svm-scale```, ```svm-predict``` - LIBSVM tools (Library for Support Vector Machines). Available through [conda](https://anaconda.org/conda-forge/libsvm) or [GitHub](https://github.com/cjlin1/libsvm).
* ```parallel``` - Tool for executing jobs in parallel using one or more computers. Available through [conda](https://anaconda.org/conda-forge/parallel).

## Custom PERL scripts

These custom scripts are provided as part of this GitHub repository. The scripts must be made available from the command line in your system (for example, by adding the scripts directory to ```$PATH``` environment variable).

* ```radarseq-strata.pl``` - Identify sets of non-overalling sequencing reads.
* ```radarseq-detect.pl``` - Detect patches of modified bases.
* ```radarseq-feat.pl``` - Extract patch features for subsequent SVM classification.

# Computational workflow

Processing PacBio sequencing data in the RADAR-seq computational workflow proceeds through a series of steps, where output of one step serves as input for the next step. As a result, the RADAR-seq workflow produces a table (```radarseq_patches.csv```) providing the genomic strand, start and stop location of the detected patch, and the DNA sequence of the detected patch.

The ```scripts/``` directory contains an example workflow (```workflow.sh```) that can be executed from the command line. Before executing the example workflow, make sure that the ```$RADARSEQ_DIR``` environment variable in ```workflow.sh``` points to the RADAR-seq installation directory on your computer.

The computational steps are described below.

## Example PacBio BAM sequencing data

For the purpose of this example, you can download PacBio sequencing data for *Thermococcus kodakarensis* RNaseH2 knockout libraries.

```
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra70/SRZ/019043/SRR19043774/lib26.bam
```

For convenience, we define a set of environment variables to define location of the input data files and the output directory. Please update these variables according to your computational environment.

Input PacBio BAM subreads file:
```
subreads=$PWD/lib26.bam
```

Reference FASTA file:
```
reference=$RADARSEQ_DIR/references/tko.fasta
```

Output directory:
```
rundir=$PWD/radarseq_output
```

This repository comes with the reference FASTA files for *Escherichia coli* strain MG1655 (```ecoli.fasta```) and *Thermococcus kodakarensis* (```tko.fasta```).

## Demultiplexing (optional)

Multiplexed PacBio sequencing data need to be demultiplexed first. This can be done with ```lima``` command-line tool. See usage example below:

```
lima --same --split-bam subreadset.xml barcodes.fasta movie.bam
```

```subreadset.xml``` is the PacBio XML dataset file (or PacBio BAM subreads file), ```barcodes.fasta``` is the list of FASTA barcodes used in multiplexed sequencing, and ```movie.bam``` is the output file.

This tool will produce a separate PacBio BAM subreads file for each barcode pair. Users have to identify the correct output PacBio BAM subreads file and use as input for the subsequent processing steps below.

More information on PacBio multiplexing kits can be found at https://www.pacb.com/products-and-services/consumables/multiplexing-kits/.

We skip this in the example workflow since the downloaded PacBio sequencing data is already demultiplexed.

## Mapping

PacBio sequencing reads are aligned to the reference genome with ```pbmm2``` command:

```
### create output directory for mapped reads
jobdir=$rundir/01-mapping
mkdir -p $jobdir
cd $jobdir

### align reads
TMPDIR=$PWD pbmm2 align $reference $subreads aligned_reads.bam \
    --sort \
    --min-concordance-perc 75.0 \
    --min-length 500 \
    --sample \"\" \
    --num-threads 16 \
    --log-level TRACE \
    --log-file aligned_reads.log

### index reads
pbindex aligned_reads.bam

### find non-overlapping reads
radarseq-strata.pl --prefix aligned_reads aligned_reads.bam
```

This step produces ```aligned_reads.bam```, which contains sequencing reads aligned to the reference genome. An additional file is ```aligned_reads.strata.csv```, which defines sets of non-overlapping sequencing reads to speed up modification detection in the next step.

## Modification detection

In this step, single-molecule reads are iteratively extracted from ```aligned_reads.bam``` and PacBio ```ipdSummary``` tool is run to output PacBio CSV modification detection data. The modification files are organized into a set of directories. The iterative processing is done by the ```mods-parallel.sh``` script. To speed up calculations, the ```parallel``` utility is used to run several jobs simultaneously. This step is computationally demanding and takes ~9 hours on a single computer with 32 CPU cores.

```
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
seq 1 $num_strata | parallel -j 16 -k $RADARSEQ_DIR/scripts/mods-parallel.sh '{}' $strata_file $reference $aligned_reads
```

As an alternative, when a computer cluster is available, a user can submit modificatino detection jobs using an example script as follows:

```
### if cluster envrionment is available, use the script below to submit array job on the cluster
# qsub -t 1-$num_strata -v strata_file=$strata_file,reference=$reference,aligned_reads=$aligned_reads $RADARSEQ_DIR/scripts/mods-cluster.sh

```

In our hands, a computer cluster with 300 CPU cores completes all modification detection in ~1.5 hours.

## Patch detection

In the final step, the custom scripts are used to identify stretches of modified bases (patches) based on the dereived modification data in the previous step.

```
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
```

The derived patches are scored using a trained SVM model.
```
### filter & convert patches
radarseq-feat.pl patches.csv > patches_svm_features.txt

svm-scale -r $RADARSEQ_DIR/svm_model/model.scale patches_svm_features.txt > patches_svm_features_scaled.txt
svm-predict patches_svm_features_scaled.txt $RADARSEQ_DIR/svm_model/model.weights patches_svm_prediction.tmp

### merge patch features & SVM predictions
(echo "SVM_Prediction" && cat patches_svm_prediction.tmp) > patches_svm_prediction.txt
paste -d, patches.csv patches_svm_prediction.txt > patches_scored.csv
```

Then the resulting table providing the genomic strand, start and stop location of the detected patch, and the DNA sequence of the detected patch is produced.
```
### generate final list
radarseq-table.pl patches_scored.csv > radarseq_patches.csv
```

## Additional notes

# Citing

* Kelly M. Zatopek, Vladimir Potapov, Jennifer L. Ong, and Andrew F. Gardner (2022). Utilizing RADAR-seq to detect and quantitate DNA damage on a genome-wide scale. (submitted to Current Protcols).

* Kelly M. Zatopek, Vladimir Potapov, Lisa L. Maduzia, Ece Alpaslan, Lixin Chen, Thomas C. Evans Jr., Jennifer L. Ong, Laurence M. Ettwiller, Andrew F. Gardner (2019). RADAR-seq: A RAre DAmage and Repair sequencing method for detecting DNA damage on a genome-wide scale. DNA repair (Amst), 80:36-44. doi: [10.1016/j.dnarep.2019.06.007](https://dx.doi.org/10.1016/j.dnarep.2019.06.007).
