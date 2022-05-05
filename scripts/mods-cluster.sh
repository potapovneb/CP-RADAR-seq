#!/bin/bash -l

#$ -S /bin/bash
#$ -N ipdSummary
#$ -pe smp 2
#$ -o /dev/null
#$ -j yes
#$ -cwd

### current processed stratum
stratum=${SGE_TASK_ID}

### define job directory
dirnum=`echo "$stratum / 1000 + 1" | bc`
jobdir=`printf %s/modfiles/%04i $PWD $dirnum`

prefix=`printf %06i $stratum`
whitelist=$prefix.txt
subreads=$prefix.bam
modfile=$prefix.csv

mkdir -p $jobdir
cd $jobdir

### extract list of ZMW's
cat $strata_file | egrep ",$stratum$" | cut -d , -f 2 > $whitelist
cat $whitelist | wc -l

zmwfilter --include $whitelist $aligned_reads $subreads
samtools index $subreads

### detect modifications
ipdSummary $subreads \
    --log-level INFO \
    --numWorkers 2 \
    --reference $reference \
    --csv $modfile \
    --alignmentSetRefWindows > /dev/null 2>& 1

rm $prefix.bam $prefix.bam.bai $prefix.txt

gzip $prefix.csv
