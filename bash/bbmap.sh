#!/bin/bash
#===============================================================================
# Name         : s15_bbmap.sh
# Description  : This script will run bbmap to estimate coverage for an assembly
# Usage        : sbatch s15_bbmap.sh
# Author       : Luke Diorio-Toth, ldiorio-toth@wustl.edu
# Version      : 1.2
# Created On   : Fri Apr 17 15:40:47 CDT 2020
# Modified On  : 2023-04-28 by Emily Benedict, ebenedict@wustl.edu
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=1-0:00:00 # days-hh:mm:ss
#SBATCH --job-name=bbmap
#SBATCH --array=1-140%20
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm_out/bbmap/z_bbmap_%a_%A.out
#SBATCH --error=slurm_out/bbmap/z_bbmap_%a_%A.out

eval $( spack load --sh bbmap@38.63 )

basedir="$PWD"
readsdir="${basedir}/d01_clean_reads"
refdir="${basedir}/d03_spades"
outdir="${basedir}/d04_assembly_qc/bbmap"

# Read in the slurm array task
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/seq_list.txt`

# Make output directory
mkdir -p ${outdir}/${sample}

set -x
time bbwrap.sh \
    in1=${readsdir}/${sample}_FW_clean.fastq.gz,null \
    in2=${readsdir}/${sample}_RV_clean.fastq.gz,null \
    ref=${refdir}/${sample}/scaffolds.fasta \
    nodisk=t \
    covstats=${outdir}/${sample}/scaffold_covg_stats.txt \
    t=${SLURM_CPUS_PER_TASK} \
    append &> ${outdir}/${sample}/summ_covg_stats.txt
RC=$?
set +x

if [ $RC -eq 0 ]
then
    echo "Job completed successfully"
else
    echo "Error Occurred in ${sample}!"
    exit $RC
fi
