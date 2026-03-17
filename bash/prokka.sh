
#!/bin/bash
#===============================================================================
# File Name    : prokka.sh
# Description  : Annotation of assembled bacterial genomes or contigs
# Usage        : sbatch prokka.sh
#===============================================================================
#Submission script for HTCF
#SBATCH --time=1-00:00:00 # days-hh:mm:ss
#SBATCH --job-name=prokka
#SBATCH --array=1-140%20
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm_out/prokka/z_prokka_%A_%a.out
#SBATCH --error=slurm_out/prokka/z_prokka_%A_%a.out

eval $( spack load --sh prokka@1.14.6 )

basedir="$PWD"
indir="${basedir}/d03_spades"
outdir="${basedir}/d06_prokka"

mkdir -p ${outdir}
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/seq_list.txt`

set -x
time prokka ${indir}/${sample}/scaffolds.fasta \
            --outdir ${outdir}/${sample} \
            --prefix ${sample} \
            --mincontiglen 200 \
            --addgenes \
            --force \
            --cpus ${SLURM_CPUS_PER_TASK}
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred in ${sample}!"
  exit $RC
fi
