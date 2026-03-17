#!/bin/bash

#===============================================================================
# Name         : raxml.sh
# Description  : Infers a maximum-likelihood phylogenetic tree
#                from alignments of nucleotide or protein sequences.
# Usage        : sbatch raxml.sh
#===============================================================================

#SBATCH --job-name=all_raxml
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --output=slurm_out/raxml/z_raxml_%A.out
#SBATCH --error=slurm_out/raxml/z_raxml_%A.out

eval $( spack load --sh raxml+pthreads@8.2.12 )

basedir="$PWD"
indir="${basedir}/d16_panaroo_bakta/all_external_and_internal"
outdir="${indir}/d14_raxml"

mkdir -p ${outdir}

set -x
time raxmlHPC-PTHREADS \
    -s ${indir}/core_gene_alignment.aln \
    -w ${outdir} \
    -n raxml_core_genome_tree \
    -m GTRGAMMA \
    -f a \
    -T ${SLURM_CPUS_PER_TASK} \
    -N 100 \
    -p 12345 \
    -x 54321
RC=$?
set +x

if [ $RC -eq 0 ]
then
  echo "Job completed successfully"
else
  echo "Error occurred!"
  exit $RC
fi
