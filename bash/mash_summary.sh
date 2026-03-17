#!/bin/bash
#===============================================================================
# Name         : s_cat_flye_info_files.sh
# Usage        : sbatch s_cat_flye_info_files.sh (run after s_extract_flye_info.sh)
#===============================================================================
#
#Submission script for HTCF
#SBATCH --time=0-6:00:00 # days-hh:mm:ss
#SBATCH --job-name=mash_cat
#SBATCH --mem=2G
#SBATCH --cpus-per-task=4
#SBATCH --output=slurm_out/mash/z_combine_mash_%a.out
#SBATCH --error=slurm_out/mash/z_combine_mash_%a.out

#set up directories
basedir="$PWD"
mash_in="/scratch/gdlab/hall.l.r/projects/crescent/other_india_samples/d05_mash"
outdir="${basedir}/d05_mash_external"

mash_comb_out="${outdir}/CRESCENT_external_mash.txt"
header="Sample\tMash Identity\tShared Hashes\tP-value\tQuery-ID\tContaminant Hashes\n"
echo -n -e ${header} >> ${mash_comb_out}

for file in ${mash_in}/*.tab; do
        # Extract the base name without extension (contig name)
    sample_name=$(basename "$file" ".tab")
    # Extract only fields with a shared hash greater than or equal to 700
    while IFS= read -r line; do
        echo -n -e "${sample_name}\t${line}\n" >> ${mash_comb_out}
    done < ${file}
done

echo "Done combining all mash output files"
