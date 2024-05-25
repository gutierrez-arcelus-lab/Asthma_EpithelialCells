#!/bin/bash
#SBATCH --partition=bch-compute-pe
#SBATCH -n 15
#SBATCH -N 1
#SBATCH --mem=2GB
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --time=48:00:00
#SBATCH --array=1-22
bedfile=$1

source /programs/biogrids.shrc
module load anaconda2
cd 	ldsc

fileprefix=${bedfile%".bed"}
chromnumber=$SLURM_ARRAY_TASK_ID
echo $chromnumber
annot_command="python.pybedtools make_annot.py --bed-file $bedfile --bimfile 1000G_EUR_Phase3/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chromnumber}.bim --annot-file $fileprefix.${chromnumber}.annot"
echo $annot_command
$annot_command
ld_command="python ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${chromnumber} --ld-wind-cm 1 --annot $fileprefix.${chromnumber}.annot --thin-annot --out $fileprefix.${chromnumber} --print-snps hapmap3_snps/hm.${chromnumber}.snp"
echo $ld_command
$ld_command

annot=$((ls ${fileprefix}.${chromnumber}.annot >> /dev/null 2>&1 && echo TRUE) || echo FALSE)
if [ $annot="TRUE" ]
then
	gzip ${fileprefix}.${chromnumber}.annot
fi
