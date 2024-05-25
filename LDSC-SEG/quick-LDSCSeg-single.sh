#!/bin/bash
while getopts s:l:o:n: option
do
	case "${option}"
	in
	s) sumstats=${OPTARG};;
	l) ldcts_file=${OPTARG};;
	o) output_path=${OPTARG};;
	n) name=${OPTARG};;
	esac
done

prefix0=$(ls $sumstats | tr '/' '\n' | tail -n1  )
prefix1=${prefix0%".sumstats.gz"}
prefix2=${prefix1#"PASS_"}
output=$output_path/$prefix2-$name
echo "Running: $prefix2-$name"

cd ldsc
source /programs/biogrids.shrc
module load anaconda2

python ldsc.py \
        --h2-cts "$sumstats" \
        --ref-ld-chr baseline_v1.2/baseline. \
        --out  $output \
        --ref-ld-chr-cts $ldcts_file \
        --w-ld-chr weights_hm3_no_hla/weights. \
        --print-coefficients
