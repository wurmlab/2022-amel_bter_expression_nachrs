#!/usr/bin/zsh

# kallisto-0.46.2

kallisto_ref=kallisto_index/<assembly>_cdna.idx
ln -s ~/scratch/kallisto_alignments tmp
mkdir -p tmp/kallisto_alignments

for file in ../input/*.fastq.gz; do
    echo "Running with ${file}"
    output=tmp/kallisto_alignments/$(basename $file .fastq.gz)
    echo "kallisto quant -i ${kallisto_ref} --bias --bootstrap-sample=100 \
    --output-dir=${output} --single --fragment-length=300 --sd=20 \
    --threads=2 --rf-stranded ${file}" >> ./kallisto_commands.sh
done
