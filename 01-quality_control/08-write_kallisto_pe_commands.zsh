#!/usr/bin/zsh

# kallisto-0.46.2

kallisto_ref=kallisto_index/<assembly>_cdna.idx
ln -s ~/scratch/kallisto_alignments tmp
mkdir -p tmp/kallisto_alignments

for file in ../input/*_1.fastq.gz; do
    identifier=$(basename $file _1.fastq.gz)
    echo "Running with ${identifier}"
    output=tmp/kallisto_alignments/$identifier
    echo "kallisto quant -i ${kallisto_ref} --bias --bootstrap-sample=100 \
    --output-dir=${output} --threads=2 \
    ../input/$identifier"_1.fastq.gz" \../input/$identifier"_2.fastq.gz"" >> \
    ./kallisto_commands.sh
done
