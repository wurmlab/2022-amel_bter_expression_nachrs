#!/usr/bin/zsh

ln -s ~/scratch/star_alignments tmp
mkdir -p tmp/star_alignments

for file in ../input/*.fastq.gz; do
    echo "Running with ${file}"
    output=tmp/star_alignments/$(basename $file .fastq.gz)
    echo "STAR --runThreadN 2 --genomeDir ./star_index/ --readFilesCommand zcat \
    --outFileNamePrefix=${output}. --outTmpDir=${output} \
    --readFilesIn $file" >> ./star_commands.sh
done
