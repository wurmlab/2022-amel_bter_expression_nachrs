#!/usr/bin/zsh

ln -s ~/scratch/star_alignments tmp
mkdir -p tmp/star_alignments

for file in ../input/*_1.fastq.gz; do
    identifier=$(basename $file _1.fastq.gz)
    echo "Running with ${identifier}"
    output=tmp/star_alignments/$identifier
    echo "STAR --runThreadN 2 --genomeDir ./star_index/ --readFilesCommand zcat \
    --outFileNamePrefix=${output}. --outTmpDir=${output} \
    --readFilesIn ../input/$identifier"_1.fastq.gz" ../input/$identifier"_2.fastq.gz" " >> ./star_commands.sh
done
