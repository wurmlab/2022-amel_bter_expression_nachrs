#!/usr/bin/zsh

module load samtools/1.9

parallel -j 8 'samtools view -@ 4 -bS {} | \
samtools sort -@ 4 - -o {.}.sorted.bam' ::: *.sam > samtools.log 2>&1 &!
