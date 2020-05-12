#!/usr/bin/zsh

module load fastqc/0.11.9
module load parallel/20170422

parallel -j 8 'fastqc {}' ::: ../input/*.fastq.gz > fastqc.log 2>&1 &!
