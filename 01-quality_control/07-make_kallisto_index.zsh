#!/usr/bin/zsh

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/<species>/<species>.<assembly>.cdna.all.fa.gz
gzip -d <species>.<assembly>.cdna.all.fa.gz

kallisto index -i <assembly>_cdna.idx <species>.<assembly>.cdna.all.fa > kallisto_index.log 2>&1 &!
