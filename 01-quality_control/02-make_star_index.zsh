#!/usr/bin/zsh

module load star/2.7.0f

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/<species>/dna/<species>.<assembly>.dna.toplevel.fa.gz
gzip -d <species>.<assembly>.dna.toplevel.fa.gz

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/<species>/<species>.<assembly>.<version>.gtf.gz
gzip -d <species>.<assembly>.<version>.gtf.gz

STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ./star_index \
--genomeFastaFiles ./<species>.<assembly>.dna.toplevel.fa \
--sjdbGTFfile ./<species>.<assembly>.<version>.gtf \
--sjdbOverhang <max(ReadLength)-1>
