#!/usr/bin/zsh

module load parallel/20170422
module load qualimap/2.2.1

parallel -j 8 'qualimap rnaseq -bam {} \
-gtf <species>.<assembly>.<version>.gtf \
-outdir results/{/.}.qualimap --sequencing-protocol non-strand-specific' ::: \
../star_alignments/*.bam > qualimap.log 2>&1 &!
