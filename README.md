# Gene expression of nAChR subunit genes in _Apis mellifera_ and _Bombus terrestris_
Workflow and code for the comparative analysis of gene expression of nAChRs subunit genes in *Apis mellifera* and *Bombus terrestris*.

**Disclaimer**: The materials in this repository enable the reproducibility of our results. Except where noted, the code herein is specific to this project and to the Apocrita HPC Cluster at QMUL.

* Federico Lopez-Osorio (f.lopez-osorio at qmul.ac.uk)
* Alicja Witwicka (a.witwicka at qmul.ac.uk)
* Yannick Wurm (y.wurm at qmul.ac.uk)

## Case studies
We included data from the following case studies:
| Study                    | Species             | Caste        | Tissue     | ENA                                                          | DOI                        |
| ------------------------ | ------------------- | ------------ | ---------- | ------------------------------------------------------------ | -------------------------- |
| Christen *et al.* 2018   | *Apis mellifera*    | Worker       | Brain      | [PRJNA450296](https://www.ebi.ac.uk/ena/browser/view/PRJNA450296) | 10.1021/acs.est.8b01801    |
| Manfredini *et al.* 2015 | *Apis mellifera*    | Queen        | Brain      | [PRJNA275154](https://www.ebi.ac.uk/ena/browser/view/PRJNA275154) | 10.1186/s12864-015-1750-7  |
| Liberti *et al.* 2015    | *Apis mellifera*    | Queen        | Brain      | [PRJNA524311](https://www.ebi.ac.uk/ena/browser/view/PRJNA524311) | 10.7554/eLife.45009        |
| Jasper *et al.* 2015     | *Apis mellifera*    | Worker       | Multiple   | [PRJNA243651](https://www.ebi.ac.uk/ena/browser/view/PRJNA243651) | 10.1093/molbev/msu292      |
| Colgan *et al*. 2019     | *Bombus terrestris* | Queen Worker | Head       | [PRJNA508397](https://www.ebi.ac.uk/ena/browser/view/PRJNA508397) | 10.1111/mec.15047          |
| Harrison *et al*. 2015   | *Bombus terrestris* | Queen Worker | Whole body | [PRJEB9366](https://www.ebi.ac.uk/ena/browser/view/PRJEB9366) | 10.1111/mec.13215          |

## Load modules
```
module load fastqc/0.11.9
module load parallel/20170422
module load qualimap/2.2.1
module load samtools/1.9
module load star/2.7.0f
```

## 01-quality_control
FastQC
```
parallel -j 10 'fastqc {}' ::: ../input/*.fastq.gz > YYYY-MM-DD-fastqc.log 2>&1 &!
```

STAR
```
# Single-end
STAR --runThreadN 2 --genomeDir ./star_index/ --readFilesCommand zcat \
--outFileNamePrefix=sample_id. --outTmpDir=outdir \
--readFilesIn ../input/sample_id.fastq.gz > star.log 2>&1 &!

# Paired-end
STAR --runThreadN 2 --genomeDir ./star_index/ --readFilesCommand zcat \
--outFileNamePrefix=sample_id. --outTmpDir=outdir \
--readFilesIn ../input/sample_id_1.fastq.gz ../input/sample_id_2.fastq.gz > star.log 2>&1 &!
```

Convert SAM to BAM
```
parallel -j 10 'samtools view -@ 4 -bS {} | \
samtools sort -@ 4 - -o {.}.sorted.bam' ::: *.sam > samtools.log 2>&1 &!
```

Qualimap
```
# Non-strand-specific
parallel -j 10 'qualimap rnaseq -bam {} \
-gtf <species>.<assembly>.<version>.gtf.gz \ # Ensembl Metazoa 46
-outdir results/{/.}.qualimap --sequencing-protocol non-strand-specific ::: \
../YYYY-MM-DD-star/star_alignments/*.bam > qualimap.log 2>&1 &!
```

## 02-kallisto_quantification
kallisto
```
# Single-end and non-stranded
kallisto quant -i <version>_cDNA.idx --bias --bootstrap-sample=100 \
--output-dir=sample_id --single --fragment-length=300 --sd=20 \
--threads=2 sample_id.fastq.gz" >

# Single-end and stranded (e.g., strands-pecific reads, first read reverse)
kallisto quant -i <version>_cDNA.idx --bias --bootstrap-sample=100 \
--output-dir=sample_id --single --fragment-length=300 --sd=20 \
--threads=2 --rf-stranded sample_id.fastq.gz

# Paired-end and non-stranded
kallisto quant -i ${kallisto_ref} --bias --bootstrap-sample=100 \
--output-dir=sample_id --threads=2 \
sample_id_1.fastq.gz sample_id_2.fastq.gz

```

## 03-summary_report
MultiQC
```
# Generate a single MultiQC report including results from all analyses above
multiqc .
```
