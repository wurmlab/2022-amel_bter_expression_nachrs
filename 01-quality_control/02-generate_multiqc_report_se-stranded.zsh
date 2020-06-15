#!/bin/zsh

###############################################################################
# Script name: 02-generate_multiqc_report_se-stranded.zsh
# Description: use MultiQC to aggregate FastQC, STAR, and Qualimap results 
# generated with single-end, strand-specific RNA-seq data
# Author(s): Federico Lopez-Osorio, Alicja Witwicka
###############################################################################

# Load modules
module load fastqc/0.11.9
module load parallel/20170422
module load star/2.7.0f
module load samtools/1.9

SPECIES="Bter"
DATASET="2019_colgan_queen-worker_head"
DATE="2020-06-10"
BASEDIR="/data/archive/archive-SBCS-WurmLab/awitwicka/2020_amel-bter-expression-nachrs"
INPUT="/data/archive/archive-SBCS-WurmLab/awitwicka/2020_amel-bter-expression-nachrs/input"
SCRATCH="/data/scratch/btx422" # Set to user scratch

# Check that working directory exists
if [ -d "$BASEDIR" ]
then
    echo "Files will be processed in ${BASEDIR}."
else
    echo "Error: ${BASEDIR} not found. Cannot continue."
    exit 1
fi

cd ${BASEDIR}
echo 'Creating "results" directory in '${BASEDIR}
mkdir -p results
cd results

# Run FastQC, STAR, and Qualimap
mkdir ${DATASET}
cd ${DATASET}
ln -s ${INPUT}/${DATASET} input

# Run FastQC
echo "Running FastQC with data set "${DATASET}
mkdir ${DATE}-fastqc
cd ${DATE}-fastqc
parallel -j 8 'fastqc {}' ::: ../input/*.fastq.gz > fastqc.log 2>&1
mv ../input/*.{zip,html} .
echo "Finished running FastQC"
cd ..

# Run STAR
echo "Running STAR with data set "${DATASET}
mkdir ${DATE}-star
cd ${DATE}-star
echo "Creating symbolic link for bumble bee index"
ln -s ${BASEDIR}/tmp/Bter_star_index star_index

mkdir ${SCRATCH}/${DATE}
ln -s ${SCRATCH}/${DATE} tmp
mkdir -p tmp/star
for file in ../input/*.fastq.gz; do
    echo "Writing STAR commands for ${file}"
    output=tmp/star/$(basename $file .fastq.gz)
    echo "STAR --runThreadN 2 --genomeDir ./star_index/ --readFilesCommand zcat \
    --outFileNamePrefix=${output}. --outTmpDir=${output} \
    --readFilesIn $file" >> ./star_commands.sh
done
parallel -j 8 < star_commands.sh > star_commands.log 2>&1
echo "Finished running STAR commands"

# Convert SAM files to sorted BAM
echo "Creating BAM files with data set "${DATASET}
parallel -j 8 'samtools view -@ 4 -bS {} | \
samtools sort -@ 4 - -o {.}.sorted.bam' ::: tmp/star/*.sam > samtools.log 2>&1
echo "Finished creating BAM files"
# Remove SAM files
rm tmp/star/*.sam
cd ..

# Run Qualimap
echo "Running Qualimap with data set "${dataset}
module unload fastqc/0.11.9
module unload java
module load qualimap/2.2.1

mkdir ${DATE}-qualimap
cd ${DATE}-qualimap
ln -s ${SCRATCH}/${DATE}/star star
GTF=${BASEDIR}/tmp/Bombus_terrestris.Bter_1.0.47.gtf
export GTF

parallel -j 8 'qualimap rnaseq -bam {} -gtf "$GTF" \
-outdir results/{/.}.qualimap \
--sequencing-protocol strand-specific-reverse' ::: \
star/*.bam > qualimap.log 2>&1
echo "Finished running Qualimap"
cd ..

# Run kallisto
echo "Running kallisto with data set "${DATASET}
mkdir ${DATE}-kallisto
cd ${DATE}-kallisto
ln -s ${BASEDIR}/softw/kallisto/kallisto kallisto
ln -s ${SCRATCH}/${DATE} tmp
mkdir -p tmp/kallisto
echo "Creating symbolic link for bumble bee index"
ln -s ${BASEDIR}/tmp/Bter_kallisto_index/Bter_kallisto.idx kallisto_index
colonies_castes=(C61_queen C48_queen C34_queen C06_queen C06_worker C61_worker C48_worker C34_worker)

for sample in $colonies_castes; do
    echo "Running with ${sample}"
    fastqs=$(find ../input/ -name "*${sample}*.fastq.gz" | tr "\n" " ")
    output=tmp/kallisto/${sample}"_head"
    echo "./kallisto quant -i kallisto_index --bias --bootstrap-sample=100 \
    --output-dir=${output} --single --fragment-length=300 --sd=20 \
    --rf-stranded --threads=2 ${fastqs}" >> ./kallisto_commands.sh
done

parallel -j 8 < kallisto_commands.sh > kallisto.log 2>&1
echo "Finished running kallisto commands"
mkdir results
mv tmp/kallisto/* results/
rm -r tmp
cd ..

# Delete BAM files in STAR results
rm ${BASEDIR}/results/${DATASET}/${DATE}-star/tmp/star/*.bam
# Move other STAR output files to STAR results directory
cd ${BASEDIR}/results/${DATASET}/${DATE}-star
mkdir results
mv tmp/star/* results/
# Remove temporary directories
rm -r tmp ${SCRATCH}/${DATE}

# Generate MultiQC report
cd ${BASEDIR}/results/${DATASET}
multiqc .
# Rename MultiQC report and data
mv multiqc_report.html ${DATASET}-multiqc_report.html
mv multiqc_data ${DATASET}-multiqc_data
