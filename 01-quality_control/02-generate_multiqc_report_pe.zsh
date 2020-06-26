#!/bin/zsh

###############################################################################
# Script name: 02-generate_multiqc_report_pe.zsh
# Description: use MultiQC to aggregate FastQC, STAR, and Qualimap results 
# generated with paired-end, unstranded RNA-seq data
# Author(s): Federico Lopez-Osorio, Alicja Witwicka
###############################################################################

# Run script within "results" directory

# Load modules
module load fastqc/0.11.9
module load parallel/20170422
module load star/2.7.0f
module load samtools/1.9

SPECIES="Bter" # Set either "Amel" or "Bter"
DATE="2020-06-16"
SCRATCH="/data/scratch/btx422" # Set to user scratch

# Select data sets for analysis
if [ ${SPECIES} = "Amel" ]
then
    DATASETS=(
        "2015-jasper_worker_tissues"
        )
    echo "Processing honey bee data sets: " "${DATASETS[@]}"
else
    DATASETS=(
        "2017-manfredini_queen_brain"
        "2019-porath_worker_brain"
        )
    echo "Processing bumble bee data sets: " "${DATASETS[@]}"
fi

# Create symbolic links for "input", "softw", and "tmp" directories
ln -sfn ../input .
ln -sfn ../softw .
ln -sfn ../tmp .

# Run FastQC, STAR, and Qualimap
for dataset in "${DATASETS[@]}"; do
    mkdir "${dataset}"
    cd "${dataset}" || exit

    # Run FastQC
    echo "Running FastQC with data set ""${dataset}"
    mkdir ${DATE}-fastqc
    cd ${DATE}-fastqc || exit
    parallel -j 8 'fastqc {}' ::: ../../input/"${dataset}"/*.fastq.gz > fastqc.log 2>&1
    mv ../../input/"${dataset}"/*.{zip,html} .
    echo "Finished running FastQC"
    cd ..

    # Run STAR
    echo "Running STAR with data set ""${dataset}"
    mkdir ${DATE}-star
    cd ${DATE}-star || exit
    
    if [ ${SPECIES} = "Amel" ]
    then
        echo "Creating symbolic link for honey bee index"
        ln -s ../../tmp/Amel_star_index star_index
    else
        echo "Creating symbolic link for bumble bee index"
        ln -s ../../tmp/Bter_star_index star_index
    fi

    mkdir ${SCRATCH}/${DATE}
    ln -s ${SCRATCH}/${DATE} tmp
    mkdir -p tmp/star
    for file in ../../input/"${dataset}"/*_1.fastq.gz; do
        identifier=$(basename "$file" _1.fastq.gz)
        echo "Writing STAR commands for sample ${identifier}"
        output=tmp/star/$(basename "$file" _1.fastq.gz)
        echo "STAR --runThreadN 2 --genomeDir ./star_index/ --readFilesCommand zcat \
        --outFileNamePrefix=${output}. --outTmpDir=${output} \
        --readFilesIn ../../input/${dataset}/${identifier}"_1.fastq.gz" ../../input/${dataset}/${identifier}"_2.fastq.gz"" >> \
        ./star_commands.sh
    done
    parallel -j 8 < star_commands.sh > star_commands.log 2>&1
    echo "Finished running STAR commands"

    # Convert SAM files to sorted BAM
    echo "Creating BAM files with data set ""${dataset}"
    parallel -j 8 'samtools view -@ 4 -bS {} | \
    samtools sort -@ 4 - -o {.}.sorted.bam' ::: tmp/star/*.sam > samtools.log 2>&1
    echo "Finished creating BAM files"
    rm tmp/star/*.sam
    cd ..

    # Run Qualimap
    echo "Running Qualimap with data set ""${dataset}"
    module unload fastqc/0.11.9
    module unload java
    module load qualimap/2.2.1

    mkdir ${DATE}-qualimap
    cd ${DATE}-qualimap || exit
    ln -s ${SCRATCH}/${DATE}/star star
    if [ ${SPECIES} = "Amel" ]
    then
        GTF="../../tmp/Apis_mellifera.Amel_4.5.47.gtf"
        export GTF
    else
        GTF="../../tmp/Bombus_terrestris.Bter_1.0.47.gtf"
        export GTF
    fi
    parallel -j 8 'qualimap rnaseq -bam {} -gtf "$GTF" --paired \
    -outdir results/{/.}.qualimap \
    --sequencing-protocol non-strand-specific' ::: \
    star/*.bam > qualimap.log 2>&1
    echo "Finished running Qualimap"
    rm star
    cd ..

    # Run kallisto
    echo "Running kallisto with data set ""${dataset}"
    mkdir ${DATE}-kallisto
    cd ${DATE}-kallisto || exit
    ln -s ../../softw/kallisto/kallisto kallisto
    ln -s ${SCRATCH}/${DATE} tmp
    mkdir -p tmp/kallisto

    if [ ${SPECIES} = "Amel" ]
    then
        echo "Creating symbolic link for honey bee index"
        ln -s ../../tmp/Amel_kallisto_index/Amel_kallisto.idx kallisto_index
    else
        echo "Creating symbolic link for bumble bee index"
        ln -s ../../tmp/Bter_kallisto_index/Bter_kallisto.idx kallisto_index
    fi

    for file in ../../input/"${dataset}"/*_1.fastq.gz; do
        identifier=$(basename "$file" _1.fastq.gz)
        echo "Writing kallisto commands for sample ${identifier}"
        output=tmp/kallisto/${identifier}
        echo "./kallisto quant -i kallisto_index --bias --bootstrap-sample=100 \
        --output-dir=${output} --threads=2 \
        ../../input/${dataset}/${identifier}"_1.fastq.gz" ../../input/${dataset}/${identifier}"_2.fastq.gz"" >> \
        ./kallisto_commands.sh
    done
    parallel -j 8 < kallisto_commands.sh > kallisto.log 2>&1
    echo "Finished running kallisto commands"
    mkdir results
    mv tmp/kallisto/* results/
    rm -r tmp
    cd ..

    # Delete BAM files in STAR results
    rm ${DATE}-star/tmp/star/*.bam
    # Move other STAR output files to STAR results directory
    cd ${DATE}-star || exit
    mkdir results
    mv tmp/star/* results/
    # Remove temporary directories
    rm -r tmp ${SCRATCH:?}/${DATE}
    cd ..

    # Generate MultiQC report
    multiqc .
    # Rename MultiQC report and data
    mv multiqc_report.html "${dataset}"-multiqc_report.html
    mv multiqc_data "${dataset}"-multiqc_data
    cd ..

    module unload qualimap/2.2.1
    module unload java
    module load fastqc/0.11.9
done
