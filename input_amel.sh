# I got cDNA A. mellifera genome from ENSEMBL Metazoa:                                                                                                                                                   
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-46/fasta/apis_mellifera/cdna/Apis_mellifera.Amel_4.5.cdna.all.fa.gz

# I download data for different tissues for foragers and queens from Jasper et al. 2015                                                                                                                  
# I create a temporary (tmp) directory in scratch                                                                                                                                                        
ln -s ~/scratch tmp
# And copy there all the previously download files from                                                                                                                                                  
# https://www.ebi.ac.uk/ena/data/view/PRJNA243651                                                                                                                                                        

# All the data used in this project is now in:                                                                                                                                                           
~/db/rna/reads/Apis_mellifera
# That is: Manfredini et al. 2015, Liberti et al. 2015 (queens), Christen et al. 2018, Jasper et al. 2015 (workers).                                                                                     
                                                                                                                                                                                
##FastQC:                                                                                                                                                                                                

## required:                                                                                                                                                                                             
module load fastqc/0.11.9

### run quantifications                                                                                                                                                                                  
for name in *.fastq.gz;
do
fastqc "$name" -t 40;
done

##STAR:

##required:
module load star/2.7.0f

#Create index:
mkdir star_amel_index

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir star_amel_index \
--genomeFastaFiles ../input/Apis_mellifera.Amel_4.5.dna.toplevel.fa
--sjdbGTFfile ../input/Apis_mellifera.Amel_4.5.47.gtf
--sjdbOverhang 100

 #paired end:
for FILE in ./*.1.fastq.gz
 do
    SRR_ID=$(basename $FILE .1.fastq.gz)
    STAR --genomeDir ./star_amel_index \
     --runThreadN 14 \
     --readFilesIn ./${SRR_ID}".1.fastq.gz" ./${SRR_ID}".2.fastq.gz" \
     --outFileNamePrefix ./STAR_UTPUT/$SRR_ID \
     --readFilesCommand gunzip -c \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard > ${SRR_ID}_star.log 
 done &!

 #single end:
for FILE in ./*.fastq.gz
 do
    SRR_ID=$(basename $FILE .fastq.gz)
    STAR --genomeDir ./star_amel_index \
     --runThreadN 8 \
     --readFilesIn ./${SRR_ID}".fastq.gz" \
     --outFileNamePrefix ./STAR_UTPUT/$SRR_ID \
     --readFilesCommand gunzip -c \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMunmapped Within \
     --outSAMattributes Standard > ${SRR_ID}_star.log 
 done &!


##QualiMap:

##required: 
module load qualimap/2.2.1

#I quantify with salmon to assess if data are strand-specific or not for a single sample from each study:# I ran it only for a single sample for each study included in the analyses.
salmon quant -t ./Apis_mellifera.Amel_4.5.dna.toplevel.fa -l A -a sample_1Aligned.sortedByCoord.out.bam -o salmon_quant 
#all resulted as unstranded 

#paired-end:
for FILE in ../STAR/STAR_OUTPUT/jasper_2015/*.sortedByCoord.out.bam
do
    BAM_ID=$(basename $FILE .sortedByCoord.out.bam)
    qualimap rnaseq -bam $FILE -gtf Apis_mellifera.Amel_4.5.47.gtf -pe -outdir ../tm/qualimap/jasper_2015/$BAM_ID > ${BAM_ID}_star.log 
done

#single-end:
for FILE in ../STAR/STAR_SAM_OUTPUT/christen_2018/*.sortedByCoord.out.bam
do
    BAM_ID=$(basename $FILE .sortedByCoord.out.bam)
    qualimap rnaseq -bam $FILE -gtf Apis_mellifera.Amel_4.5.47.gtf -outdir ../tm/qualimap/christen_2018/$BAM_ID > ${BAM_ID}_star.log 
done

##Kallisto

# Define reference file for Kallisto
kallisto_ref=tm/Amel_index

## Generate kallisto-index:
./kallisto index -i $kallisto_ref input/Apis_mellifera.Amel_4.5.cdna.all.fa.gz 

### Create a directory for the results 
mkdir tmp/kallisto_output

#single-end:
for FILE in *.fastq.gz
 do
    OUTPUT_ID=$(basename $FILE .fastq.gz)
    ./kallisto quant --index=tm/Amel_index --bias --bootstrap-sample=100 \
    --output-dir=kallisto_output/$OUTPUT_ID \
    --threads=10 --single -l 300 -s 20  input/renamed/${OUTPUT_ID}.fastq.gz > ${OUTPUT_ID}__kallisto_run.log 
 done

 #paired-end:
 for FILE in *.1.fastq.gz
 do
     OUTPUT_ID=$(basename $FILE .1.fastq.gz)
    ./kallisto quant --index=tm/Amel_index --bias --bootstrap-sample=100 \
    --output-dir=kallisto_output/$OUTPUT_ID \
    --threads=10 input/renamed/${OUTPUT_ID}".1.fastq.gz" input/renamed/${OUTPUT_ID}".2.fastq.gz" > ${OUTPUT_ID}_kallisto_run.log 
 done