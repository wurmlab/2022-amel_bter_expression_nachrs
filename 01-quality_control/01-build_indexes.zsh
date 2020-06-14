#!/bin/zsh

###############################################################################
# Script name: 01-build_indexes.zsh
# Description: build STAR and kallisto indexes using Ensembl data
# Author(s): Federico Lopez-Osorio, Alicja Witwicka
###############################################################################

# Load modules
module load star/2.7.0f

SPECIES=(Amel Bter) 
BASEDIR="/data/archive/archive-SBCS-WurmLab/awitwicka/2020_amel-bter-expression-nachrs"

# Check that working directory exists
if [ -d "$BASEDIR" ]
then
	echo "Files will be processed in ${BASEDIR}."
else
	echo "Error: ${BASEDIR} not found. Cannot continue."
	exit 1
fi

cd ${BASEDIR}
mkdir tmp
cd tmp

# Build STAR indexes
echo "Downloading Apis mellifera genome and annotation from Ensembl"
wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/apis_mellifera/dna/Apis_mellifera.Amel_4.5.dna.toplevel.fa.gz
gzip -d Apis_mellifera.Amel_4.5.dna.toplevel.fa.gz
ln -s Apis_mellifera.Amel_4.5.dna.toplevel.fa Amel.dna.fa

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/apis_mellifera/Apis_mellifera.Amel_4.5.47.gtf.gz
gzip -d Apis_mellifera.Amel_4.5.47.gtf.gz
ln -s Apis_mellifera.Amel_4.5.47.gtf Amel.gtf

echo "Downloading Bombus terrestris genome and annotation from Ensembl"
wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/dna/Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz
gzip -d Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz
ln -s Bombus_terrestris.Bter_1.0.dna.toplevel.fa Bter.dna.fa

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/bombus_terrestris/Bombus_terrestris.Bter_1.0.47.gtf.gz
gzip -d Bombus_terrestris.Bter_1.0.47.gtf.gz
ln -s Bombus_terrestris.Bter_1.0.47.gtf Bter.gtf

echo "Building STAR indexes"
for species in ${SPECIES}; do
	mkdir ${species}_star_index
	cd ${species}_star_index
	STAR --runThreadN 8 --runMode genomeGenerate \
	--genomeDir . \
	--genomeFastaFiles ../${species}.dna.fa \
	--sjdbGTFfile ../${species}.gtf \
	--sjdbOverhang 100 > ${species}_star_index.log 2>&1
	cd ..
done
echo "Finished building STAR indexes"

# Build kallisto indexes
echo "Downloading Apis mellifera transcriptome from Ensembl"
wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/apis_mellifera/cdna/Apis_mellifera.Amel_4.5.cdna.all.fa.gz
gzip -d Apis_mellifera.Amel_4.5.cdna.all.fa.gz
ln -s Apis_mellifera.Amel_4.5.cdna.all.fa Amel.cdna.fa

echo "Downloading Bombus terrestris transcriptome from Ensembl"
wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/cdna/Bombus_terrestris.Bter_1.0.cdna.all.fa.gz
gzip -d Bombus_terrestris.Bter_1.0.cdna.all.fa.gz
ln -s Bombus_terrestris.Bter_1.0.cdna.all.fa Bter.cdna.fa

# Build kallisto index
ln -s ${BASEDIR}/softw/kallisto/kallisto kallisto
echo "Building kallisto indexes"
for species in $SPECIES; do
	mkdir ${species}_kallisto_index
	cd ${species}_kallisto_index
	../kallisto index -i ${species}_kallisto.idx ../${species}.cdna.fa > \
	${species}_kallisto_index.log 2>&1
	cd ..
done
echo "Finished building kallisto indexes"
