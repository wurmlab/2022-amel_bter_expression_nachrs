#!/bin/zsh

###############################################################################
# Script name: 00-create_data_symlinks.zsh
# Description: create symbolic links with reads on "db"
# Author(s): Federico Lopez-Osorio, Alicja Witwicka
###############################################################################

BASEDIR="/data/archive/archive-SBCS-WurmLab/awitwicka/2020-amel_bter_expression_nachrs"
AMELDATASETS=(
    "2018-christen_worker_brains"
    "2020-01-Jasper_2015"
    "2015-liberti_queen_brains"
    "2015-manfredini_queen_brains"
    )
BTERDATASETS=(
    "2019_colgan_queen-worker_head"
    "2015_harrison_queen_worker_whole_body"
    "2017_manfredini_queen_brain"
    "2019_porath_worker_brain"
    )

# Check that working directory exists
if [ -d "$BASEDIR" ]
then
    echo "Files will be processed in ${BASEDIR}."
else
    echo "Error: ${BASEDIR} not found. Cannot continue."
    exit 1
fi

cd ${BASEDIR}
mkdir input
cd input

echo "Creating symbolic links for honey bee data sets: ${AMELDATASETS}"
for dataset in ${AMELDATASETS}; do
    echo "Creating symbolic links for dataset ${dataset}"
    mkdir ${dataset}
    cd ${dataset}
    ln -s /data/SBCS-WurmLab/archive/db/rna/reads/Apis_mellifera/${dataset}/renamed/*.fastq.gz .
    cd ..
done

echo "Creating symbolic links for bumble bee data sets: ${BTERDATASETS}"
for dataset in ${BTERDATASETS}; do
    echo "Creating symbolic links for dataset ${dataset}"
    mkdir ${dataset}
    cd ${dataset}
    ln -s /data/SBCS-WurmLab/archive/db/rna/reads/Bombus_terrestris/${dataset}/*.fastq.gz .
    cd ..
done
