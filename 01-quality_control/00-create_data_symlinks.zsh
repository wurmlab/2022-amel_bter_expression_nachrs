#!/bin/zsh

###############################################################################
# Script name: 00-create_data_symlinks.zsh
# Description: create symbolic links with reads on "db"
# Author(s): Federico Lopez-Osorio, Alicja Witwicka
###############################################################################

# Run script within "input" directory

AMELDATASETS=(
    "2015-jasper_worker_tissues"
    "2015-liberti_queen_brains"
    "2015-manfredini_queen_brains"
    "2018-christen_worker_brains"
    )
BTERDATASETS=(
    "2019-colgan_queen_worker_head"
    "2015-harrison_queen_worker_whole_body"
    "2017-manfredini_queen_brain"
    "2019-porath_worker_brain"
    )

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
