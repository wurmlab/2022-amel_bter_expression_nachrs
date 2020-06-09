#!/usr/bin/zsh

module load parallel/20170422
# kallisto-0.46.2

parallel -j 8 < kallisto_commands.sh > kallisto_commands.log 2>&1 &!
