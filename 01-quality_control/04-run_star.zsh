#!/usr/bin/zsh

module load parallel/20170422
module load star/2.7.0f

parallel -j 8 < star_commands.sh > star_commands.log 2>&1 &!
