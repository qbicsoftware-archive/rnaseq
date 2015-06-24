#!/bin/sh
#PBS -q cfc
#PBS -A qbic
#PBS -l nodes=1:ppn=10:cfc
#PBS -l walltime=40:00:00
# properties = {properties}

set -e

module load bio/fastqc/0.10
module load qbic/anaconda
module load qbic/htseq/0.6.1p2
module load qbic/tophat

{exec_job}
exit 0
