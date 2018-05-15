#!/bin/bash

logdir="log/rscript/process_blueprint-qtl"
mkdir -p $logdir

qsub -P compbio_lab -o $logdir -cwd -V -l h_vmem=64g -l h_rt=6:00:00 -b y -j y -N blueprint_cis_eqtl ./step2.sh blueprint-qtl/mono_K4ME1_log2rpm_peer_10_all_summary.txt.gz processed-qtl/mono_K4ME1
qsub -P compbio_lab -o $logdir -cwd -V -l h_vmem=96g -l h_rt=12:00:00 -b y -j y -N blueprint_cis_eqtl ./step2.sh blueprint-qtl/mono_K27AC_log2rpm_peer_10_all_summary.txt.gz processed-qtl/mono_K27AC

