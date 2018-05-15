#!/bin/bash

source /broad/tools/scripts/useuse
reuse R-3.4

LD="ldblocks/fourier_ls-all.eur.bed.txt"
QTLfile=$1
outdir=$2

mkdir -p $outdir
Rscript process.blueprint-qtl.R $LD $QTLfile $outdir

