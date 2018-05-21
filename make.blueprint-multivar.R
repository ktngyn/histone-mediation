#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
argv <- commandArgs(trailingOnly = TRUE)

if(length(argv) != 3) {
    q()
}

ldidx <- as.integer(argv[1]) # e.g., 9
univar.input.file <- argv[2] # e.g., 'qtls/mono_K27AC/9_qtl.txt.gz'
out.file <- argv[3]          # e.g., 'temp.txt.gz'

################################################################
## Convert univariate to multivariate statistics
library(dplyr)
library(readr)
library(zqtl)
library(Matrix)
library(methods)

source('util.R')

################################################################

ldfile <- 'ldblocks/fourier_ls-all.eur.bed.txt' # European LD

################################################################
temp.dir <- 'temp/' %&% out.file
dir.create(temp.dir, recursive = TRUE)

ld.tab <- read.ldfile(ldfile)

chr.input <- ld.tab[ldidx, 'chr']
lb.input <- ld.tab[ldidx, 'lb']
ub.input <- ld.tab[ldidx, 'ub']

bp.geno <- 'genotypes/_EGAZ00001235598_release_vcf_06092016_.BPWP10_13_12_15_chr' %&% chr.input
plink.bp <- subset.plink(bp.geno, chr.input, lb.input, ub.input, temp.dir)

kg.geno <- '1kg_eur/chr' %&% chr.input
plink.kg <- subset.plink(kg.geno, chr.input, lb.input, ub.input, temp.dir)

.plink <- match.plink(plink.kg, plink.bp)

if(is.null(.plink) || is.null(.plink$qtl)) {
    log.msg('Failed to reconcile two plink filesets')
    q()
}

plink.bp <- .plink$qtl

take.multivar <- function(.med, plink, univar.tab) {

    .univar <- take.matched(.med, plink$BIM, univar.tab) %>%
        mutate(qtl.se = qtl.beta / qtl.z) %>%
            select(-missing, -qtl.a1, -qtl.a2)

    .effect <- matrix(.univar$qtl.beta, ncol = 1)
    .se <- matrix(.univar$qtl.se, ncol = 1)
    .xx <- plink$BED[, .univar$x.pos, drop = FALSE]

    n <- plink.bp$FAM %>% nrow()

    vb.opt <- list(vbiter = 5000, tol = 1e-8, eigen.tol = 1e-2, gammax = 1e4,
                   pi = -1, tau = -4, do.hyper = FALSE,
                   do.rescale = FALSE, do.stdize = FALSE,
                   rate = 1e-2, decay = -1e-2, print.interv = 1000)
    
    zqtl.out <- fit.zqtl(effect = .effect, effect.se = .se, X = .xx, n = n, options = vb.opt)

    ret <- effect2tab(zqtl.out$param) %>%
        select(-x.col, -y.col) %>%
            rename(multi.beta = theta, multi.se = theta.se, multi.lodds = lodds)

    ret <- bind_cols(.univar, ret) %>% mutate(med.id = .med) %>%
        mutate(qtl.z = signif(qtl.z, 2), qtl.se = signif(qtl.se, 2), qtl.beta = signif(qtl.beta, 2))
}

univar.tab <- read_tsv(univar.input.file)

med.vars <- univar.tab %>% select(med.id) %>% unlist() %>%
    unique()

result <- lapply(med.vars, take.multivar, plink = plink.bp, univar.tab = univar.tab) %>%
    bind_rows()

write_tsv(result, out.file)
