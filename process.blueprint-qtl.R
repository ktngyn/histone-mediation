#!/usr/bin/env Rscript

# go from:
# chr:pos_ref_alt, rsid, phenotypeID, p.value, beta, Bonferroni.p.value, FDR, alt_allele_frequency, std.error_of_beta
# 10:90164_C_G rs141504207 10:1016421:1019322 7.050e-01 0.08241 1.0 6.222e-01 0.08 0.2177
 
# to: 
#chr	rs	snp.loc	med.id	qtl.a1	qtl.a2	qtl.beta	qtl.z
#1	rs13303002	1894338	ENSG00000008128.18@1	C	A	0.05796	1.302

# and split by LD block

argv <- commandArgs(trailingOnly = TRUE)

options(stringsAsFactors=FALSE)

library(stringr)

ldfile=argv[1] # ldblocks/fourier_ls-all.eur.bed.txt
qtlfile=argv[2] # blueprint-qtl/mono_K27AC_log2rpm_peer_10_all_summary.txt.gz
outdir=argv[3] # processed-qtl/mono_K27AC

ld <- read.table(ldfile,sep="\t",header=TRUE,stringsAsFactors=FALSE)
ld$chrnum= substr(trimws(ld$chr),4,5)

qtl <- read.table(qtlfile,sep=" ",header=FALSE,stringsAsFactors=FALSE)
qtl$chr <- str_split_fixed(qtl$V1,":",2)[,1]
parts <- str_split_fixed(str_split_fixed(qtl$V1,":",2)[,2],"_",3)
qtl$pos <- as.numeric(parts[,1])
qtl$ref <- parts[,2]
qtl$alt <- parts[,3]
qtl$z <- qnorm(qtl$V4/2,lower.tail=FALSE)*sign(qtl$V5)
outqtl <- qtl[,c(10,2,11,3,12,13,5,14)]
colnames(outqtl) <- c("chr","rs","snp.loc","med.id","qtl.a1","qtl.a2","qtl.beta","qtl.z")

sapply(1:dim(ld)[1],function(x){
	filename=paste0(outdir,"/",rownames(ld)[x],"_qtl.txt.gz")
	curr=outqtl[outqtl$chr==ld[x,4]&outqtl[,3]>ld[x,2]&outqtl[,3]<=ld[x,3],]
	curr=curr[order(curr[,3]),]
	write.table(curr,gzfile(filename),col.names=TRUE,row.names=FALSE,quote=FALSE,sep="\t")
})


