options(stringsAsFactors = FALSE)

log.msg <- function(...) {
    ss <- as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
}

`%&%` <- function(a, b) paste(a, b, sep = '')

read.ldfile <- function(ldfile) {
    require(dplyr)
    ld.tab <- read.table(ldfile, sep = '\t', header = TRUE) %>%
        dplyr::rename(lb = `start`, ub = `stop`) %>%
            dplyr::mutate(chr = gsub(chr, pattern = 'chr', replacement = '')) %>%
                dplyr::mutate(chr = as.integer(chr), lb = as.integer(lb), ub = as.integer(ub))
}

subset.plink <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {
    require(zqtl)
    require(dplyr)

    .error <- function(e) {
        log.msg('No QTL here!\n')
        return(NULL)
    }

    .subset <- function(plink.hdr, chr, plink.lb, plink.ub, temp.dir) {

        chr.num <- gsub(pattern = 'chr', replacement = '', chr) %>% as.integer()
        plink.cmd <- sprintf('./bin/plink --bfile %s --make-bed --geno 0.05 --maf 0.05 --chr %d --from-bp %d --to-bp %d --out %s', plink.hdr, chr.num, plink.lb, plink.ub, (temp.dir %&% '/plink'))
        system(plink.cmd)

        plink <- read.plink((temp.dir %&% '/plink'))
        colnames(plink$BIM) <- c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
        colnames(plink$FAM) <- c('fam', 'iid', 'father', 'mother', 'sex.code', '.pheno')
        plink$FAM <- plink$FAM %>% mutate(iid = sapply(iid, gsub, pattern = 'GTEX-', replacement = ''))

        if(any(is.logical(plink$BIM$plink.a1))) {
            plink$BIM$plink.a1 <- 'T'
        }

        if(any(is.logical(plink$BIM$plink.a2))) {
            plink$BIM$plink.a2 <- 'T'
        }
        return(plink)
    }

    plink <- tryCatch(.subset(plink.hdr, chr, plink.lb, plink.ub, temp.dir),
                      error = .error)
    return(plink)
}

take.matched <- function(.med, x.bim, univar.tab) {
    require(dplyr)

    .univar <- univar.tab %>% dplyr::filter(med.id == .med)

    if('x.pos' %in% colnames(x.bim)) {
        ret <- x.bim %>% left_join(.univar)
    } else {
        ret <- x.bim %>% mutate(x.pos = 1:n()) %>% left_join(.univar)
    }

    ret.matched <- ret %>% na.omit() %>%
        dplyr::filter((plink.a1 == qtl.a1 & plink.a2 == qtl.a2) |
                   (plink.a1 == qtl.a2 & plink.a2 == qtl.a1))

    ret.matched <- ret.matched %>%
        mutate(flip = if_else(qtl.a1 == plink.a1, 1.0, -1.0)) %>%
            mutate(qtl.beta = flip * qtl.beta, qtl.z = flip * qtl.z) %>%
                dplyr::select(-flip)

    ret <- ret.matched %>%
        arrange(snp.loc)
}

mat2tab <- function(mat, val.name) {
    require(dplyr)
    require(tidyr)
    ret <- as.matrix(mat) %>% as.data.frame()
    col.names <- 1:ncol(ret)
    colnames(ret) <- col.names
    ret <- ret %>%
        mutate(x.col = 1:n()) %>%
            gather_(key_col= 'y.col', value_col = val.name, col.names)
    return(ret)
}

effect2tab <- function(param) {
    require(dplyr)
    require(tidyr)
    .theta <- mat2tab(param$theta, val.name = 'theta')
    .theta.var <- mat2tab(param$theta.var, val.name = 'theta.var')
    .lodds <- mat2tab(param$lodds, val.name = 'lodds')
    ret <- .theta %>%
        left_join(.theta.var) %>%
            left_join(.lodds) %>%
                mutate(theta.se = sqrt(theta.var)) %>%
                    select(-theta.var) %>%
                        mutate(x.col = as.integer(x.col)) %>%
                            mutate(y.col = as.integer(y.col))
    ret <- ret %>%
        mutate(theta = signif(theta, 2),
               theta.se = signif(theta.se, 2),
               lodds = signif(lodds, 2))
}

match.plink <- function(plink.gwas, plink.qtl) {

    if(is.null(plink.gwas)) return(NULL)
    if(is.null(plink.qtl)) return(NULL)

    ret.gwas <- plink.gwas
    ret.qtl <- plink.qtl

    gwas.bim <- plink.gwas$BIM %>%
        mutate(gwas.x.pos = 1:n()) %>%
            rename(gwas.plink.a1 = plink.a1,
                   gwas.plink.a2 = plink.a2) %>%
                       select(-missing)

    qtl.bim <- plink.qtl$BIM %>%
        mutate(qtl.x.pos = 1:n()) %>%
            rename(qtl.plink.a1 = plink.a1,
                   qtl.plink.a2 = plink.a2,
                   qtl.rs = rs) %>%
                       select(-missing)

    bim.matched <- gwas.bim %>%
        left_join(qtl.bim) %>%
            na.omit()

    if(nrow(bim.matched) < 1) return(NULL)

    bim.matched <- bim.matched %>%
        dplyr::filter(((gwas.plink.a1 == qtl.plink.a1) & (gwas.plink.a2 == qtl.plink.a2)) |
                          ((gwas.plink.a2 == qtl.plink.a1) & (gwas.plink.a1 == qtl.plink.a2))) %>%
                              arrange(chr, snp.loc)

    if(nrow(bim.matched) < 1) return(NULL)

    ret.gwas$BIM <- ret.gwas$BIM[bim.matched$gwas.x.pos, , drop = FALSE]
    ret.gwas$BED <- ret.gwas$BED[ , bim.matched$gwas.x.pos, drop = FALSE]

    ret.qtl$BIM <- ret.qtl$BIM[bim.matched$qtl.x.pos, , drop = FALSE]
    ret.qtl$BED <- ret.qtl$BED[ , bim.matched$qtl.x.pos, drop = FALSE]

    flip.tab <- ret.gwas$BIM %>% mutate(gwas.x.pos = 1:n()) %>%
        left_join(ret.qtl$BIM %>% mutate(qtl.x.pos = 1:n()),
                  by = c('chr', 'snp.loc'),
                  suffix = c('.gwas', '.qtl')) %>%                      
                      filter(plink.a1.gwas != plink.a1.qtl)

    ret.qtl$BIM[flip.tab$qtl.x.pos, ] <- ret.gwas$BIM[flip.tab$gwas.x.pos, ]

    flip.bed <- ret.qtl$BED[, flip.tab$qtl.x.pos]
    zero.idx <- flip.bed <= 0.5
    two.idx <- flip.bed >= 1.5
    flip.bed[two.idx] <- 0
    flip.bed[zero.idx] <- 2
    ret.qtl$BED[, flip.tab$qtl.x.pos] <- flip.bed

    return(list(gwas = ret.gwas, qtl = ret.qtl))
}
