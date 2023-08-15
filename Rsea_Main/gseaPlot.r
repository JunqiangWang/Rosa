# Functions for gsea analysis
# TODO
# Add ledge

#' Gene Set Enrichment Analysis
#'
#' This function performs gene eset enrichment analysis and plots
#'
#' @param signature Numeric vector containing the gene expression signature
#' @param geneset Vector of character strings indicating the genes in the geneset, or named numeric vector for 2-tail analysis
#' @param score Number indicating the exponent score for GSEA
#' @param twoTails Logical, whether the ES should be computed as the different between each tail ES
#' @param pout Logical, whether a plot should be generated
#' @param per Integer indicating the number of permutations for p-value estimation
#' @param alternative Character string indicating the tail to test for statistical significance
#' @param colSig Vector indicating the colors for the signature, including hue for negative values, hue for positive values, mate value and gama
#' @param colHit Vector indicating the colors for the hits, including hue for negative values, hue for positive values, mate value and gama
#' @param ylim Optional numeric vector indicating the limits for the y-axis
#' @param axes Logical, whether axis should be ploted
#' @param xlab Character string indicating the label for the x-axis
#' @param ylab Character string indicating the label for the y-axis
#' @param lwd Number indicating the line width for the enrichment plot
#' @param maxhit Integer indicating the maximum number of hits to plot, 0 means all
#' @param ... Additional parameters to pass to plot function
#' @return List of results and plot
#' @export
gsea <- function(signature, geneset, score=1, twoTails=FALSE, pout=TRUE, per=0, alternative=c("two.sided", "greater", "less"), colSig=c(.45, .15, .3, 1), colHit=c(.58, .05, .1, 2), ylim=NULL, axes=TRUE, xlab="Signature", ylab="ES", main="", lwd=1, maxhit=0, ...) {
    alternative <- match.arg(alternative)
    if (is.null(names(geneset))) {
        tmp <- geneset
        geneset <- rep(1, length(geneset))
        names(geneset) <- tmp
    }
    if (!is.null(ncol(signature))) signature <- signature[, 1]
    if (prod(range(signature))>=0) { #One tail enrichment
        twoTails <- FALSE
        signature <- sort(abs(signature), decreasing=TRUE)
        es1 <- gsea.es(signature, which(names(signature) %in% names(geneset)), score)
        pos <- which.max(abs(es1))
        es <- es1[pos]
        # LEDGE
        if (es<0) ledgesig <- signature[pos:length(signature)]
        else ledgesig <- rev(signature[1:pos])
        ledge <- ledgesig[names(ledgesig) %in% names(geneset)]
        nes <- list(nes=NULL, p.value=NULL)
        if (per>0) {
            nesnull <- sapply(1:per, function(i, signature, setsize, score) {
                es1 <- gsea.es(signature, sample(length(signature), setsize), score)
                es1[which.max(abs(es1))]
            }, signature=signature, setsize=length(which(names(signature) %in% names(geneset))), score=score)
            nes <- aecdf(nesnull, symmetric=TRUE)(es, alternative)
        }
        if (pout) {
            if (is.null(ylim)) ylim <- c(min(0, min(es1)), max(0, max(es1)))
            if (length(ylim)==1) ylim <- c(-abs(ylim), abs(ylim))
            wd <- diff(ylim)
            base.bar <- ylim[2]+wd*.05
            top.bar <- ylim[2]+wd*.15
            base.sig <- ylim[1]-wd*.15
            top.sig <- ylim[1]-wd*.05
            plot(0, 0, type="n", xlim=c(0, length(signature)+1), ylim=c(base.sig, top.bar), axes=axes, xlab=xlab, ylab=ylab, main=main, font.main=1, ...)
            abline(h=0, col="grey")
            pos <- which.max(abs(es1))
            lines(c(pos, pos), c(0, es1[pos]), col="grey")
            lines(es1, lwd=2)
            addDensityBar(which(names(signature) %in% names(geneset)), length(signature), col=colHit[1], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
            addSigBar(signature, colSig[1], base.sig, top.sig, colSig[3], colSig[4])
            if (!is.null(nes$nes)) {
                axis(3, length(signature)*.7, paste("NES: ", round(nes$nes, 2), ", p:", sep=""), hadj=1, tick=FALSE)
                axis(3, length(signature)*.71, niceExponent(nes$p.value), hadj=0, tick=FALSE)
            }
        }
        return(list(es=es, nes=nes$nes, p.value=nes$p.value, ledge=ledge, ledgesig=ledgesig))
    }
    if (prod(range(geneset))>=0) { #1-tail gsea on a +- signature
        signature <- sort(signature)
        es1 <- -gsea.es(signature, which(names(signature) %in% names(geneset)), score)
        # Check whether the analysis is 1 or 2 tails !!!!
        if (twoTails) {
            pos <- sort(c(which.min(es1), which.max(es1)))
            es <- diff(range(es1))
            # LEDGE
            ledgesig <- signature[c(1:pos[1], pos[2]:length(signature))]
            ledge <- ledgesig[names(ledgesig) %in% names(geneset)]
            nes <- list(nes=NULL, p.value=NULL)
            if (per>0) {
                nesnull <- sapply(1:per, function(i, signature, setsize, score) {
                    es1 <- -gsea.es(signature, sample(length(signature), setsize), score)
                    diff(range(es1))
                }, signature=signature, setsize=length(which(names(signature) %in% names(geneset))), score=score)
                nes <- aecdf(nesnull, symmetric=TRUE)(es, alternative)
            }
            if (pout) {
                if (is.null(ylim)) ylim <- c(min(0, min(es1)), max(0, max(es1)))
                if (length(ylim)==1) ylim <- c(-abs(ylim), abs(ylim))
                wd <- diff(ylim)
                base.bar <- ylim[2]+wd*.05
                top.bar <- ylim[2]+wd*.15
                base.sig <- ylim[1]-wd*.15
                top.sig <- ylim[1]-wd*.05
                plot(0, 0, type="n", xlim=c(0, length(signature)+1), ylim=c(base.sig, top.bar), axes=axes, xlab=xlab, ylab=ylab, main=main, font.main=1, cex.main=1, ...)
                abline(h=0, col="grey")
                pos <- which.max(es1)
                lines(c(pos, pos), c(0, es1[pos]), col="grey")
                pos <- which.min(es1)
                lines(c(pos, pos), c(0, es1[pos]), col="grey")
                lines(es1, lwd=2)
                addDensityBar(which(names(signature) %in% names(geneset)), length(signature), col=colHit[1], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
                addSigBar(signature, colSig[1:2], base.sig, top.sig, colSig[3], colSig[4])
                if (!is.null(nes$nes)) {
                    axis(3, length(signature)*.7, paste("NES: ", round(nes$nes, 2), ", p:", sep=""), hadj=1, tick=FALSE)
                    axis(3, length(signature)*.71, niceExponent(nes$p.value), hadj=0, tick=FALSE)
                }
            }
            return(list(es=es, nes=nes$nes, p.value=nes$p.value, ledge=ledge, ledgesig=ledgesig))
        }
        else {
            pos <- which.max(abs(es1))
            es <- es1[pos]
            # ledge
            if (es<0) ledgesig <- signature[1:pos]
            else ledgesig <- rev(signature[pos:length(signature)])
            ledge <- ledgesig[names(ledgesig) %in% names(geneset)]
            nes <- list(nes=NULL, p.value=NULL)
            if (per>0) {
                nesnull <- sapply(1:per, function(i, signature, setsize, score) {
                    es1 <- -gsea.es(signature, sample(length(signature), setsize), score)
                    es1[which.max(abs(es1))]
                }, signature=signature, setsize=length(which(names(signature) %in% names(geneset))), score=score)
                nes <- aecdf(nesnull, symmetric=TRUE)(es, alternative)
            }
            if (pout) {
                if (is.null(ylim)) ylim <- c(min(0, min(es1)), max(0, max(es1)))
                if (length(ylim)==1) ylim <- c(-abs(ylim), abs(ylim))
                wd <- diff(ylim)
                base.bar <- ylim[2]+wd*.05
                top.bar <- ylim[2]+wd*.15
                base.sig <- ylim[1]-wd*.15
                top.sig <- ylim[1]-wd*.05
                plot(0, 0, type="n", xlim=c(0, length(signature)+1), ylim=c(base.sig, top.bar), axes=axes, xlab=xlab, ylab=ylab, main=main, font.main=1, cex.main=1, ...)
                abline(h=0, col="grey")
                pos <- which.max(abs(es1))
                lines(c(pos, pos), c(0, es1[pos]), col="grey")
                lines(es1, lwd=2)
                addDensityBar(which(names(signature) %in% names(geneset)), length(signature), col=colHit[1], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
                addSigBar(signature, colSig[1:2], base.sig, top.sig, colSig[3], colSig[4])
                if (!is.null(nes$nes)) {
                    axis(3, length(signature)*.7, paste("NES: ", round(nes$nes, 2), ", p:", sep=""), hadj=1, tick=FALSE)
                    axis(3, length(signature)*.71, niceExponent(nes$p.value), hadj=0, tick=FALSE)
                }
            }
            return(list(es=es, nes=nes$nes, p.value=nes$p.value, ledge=ledge, ledgesig=ledgesig))
        }
    }
    # 2-tails GSEA
    es1 <- gsea2.es(signature, geneset, score)
    pos <- which.max(abs(es1$es))
    es <- es1$es[pos]
    # ledge
    if (es>0) {
        tmp <- es1$rlist[1:pos]
        ledge <- tmp[names(tmp) %in% names(geneset)]
        ledge <- ledge * sign(geneset)[match(names(ledge), names(geneset))]
        ledgesig=1
        #        ledgesig <- signature[c(1:max(which(names(signature) %in% names(ledge[ledge<0]))), min(which(names(signature) %in% names(ledge[ledge>0]))):length(signature))]
    }
    else {
        tmp <- rev(es1$rlist[pos:length(es1$rlist)])
        ledge <- tmp[names(tmp) %in% names(geneset)]
        ledge <- ledge * sign(geneset)[match(names(ledge), names(geneset))]
        ledgesig=1
        #        ledgesig <- signature[c(1:max(which(names(signature) %in% names(ledge[ledge>0]))), min(which(names(signature) %in% names(ledge[ledge<0]))):length(signature))]
    }
    nes <- list(nes=NULL, p.value=NULL)
    if (per>0) {
        nesnull <- sapply(1:per, function(i, signature, geneset, score) {
            names(geneset) <- sample(names(signature), length(geneset))
            es1 <- gsea2.es(signature, geneset, score)
            es1$es[which.max(abs(es1$es))]
        }, signature=signature, geneset=geneset, score=score)
        nes <- aecdf(nesnull, symmetric=TRUE)(es, alternative)
    }
    if (pout) {
        signature <- sort(signature)
        es1 <- gsea.es(signature, which(names(signature) %in% names(geneset[geneset<0])), score)
        es2 <- gsea.es(signature, which(names(signature) %in% names(geneset[geneset>=0])), score)
        if (is.null(ylim)) ylim <- range(c(es1, es2))
        if (length(ylim)==1) ylim <- c(-abs(ylim), abs(ylim))
        wd <- diff(ylim)
        base.sig <- ylim[1]-wd*.15
        top.sig <- ylim[1]-wd*.05
        plot(0, 0, type="n", xlim=c(0, length(signature)+1), ylim=c(base.sig, ylim[2]+wd*.25), axes=axes, xlab=xlab, ylab=ylab, main=main, font.main=1, cex.main=1, ...)
        abline(h=0, col="grey")
        pos <- which.max(abs(es1))
        if (colHit[1]>1) {
            lines(c(pos, pos), c(0, es1[pos]), col="grey")
            lines(es1, lwd=2)
        }
        else {
            lines(c(pos, pos), c(0, es1[pos]), col=hsv(colHit[1], .5, 1-colHit[3]))
            lines(es1, lwd=2, col=hsv(colHit[1], 1, 1-colHit[3]))
        }
        pos <- which.max(abs(es2))
        if (colHit[2]>1) {
            lines(c(pos, pos), c(0, es2[pos]), col="grey")
            lines(es2, lwd=2)
        }
        else {
            lines(c(pos, pos), c(0, es2[pos]), col=hsv(colHit[2], .5, 1-colHit[3]))
            lines(es2, lwd=2, col=hsv(colHit[2], 1, 1-colHit[3]))
        }
        base.bar <- ylim[2]+wd*.05
        top.bar <- ylim[2]+wd*.15
        addDensityBar(which(names(signature) %in% names(geneset[geneset<0])), length(signature), col=colHit[1], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
        base.bar <- ylim[2]+wd*.15
        top.bar <- ylim[2]+wd*.25
        addDensityBar(which(names(signature) %in% names(geneset[geneset>=0])), length(signature), col=colHit[2], base.bar=base.bar, top.bar=top.bar, mate=colHit[3], gama=colHit[4], lwd=lwd, maxhit=maxhit)
        addSigBar(signature, colSig[1:2], base.sig, top.sig, colSig[3], colSig[4])
        if (!is.null(nes$nes)) {
            axis(3, length(signature)*.7, paste("NES: ", round(nes$nes, 2), ", p:", sep=""), hadj=1, tick=FALSE)
            axis(3, length(signature)*.71, niceExponent(nes$p.value), hadj=0, tick=FALSE)
        }
    }
    return(list(es=es, nes=nes$nes, p.value=nes$p.value, ledge=ledge, ledgesig=ledgesig))
}

gsea.es <- function(rlist, x, score) {
    nr <- sum(abs(rlist[x])^score)
    nh <- length(rlist)-length(x)
    es <- rep(-(1/nh),length(rlist))
    es[x] <- abs(rlist[x])^score/nr
    return(cumsum(es))
}


gsea2.es <- function(rlist, x, score) {
    x1 <- names(x)[x>=0]
    x2 <- names(x)[x<0]
    px1 <- match(x1, names(sort(rlist, decreasing=T)))
    px2 <- match(x2, names(sort(rlist, decreasing=F)))
    names(px1) <- x1
    names(px2) <- x2
    rlistr <- rank(rlist)
    rlistr <- rlistr[!(names(rlist) %in% names(x))]
    rlistr <- sort(c(rlistr, px1, px2))
    rlist <- rlist[match(names(rlistr), names(rlist))]
    x <- which(names(rlist) %in% names(x))
    nr <- sum(abs(rlist[x])^score)
    nh <- length(rlist)-length(x)
    es <- rep(-(1/nh),length(rlist))
    es[x] <- abs(rlist[x])^score/nr
    return(list(es=cumsum(es), rlist=rlist))
}




addDensityBar <- function(x, sigLength, col, base.bar, top.bar, mate=.2, gama=2, lwd=1, maxhit=0) {
    den <- density(x, adj=.1)
    if (maxhit>0 & maxhit<length(x)) {
        x <- x[unique(round(seq(1, length(x), length=maxhit)))]
    }
    den <- approx(den, xout=x)
    den$y <- den$y/max(den$y)
    den <- lapply(den, function(x, pos) x[pos], pos=order(den$y))
    if (col>1) col <- hsv(0, 0, (1-den$y^gama))
    else col <- hsv(col, den$y^gama, 1-mate*den$y)
    for (i in 1:length(x)) lines(c(den$x[i], den$x[i]), c(base.bar, top.bar), col=col[i], lwd=lwd)
    rect(0, base.bar, sigLength+1, top.bar)
}

addSigBar <- function(signature, col, base.sig, top.sig, mate=.1, gama, res=1000) {
    x <- signature/max(abs(signature))
    if (length(col)==1) col <- rep(col, 2)
    col1 <- rep(hsv(0, 0, 0), length(x))
    pos <- x<0
    if (col[1]>1) col1[pos] <- hsv(0, 0, (1-abs(x[pos])^gama))
    else col1[pos] <- hsv(col[1], abs(x[pos])^gama, 1-mate*abs(x[pos]))
    if (col[2]>1) col1[!pos] <- hsv(0, 0, (1-x[!pos])^gama)
    else col1[!pos] <- hsv(col[2], x[!pos]^gama, 1-mate*x[!pos])
    for (i in unique(round(seq(1, length(x), length=res)))) lines(c(i, i), c(base.sig, top.sig), col=col1[i])
    rect(0, base.sig, length(x)+1, top.sig)
}
