# $Id: LD.R 453 2005-11-09 17:04:02Z warnes $

# R translation of Cathy Stack's SAS macro
# Assumes 2-alleles

LD <- function(g1,...)
  UseMethod("LD",g1)

LD.data.frame <- function(g1,...)
  {
    gvars <- sapply( g1, function(x) (is.genotype(x) && nallele(x)==2) )
    if(any(gvars==FALSE))
      {
        warning("Non-genotype variables or genotype variables ",
                "with more or less than two alleles detected. ",
                "These variables will be omitted: ",
                paste( colnames(g1)[!gvars] , collapse=", " )
                )
        g1 <- g1[,gvars]
      }


    P <- matrix(nrow=ncol(g1),ncol=ncol(g1))
    rownames(P) <- colnames(g1)
    colnames(P) <- colnames(g1)

    P <- D <- Dprime <- nobs <- chisq <- p.value <- corr <- R.2 <- P

    for(i in 1:(ncol(g1)-1) )
      for(j in (i+1):ncol(g1) )
        {
          ld <- LD( g1[,i], g1[,j] )

          D      [i,j] <- ld$D
          Dprime [i,j] <- ld$"D'"
          corr   [i,j] <- ld$"r"
          R.2    [i,j] <- ld$"R^2"
          nobs   [i,j] <- ld$"n"
          chisq  [i,j] <- ld$"X^2"
          p.value[i,j] <- ld$"P-value"
        }

    retval <- list(
                   call=match.call(),
                   "D"=D,
                   "D'"=Dprime,
                   "r" = corr,
                   "R^2" = R.2,
                   "n"=nobs,
                   "X^2"=chisq,
                   "P-value"=p.value
           )

    class(retval) <- "LD.data.frame"

    retval
  }

LD.genotype <- function(g1,g2,...)
  {
    if(is.haplotype(g1) || is.haplotype(g2))
      stop("Haplotype options are not yet supported.")

    if(nallele(g1)!=2 || nallele(g2)!=2)
      stop("This function currently only supports 2-allele genotypes.")

    prop_g1 <- summary(g1)$allele.freq[,2]
    prop_g2 <- summary(g2)$allele.freq[,2]

    increase_g1 <- names(prop_g1)[which.min(prop_g1)]
    increase_g2 <- names(prop_g2)[which.min(prop_g2)]
    pG1 <- max(prop_g1, na.rm=TRUE)
    pG2 <- max(prop_g2, na.rm=TRUE)
    pg1 <- 1-pG1
    pg2 <- 1-pG2

    Dm_minus <- max(-pG1*pG2, -pg1*pg2)
    pmin <- pG1*pG2 + Dm_minus;

    Dm_plus <- min(pG1*pg2, pG2*pg1);
    pmax <- pG1*pG2 + Dm_plus;

    counts <- table(
                    allele.count(g1, increase_g1),
                    allele.count(g2, increase_g2)
                    )

    # reverse the matrix, row-wise and column-wise
#    require(magic)
#    n3x3 = arev(counts)


	# can be visualized in a diagram
    loglik <- function(pG1G2,...)
      {
        (2*counts[1,1]+counts[1,2]+counts[2,1])*log(pG1G2) +
        (2*counts[1,3]+counts[1,2]+counts[2,3])*log(pG1-pG1G2) +
        (2*counts[3,1]+counts[2,1]+counts[3,2])*log(pG2-pG1G2) +
        (2*counts[3,3]+counts[3,2]+counts[2,3])*log(1-pG1-pG2+pG1G2) +
        counts[2,2]*log(pG1G2*(1-pG1-pG2+pG1G2) + (pG1-pG1G2)*(pG2-pG1G2)) # pG1G2*pg1g2 + pG1g2*pg1G2
      }

    solution <- optimize(
                         loglik,
                         lower=pmin+.Machine$double.eps,
                         upper=pmax-.Machine$double.eps,
                         maximum=TRUE
                         )
    pG1G2 <- solution$maximum

    estD <- pG1G2 - pG1*pG2
    if (estD>0)
      estDp <- estD / Dm_plus
    else
      estDp <- estD / Dm_minus

    n <-  sum(counts)

    corr <- estD / sqrt( pG1 * pG2 * pg1 * pg2 )

    dchi <- (2*n*estD^2)/(pG1 * pg1 * pG2* pg2)
    dpval <- 1 - pchisq(dchi,1)

    retval <- list(
                   call=match.call(),
                   "D"=estD,
                   "D'"=estDp,
                   "r" = corr,
                   "R^2" = corr^2,
                   "n"=n,
                   "X^2"=dchi,
                   "P-value"=dpval
                   )

    class(retval) <- "LD"

    retval

  }


