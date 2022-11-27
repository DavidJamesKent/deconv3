##' Constructs q-dimensional basis on an interval [a, b]
##'
##' 
##' @title createBasis
##' @param rangeval Numeric vector of length 2 indicating the interval
##'     [a, b] on which the basis may be evaluated.
##' @param qdim Integer dimension of basis
##' @return A list, where `$B`, `$Bp`, and `$Bpp` evaluate the basis
##'     functions and their first two derivatives, respectively, `$G`
##'     and `$P` are inner-product matrices for `$B` and `$Bpp`,
##'     respectively, and `knots` and `h` are the knots and
##'     knot-spacing associated with the basis.
##' @author David Kent
##' @export
createBasis <- function(rangeval, qdim) {
    if( !(is.numeric(rangeval) & length(rangeval)==2) )
        stop("`rangeval` must be numeric of length 2")
    if( qdim != as.integer(qdim) )
        stop("`qdim` must be an integer value")
    if( qdim < 7 )
        stop("`qdim` must be greater than or equal to 7")

    knots <- seq(rangeval[1],rangeval[2],l=qdim+4)
    h <- knots[2] - knots[1]
    B <- function(x) splines::splineDesign(knots,x,outer.ok=TRUE)/h
    Bp <- function(x) splines::splineDesign(knots,x,outer.ok=TRUE,derivs=1)/h
    Bpp <- function(x) splines::splineDesign(knots,x,outer.ok=TRUE,derivs=2)/h

    attr(B,"rangeval") <- attr(Bp,"rangeval") <- attr(Bpp,"rangeval") <- rangeval
    attr(B,"qdim")     <- attr(Bp,"qdim")     <- attr(Bpp,"qdim") <- qdim
    attr(B,"h")        <- attr(Bp,"h")        <- attr(Bpp,"h") <- h

    ## Inner-product matrix for f
    G <- matrix(0,qdim,qdim); dqq <- dim(G)
    G[.row(dqq) == .col(dqq)] <- 151/(315*h)
    G[abs(.row(dqq) - .col(dqq)) == 1] <- 397/(1680*h)
    G[abs(.row(dqq) - .col(dqq)) == 2] <- 1/(42*h)
    G[abs(.row(dqq) - .col(dqq)) == 3] <- 1/(5040*h)

    ## Inner product matrix for s'
    Z <- matrix(0,qdim,qdim)
    Z[.row(dqq) == .col(dqq)] <- 2/(3*h^3)
    Z[abs(.row(dqq) - .col(dqq)) == 1] <- -1/(8*h^3)
    Z[abs(.row(dqq) - .col(dqq)) == 2] <- -1/(5*h^3)
    Z[abs(.row(dqq) - .col(dqq)) == 3] <- -1/(120*h^3)

        ## Inner product matrix for s''
    P <- matrix(0,qdim,qdim)
    P[.row(dqq) == .col(dqq)] <- 8/(3*h^5)
    P[abs(.row(dqq) - .col(dqq)) == 1] <- -3/(2*h^5)
    P[abs(.row(dqq) - .col(dqq)) == 2] <- 0
    P[abs(.row(dqq) - .col(dqq)) == 3] <- 1/(6*h^5)

    ####
    ## Sparse versions of the above -- only faster if qdim > 200 or so
    ####
    ## G <- Matrix::bandSparse(qdim,qdim,c(0,-1,1,-2,2,-3,3),
    ##                               diagonals=list(rep(151/(315*h),qdim),
    ##                                              rep(397/(1680*h),qdim-1),
    ##                                              rep(397/(1680*h),qdim-1),
    ##                                              rep(1/(42*h),qdim-2),
    ##                                              rep(1/(42*h),qdim-2),
    ##                                              rep(1/(5040*h),qdim-3),
    ##                                              rep(1/(5040*h),qdim-3)))
    ## P <- Matrix::bandSparse(qdim,qdim,c(0,-1,1,-3,3),
    ##                               diagonals=list(rep(8/(3*h^5),qdim),
    ##                                              rep(-3/(2*h^5),qdim-1),
    ##                                              rep(-3/(2*h^5),qdim-1),
    ##                                              rep(1/(6*h^5),qdim-3),
    ##                                              rep(1/(6*h^5),qdim-3)))

    return(list(B=B,Bp=Bp,Bpp=Bpp,G=G,P=P,Z=Z,knots=knots,h=h))
}
