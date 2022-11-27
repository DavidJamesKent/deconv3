##' computeM
##'
##' Computes \eqn{q \times q}{q × q} matrix M with entries \eqn{M_{ij} = \int (g*B_i)(g*B_j)}{Mᵢⱼ = ∫(g*Bᵢ)(g*Bⱼ)}
##' 
##' @title Computes Gram matrix of convolved B-splines
##' @param gtwid Fourier transform of g
##' @param qdim Dimension of spline basis
##' @param knots Location of knots
##' @param h Knot spacing
##' @return 
##' @author David Kent
##' @export
computeM <- function(gtwid,qdim,knots,h) {
    M <- matrix(NA,qdim,qdim); dqq <- dim(M)
    for(i in 0:(qdim-1)) {
        M[abs(.row(dqq)-.col(dqq))==i] <- stats::integrate(
            function(t) gtwid(t)^2*
                        Re(
                            btwid(t,knots[1],h)*
                            Conj(btwid(t,knots[1+i],h))
                        ),-Inf,Inf)$value/(2*pi)
    }
    return(M)
}

##' Computes the matrix N for the quadratic programming problem
##'
##' Computes \eqn{K\times q}{K × q} matrix with entries
##' \deqn{N_{ij} = \int_a^b g*B_j,}{Nᵢⱼ = ∫ₐᵇg*Bⱼ,}
##' where \eqn{a = m_i-\delta/2}{a = mᵢ-δ/2} and \eqn{b = m_i+\delta/2}{b = mᵢ+δ/2}, with \eqn{m_i}{mᵢ} the midpoint of the ith histogram bin.
##' 
##' @title computeN
##' @param basis The B-spline basis
##' @param midp The K midpoints of the histogram
##' @param gtwid The Fourier transform of the error density
##' @param K The number of histogram midpoints
##' @param qdim The dimension of the B-Spline basis
##' @return 
##' @author David Kent
##' @export
computeN <- function(basis,midp,gtwid,K,qdim) {
    h <- basis$h
    delt <- attr(midp,"delt")
    knots <- basis$knots
    s <- attr(gtwid,"s")

    tmp1 <- function(x) abs(Re(gtwid(x)*btwid(x,knots[1],h)))
    xs <- seq(0,100,l=1e3)
    lim <- xs[max(which(tmp1(xs) > 1e-4))]
    gl <- pracma::gaussLegendre(15,0,lim) #integrand seems to be even?
    N <- matrix(0,K,qdim)
    for(i in 1:K) {
        for(j in 1:qdim) {
            sig <- sqrt(4*h^2/12 + s^2) # "SD" of a rv which has density b-spline convolved with g
            if( abs(midp[i] - (knots[j] + 2*h)) < 6*sig ) {
                tmp2 <- function(t) Re(gtwid(t)*btwid(t,knots[j],h)*
                                      Conj(histogtwid(t,midp[i],delt)))
                N[i,j] <- 2*sum(tmp2(gl$x)*gl$w)/(2*pi)*delt
            }
        }
    }
    return(N)
}
