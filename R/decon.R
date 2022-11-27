##' Computes bin midpoints and densities.
##'
##' Padding is not typically useful, but in circumstances where new
##' data will be used with the same `control` object
##' (e.g. simulations), some padding can help avoid the necessity of
##' recomputing the histogram bins and associated inner product
##' matrices.
##'
##' 
##' @title computeHist
##' @param Y The data; a numeric vector
##' @param K The number of bin midpoints
##' @param pad Distance (in units of Y) beyond [min(Y),max(Y)] to compute bins
##' @param padperc Similar to `pad`, but expressed as a percentage of the range of Y
##' @return 
##' @author David Kent
##' @export
computeHist <- function(Y,K,pad=NULL,padperc=0.0) {
    ymi <- min(Y); yma <- max(Y)
    if( !is.null(pad) & !is.null(padperc) ) {
        stop("Only one of `pad` and `padperc` may be specified")
    } else if( !is.null(padperc) ) {
        if( padperc < 0 )
            stop("`padperc` must be non-negative")
        pad <- (yma-ymi)*padperc
    } else if( !is.null(pad) & pad < 0 ) {
        stop("`pad` must be non-negative")
    } else {
        pad <- 0
    }

    histmi <- ymi-pad; histma <- yma+pad
    delt <- (histma - histmi)/(K-1)
    breaks <- seq(histmi-delt/2,histma+delt/2,delt)
    yhist <- graphics::hist(Y,breaks=breaks,plot=F)
    midp <- yhist$mids
    hhat <- yhist$density
    attr(midp,"delt") <- delt
    attr(midp,"breaks") <- breaks
    attr(midp,"midpa") <- midp[1]
    attr(midp,"midpb") <- midp[K]
    return(list(hhat=hhat,midp=midp))
}

##' Computes constraints for quadratic programming problem.
##'
##' In either the constrained quadratic program, or the metric
##' projection, the PDF and positive constraints are required by
##' `quadprog` to be expressed as linear functions of the
##' coefficients, and those matrices and vectors are computed
##' here. Since the positivity constraint is only approximate, we need
##' the B-spline basis and a grid density on which to express the
##' constraint.
##'
##' 
##' @title computeConstraints
##' @param B The B-spline basis
##' @param nxeval Number of points x at which to compute B(x)
##' @return 
##' @author David Kent
##' @export
computeConstraints <- function(B,nxeval=1e3) {
    qdim <- attr(B,"qdim")
    rangeval <- attr(B,"rangeval")
    xeval <- seq(rangeval[1],rangeval[2],l=nxeval)
    Bxeval <- t(B(xeval))
    ## Always have pdf constraints
    Asum2one <- rep(1,qdim); bsum2one <- 1
    #Apos <- diag(qdim); bpos <- rep(0,qdim)
    Apos <- Bxeval; bpos <- rep(0,nxeval)
    A <- cbind(Asum2one,Apos)
    bvec <- c(bsum2one,bpos)
    meq <- 1

    return(list(A=A,bvec=bvec,meq=meq))
}

##' Computes all the details required for the estimator
##'
##' Utility function which computes the histogram grid, B-spline
##' basis, inner- and cross-product matrices, and linear constraint
##' matrices.
##' 
##' The objects computed here do *not* depend on the data except that
##' the data is used to determine the intervals on which the histogram
##' and spline basis are computed; thus, if new data is acquired
##' (e.g. through simulations), which are consistent with those
##' intervals, these objects may all be re-used.
##' 
##' @title computeControl
##' @return 
##' @author David Kent
##' @param Y The data
##' @param K The number of histogram bins
##' @param qdim The dimension of the spline space
##' @param gtwid The Fourier transform of the error density
##' @param padperc How far beyond the data the spline support should
##'     extend, as a proportion of the range of Y
##' @export
computeControl <- function(Y,K,qdim,gtwid,padperc) {
    histog <- computeHist(Y,K)
    mi <- min(Y) - padperc*diff(range(Y))
    ma <- max(Y) + padperc*diff(range(Y))
    basis <- createBasis(c(mi,ma),qdim)

    M <- computeM(gtwid,qdim,basis$knots,basis$h)
    N <- computeN(basis,histog$midp,gtwid,K,qdim)

    constraints <- computeConstraints(basis$B,nxeval=1e3)

    return(list(M=M,
                N=N,
                P=basis$P,
                midp=histog$midp,
                A=constraints$A,
                bvec=constraints$bvec,
                meq=constraints$meq,
                qdim=qdim,
                K=K,
                gtwid=gtwid,
                basis=basis))
}
##' (Re-)computes `hhat` on the `midp` grid in `control`
##'
##' 
##' @title computeData
##' @param Y The (possibly new) data
##' @param control A `control` object from computeControl()
##' @return A list containing two objects: `$hhat` is the histogram bin-heights, and `$Y` is the data.
##' @author David Kent
##' @export
computeData <- function(Y,control) {
    breaks <- attr(control$midp,"breaks")
    hhat <- graphics::hist(Y,breaks=breaks,plot=F)$density
    return(list(hhat=hhat,Y=Y))
}

##' Does solve, given `control`, `dat`, and `lambda`
##'
##' With all the pre-computed pieces from computeControl() and
##' computeData(), solve the minimization problem over the spline
##' space.
##' 
##' @title doSolve
##' @param control Control object produced by computeControl()
##' @param dat Data object produced by computeData()
##' @param lambda Penalty parameter
##' @param unconstrained If TRUE, do not apply PDF constraints
##' @param metric If TRUE, compute metric projection of unconstrained
##'     solution; otherwise, apply constraint directly in minimization
##'     problem.
##' @return 
##' @author David Kent
##' @export
doSolve <- function(control,dat,lambda,unconstrained=FALSE,metric=TRUE) {
    if( !metric ) {
        print("Doing non-metric-projection method")
        theta <- quadprog::solve.QP(control$M + lambda*control$P,
                                    dat$hhat%*%control$N,
                                    Amat=control$A,bvec=control$bvec,meq=control$meq)$solution
    } else {
        thetaunc <- solve(control$M + lambda*control$P,t(control$N)%*%dat$hhat)
        if( !unconstrained ) {
            theta <- quadprog::solve.QP(control$basis$G,
                                        t(thetaunc)%*%control$basis$G,
                                        Amat=control$A,bvec=control$bvec,
                                        meq=control$meq)$solution
        } else {
            theta <- thetaunc
        }
    }
    fhat <- function(x) as.vector(control$basis$B(x)%*%theta)
    attr(fhat,"theta") <- theta
    attr(fhat,"lambda") <- lambda
    return(fhat)
}

##' Deconvolves a density
##'
##' Given data Y = X + e, where e has density g, computes an estimate
##' of the density of X. The Fourier transform of the density of e
##' must be provided as `gtwid`.
##' 
##' @title decon
##' @param Y The data, i.e. Y = X + e
##' @param gtwid A function which returns the Fourier transform of the density of e
##' @param lambda Either "risk" to choose a penalty parameter automatically, or a positive number
##' @param K (optional) number of bins for histogram
##' @param qdim (optional) dimension of spline space
##' @param control (optional) control object from previous call to decon()
##' @param padperc (optional) How far beyond the data the spline
##'     support should extend, as a proportion of the range of Y
##' @param unconstrained (optional) If TRUE, returns unconstrained
##'     solutionmetric (optional) If TRUE, projects
##'     unconstrained solution onto PDFs. If FALSE, computes solution
##'     as constrained minimization problem.
##' @param metric Compute metric projection? If FALSE, adds constraint directly to minimization problem
##' @param verbose (optional) If TRUE, prints extra information.
##' @usage decon(Y,gtwid=NULL,lambda="risk",K=NULL,qdim=NULL,
##'     control=NULL,padperc=0.1,unconstrained=FALSE,
##'     metric=TRUE,verbose=FALSE)
##' @return A list, where `$fhat` is a function which evaluates the
##'     estimate of the density of X, `$dat` contains the data used to
##'     compute that estimate, and `$control` contains the
##'     miscellaneous objects used to compute `$fhat`. `$control` can
##'     be supplied to another call of decon() so that it does not
##'     have to be recomputed.
##' @author David Kent
##' @export
decon <- function(Y,gtwid=NULL,lambda="risk",K=NULL,qdim=NULL,
                  control=NULL,padperc=0.1,unconstrained=FALSE,
                  metric=TRUE,verbose=FALSE) {
    if( is.null(control) ) {
        if( is.null(gtwid) )
            stop("If `control` is not specified, an error density Fourier transform `gtwid` must be supplied")
        if( is.null(K) ) {
            K <- floor(diff(range(Y))/(stats::IQR(Y)*2/length(Y)^(1/3))) ## Freedman-Diaconis rule
            if(verbose)
                cat(sprintf("Automatically choosing K=%d\n",K))
        }
        if( is.null(qdim) ) {
            qdim <- max(50,K)
            if(verbose)
                cat(sprintf("Automatically choosing qdim=%d\n",qdim))
        }
        control <- computeControl(Y,K,qdim,gtwid,padperc=padperc)
        dat <- computeData(Y,control)
    } else {
        ## TODO: Check condition where K/qdim AND control are supplied
        ## If data are outside of the support of our histograms, just recompute the whole thing
        if( min(Y) < attr(control$midp,"midpa") | max(Y) > attr(control$midp,"midpb") )
            control <- computeControl(Y,control$K,control$qdim,control$gtwid,padperc=padperc)
        dat <- computeData(Y,control)
    }

    if( is.numeric(lambda) & length(lambda) == 1 ) {
        fhat <- doSolve(control,dat,lambda,unconstrained=unconstrained,metric=metric)
    } else if( lambda == "risk" ) {
        l <- chooseLambda(control,dat,verbose=verbose)
        fhat <- doSolve(control,dat,l,unconstrained=unconstrained,metric=metric)
    }
    return(list(fhat=fhat,dat=dat,control=control))
}

##' Sinc function
##'
##' Computes sinc(x)
##'
##' @param x Location to evaluation
##' @return sin(x)/x
##' @title sinc
##' @export
sinc <- function(x) ifelse(x!=0,1*sin(x)/x,1)

##' Unscaled B-Spline Fourier Transform
##' 
##' Given fourth-order integer-knot B-Splines, evaluates the Fourier transform of B_{0,4}.
##' 
##' @param t Argument of Fourier Transform
##' @title Unscaled B-Spline Fourier Transform
##' @return (complex-valued) Fourier transform of B_{0,4} evaluated at t
##' @export
unscaledtwid <- function(t) (exp(-1i*t*0.5)*sinc(0.5*t))^4

##' Scaled B-Spline Fourier Transform
##' 
##' Evaluates Fourier transform of an integer-knot B-Spline scaled by
##' h and translated so that the left-most support begins at k1.
##' 
##' @param t Argument of Fourier Transform
##' @param k1 Leftmost side of basis function
##' @param h Knot spacing
##' @title Scaled and Translated B-Spline Fourier Transform
##' @return (complex-valued) Fourier transform of scaled and translated B-Spline
##' @export
btwid <- function(t,k1,h) exp(-1i*k1*t)*unscaledtwid(t*h)

##' Fourier transform of histogram basis function
##'  
##' The basis functions of a histogram are indicator functions for the
##' bins, scaled by the bin-spacing delt. This function evaluates the
##' Fourier transform of the basis function at midpoint `midp` with
##' spacing `delt`.
##'  
##' @param t Argument of Fourier transform
##' @param midp Midpoint of histogram bin
##' @param delt Width of histogram bin
##' @title Histogram Basis Function Fourier Transform
##' @return (complex-valued) Fourier transform of histogram basis function evaluated at t
##' @export
histogtwid <- function(t,midp,delt) exp(-1i*t*(midp-delt/2))*exp(-1i*t*0.5*delt)*sinc(0.5*delt*t)

##' Returns Fourier Transform of specified error density
##'
##' For Gaussian density \eqn{g(x) = \sqrt{2\pi s^2}\exp{-x^2/2s^2}}{g(x) = √(2πs²)·exp{ x²/2s² }}, returns
##' \deqn{F[g](t) = \int e^{-itx} g(x)\,dx = \exp{-s^2t^2/2}}{F[g](t) = ∫e⁻ⁱᵗˣg(x)dx = exp(-s²t²/2)}
##' 
##' For Laplace density \eqn{g(x) = \exp{-|x|/b}/2b}{g(x) = exp{ -|x|/b }/2b}, returns
##' \deqn{F[g](t) = \int e^{-itx} g(x)\,dx = 1/(1+b^2t^2)}{F[g](t) = ∫e⁻ⁱᵗˣg(x)dx = 1/(1+b²t²)}
##'
##' For Uniform density \eqn{g(x) = I[-a,a]/2a}{g(x) = Ind[-a,a]/2a}, returns
##' \deqn{F[g](t) = \int e^{-itx}g(x)\,dx = sinc(at)}{F[g](t) = ∫e⁻ⁱᵗˣg(x)dx = sinc(at)}
##' 
##' @title getgtwid()
##' @param param Named numeric vector containing an entry named 's'
##'     for the standard deviation of a Gaussian, 'b' for the scale of
##'     the Laplace, or 'a' for the interval half-width for the
##'     Uniform.
##' @param g Character string, one of "Gaussian", "Laplace", or
##'     "Uniform"
##' @return Function which evaluates the Fourier transform of the
##'     specified density
##' @author David Kent
##' @export
getgtwid <- function(g="Gaussian",param) {
    accepts <- c("Gaussian","Laplace","Uniform")
    if( !(g %in% accepts) ) {
        stop(paste("g must be one of",paste(accepts,collapse=", ")))
    } else if( g == "Gaussian" ) {
        if( is.na(param['s']) )
            stop(sprintf("Gaussian error has exactly one parameter, standard deviation `s`, but `param` has no element named `s`."))
        s <- unname(param['s'])
        four <- function(t) exp(-s^2*t^2/2)
        attr(four,"s") <- s
        return(four)
    } else if( g == "Laplace" ) {
        if( is.na(param['b']) )
            stop(sprintf("Laplace error has exactly one parameter, scale `b`, but `param` has no element named `b`."))
        b <- unname(param['b'])
        four <- function(t) 1/(1+b^2*t^2)
        attr(four,"s") <- sqrt(2)*b
        return(four)
    } else if( g == "Uniform" ) {
        if( is.na(param['a']) )
            stop(sprintf("Uniform error has exactly one parameter, interval half-width `a`, but `param` has no element named `a`."))
        a <- unname(param['a'])
        four <- function(t) sinc(a*t)
        attr(four,"s") <- a^2/3
        return(four)
    }
}
