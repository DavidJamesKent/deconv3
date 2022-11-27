##' Computes estimated risk function
##'
##' 
##' @title bootRisk
##' @param control Object computed by `computeControl()`
##' @param dat Object computed by `computeData()`
##' @param thetal0 Estimated \eqn{\theta}{θ} from \eqn{\lambda^{(k)}}{λ⁽ᵏ⁾} for use in bootstrap
##' @param mix Uses Gramian mix[1]*G + mix[2]*M, which corresponds to
##'     a weighted average of risks in f-space and h-space,
##'     respectively
##' @return `logR` is a function for which logR(l) returns the
##'     estimated MISE of the estimator with penalty 10^l
##' @author David Kent
##' @export
bootRisk <- function(control,dat,thetal0,mix=c(1,0)) {
    n <- length(dat$Y)
    G <- control$basis$G
    M <- control$M
    P <- control$P
    N <- control$N
    midp <- control$midp
    delt <- attr(midp,"delt")
    breaks <- attr(midp,"breaks")
    phat <- as.vector(N%*%thetal0); phat <- phat/sum(phat)
    mu <- phat/delt
    V <- (diag(phat) - phat%o%phat)/(n*delt^2)

    gramian <- mix[1]*G + mix[2]*M

    logR <- Vectorize(function(l) {
        Sl <- solve(M + 10^l*P,t(N))
        Rl <- (t(thetal0)%*%(gramian)%*%thetal0 -
               2*t(thetal0)%*%(gramian)%*%Sl%*%mu +
               sum(diag(t(Sl)%*%(gramian)%*%Sl%*%V)) +
               t(mu)%*%t(Sl)%*%(gramian)%*%Sl%*%mu)[1]
        return(log10(Rl))
    })
    return(logR)
}

##' Chooses lambda
##'
##' @title chooseLambda
##' @param control Object computed by `computeControl()`
##' @param dat Object computed by `computeData()`
##' @param eps Convergence criterion
##' @param maxit Maximum number of iterations
##' @param ll0 Initial \eqn{\lambda}{λ}
##' @param verbose Output information at each iteration
##' @param mix Two-element vector describing weighting of MISE in f (first element) vs. MISE of g*f (second)
##' @return \eqn{\lambda}{λ} chosen by iterated parametric bootstrap
##' @author David Kent
##' @export
chooseLambda <- function(control,dat,eps=1e-4,maxit=100,ll0=-1,verbose=FALSE,mix=c(1,0)) {
    llkm1 <- ll0
    thetakm1 <- attr(doSolve(control,dat,10^llkm1),"theta")
    logR <- bootRisk(control,dat,thetakm1,mix=mix)
    llk <- stats::optim(-8,logR,method="Brent",lower=-8,upper=5)$par
    k <- 1
    while( abs(llk-llkm1)/abs(llk) > eps && k < maxit ) {
        k <- k+1
        llkm1 <- llk
        thetakm1 <- attr(doSolve(control,dat,10^llkm1),"theta")
        logR <- bootRisk(control,dat,thetakm1,mix=mix)
        llk <- stats::optim(-8,logR,method="Brent",lower=-8,upper=5)$par
        if(verbose) cat(sprintf("%d: log(lambda) = %2.3f\n",k,llk))
    }
    return(10^llk)
}
