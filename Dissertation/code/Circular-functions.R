#=================================================================================================
# THE JONES-PEWSEY DISTRIBUTION
#=================================================================================================
# Jones-Pewsey normalising constant: evaluated separately for reasons of efficiency
JP.NCon <- function(kappa, psi){
    
    # if kappa small, assume circular uniform distribution
    if (kappa < 0.001) {
        ncon <- 1/(2*pi) 
        return(ncon)
    } else {
        eps <- 10*.Machine$double.eps
        
        # if psi small, assume von Mises distribution
        if (abs(psi) <= eps) {
            ncon <- 1/(2*pi*I.0(kappa))
            return(ncon)
        } else {
            # numerical integration of normalising constant
            intgrnd <- function(x){ (cosh(kappa*psi)+sinh(kappa*psi)*cos(x))**(1/psi) }
            ncon <- 1/integrate(intgrnd, lower=-pi, upper=pi)$value
            return(ncon)
        } }
}

# Jones-Pewsey density function
JP.pdf <- function(theta, mu, kappa, psi, ncon){
    
    if (kappa < 0.001) {pdfval <- 1/(2*pi) ; return(pdfval)}
    else {
        # if psi small, use von Mises density
        eps <- 10*.Machine$double.eps
        if (abs(psi) <= eps) {
            pdfval <- ncon*exp(kappa*cos(theta-mu))
            return(pdfval)
        } else { 
            pdfval <- (cosh(kappa*psi)+sinh(kappa*psi)*cos(theta-mu))**(1/psi)
            pdfval <- ncon*pdfval
            return(pdfval) 
        } }
}

# Jones-Pewsey distribution function
JP.df <- function(theta, mu, kappa, psi, ncon) {
    
    eps <- 10*.Machine$double.eps
    if (theta <= eps) {
        dfval <- 0 
        return(dfval)
    } else {
        if (theta >= 2*pi-eps) {
            dfval <- 1 
            return(dfval)
        } else {
            if (kappa < 0.001) {
                dfval <- theta/(2*pi) 
                return(dfval)
            } else {
                if (abs(psi) <= eps) {
                    vMPDF <- function(x){ ncon*exp(kappa*cos(x-mu)) }
                    dfval <- integrate(vMPDF, lower = 0, upper=theta)$value
                    return(dfval)
                } else {
                    dfval <- integrate(JP.pdf, mu, kappa, psi, ncon, lower = 0, upper=theta)$value
                    return(dfval)
                } } } }
}

# Jones-Pewsey quantile function
JP.qf <- function(u, mu, kappa, psi, ncon) {
    
    eps <- 10*.Machine$double.eps
    if (u <= eps) {theta <- 0 ; return(theta)}
    
    else 
        if (u >= 1-eps) {
            theta <- 2*pi-eps
            return(theta)
        } else {
            if (kappa < 0.001) {
                theta <- u*2*pi 
                return(theta)
            } else {
                roottol <- .Machine$double.eps**(0.6)
                qzero <- function(x) {
                    y <- JP.df(x, mu, kappa, psi, ncon) - u ; return(y) }
                res <- uniroot(qzero, lower=0, upper=2*pi-eps, tol=roottol)
                theta <- res$root 
                return(theta) 
            } }
}