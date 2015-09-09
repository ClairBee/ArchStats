require(circular)

#=================================================================================================
# TESTS OF REFLECTIVE SYMMETRY
#=================================================================================================
# Pewsey's test of reflective symmetry
r.symm.test.stat <- function(data) {
    n <- length(data)
    bar <- get.moments(data)
    
    # est. variance of b_2
    var <- ((1-bar$a4)/2-(2*bar$a2)+(2*bar$a2/bar$r) * (bar$a3 + (bar$a2*(1-bar$a2)/bar$r)))/n
    
    # test statistic
    ts <- abs(bar$b2/sqrt(var))
    
    pval <- 2*pnorm(ts, mean = 0, sd = 1, lower = F)
    list(test.statistic = ts, p.val = pval)
}

# Bootstrap version of test
r.symm.test.boot <- function(data, B = 9999, alpha = 0.05, show.progress = T) {
    
    n <- length(data)
    
    absz <- r.symm.test.stat(data)$test.statistic           # observed test statistic
    
    tbar <- mean(data)
    reflection <- 2*tbar-data
    symm.data <- c(data, reflection)                        # symmetrize data
    
    # take bootstrap samples from symmetrized data, calculate test statistics
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 2:(B+1)) { 
        bootstrap.sample <- sample(symm.data, size = n, replace = TRUE)
        absz[b] <- r.symm.test.stat(bootstrap.sample)$test.statistic
        if (show.progress) {setTxtProgressBar(pb, b)}
    }
    
    # estimate p-value
    p.val = length(absz[absz >= absz[1]]) / (B + 1)
    if (show.progress) {close(pb)}
    list(p.val = p.val,
         std.error = qnorm(1-(alpha/2)) * sqrt((p.val*(1-p.val))/(B + 1)))
}
#=================================================================================================
# PARAMETER ESTIMATION
#=================================================================================================
# support function: extract string of named trigonometric moments
get.moments <- function(data) {
    
    t2 <- trigonometric.moment(data, p = 2, center = T)
    t3 <- trigonometric.moment(data, p = 3, center = T)
    t4 <- trigonometric.moment(data, p = 4, center = T)
    
    list(mu = as.numeric(mean(data)),
         r = rho.circular(data),
         b2 = t2$sin,
         b3 = t3$sin,
         b4 = t4$sin,
         a2 = t2$cos,
         a3 = t3$cos,
         a4 = t4$cos)
}

# bias-corrected sample statistics
bc.sample.statistics <- function(data) {
    
    n <- length(data)
    bar <- get.moments(data)    
    
    rho.est <- bar$r - ((1-bar$a2)/(4*n*bar$r))
    r2 <- bar$r^2
    r4 <- r2^2
    
    beta2.est <- bar$b2 + ((bar$b3/bar$r)+(bar$b2/r2)-(2*bar$a2*bar$b2/r4))/n
    
    mu.est <- (bar$mu + (bar$b2/(2*n*r2))) %% (2*pi)
    
    alpha2.est <- bar$a2 - (1-(bar$a3/bar$r)-((bar$a2*(1-bar$a2)+bar$b2*bar$b2)/r2))/n
    
    # return list of parameters
    list(mu = mu.est, rho = rho.est, beta2 = beta2.est, alpha2 = alpha2.est)   
}

# large-sample (1-alpha)% confidence interval for bias-corrected sample statistics
bc.ci.LS <- function(data, alpha = 0.05) {
    
    n <- length(data)
    ests <- bc.sample.statistics(data)
    bar <- get.moments(data)
    qval <- qnorm(1-alpha/2)
    
    r2 <- bar$r^2 ; r4 <- r2^2
    
    rho.bc <- ests$rho
    r.SE <- sqrt((1-2*r2+bar$a2)/(2*n))    
        rho.upper <- rho.bc + qval*r.SE
        rho.lower <- rho.bc - qval*r.SE
    rho <- c(estimate = rho.bc, lower = rho.lower, upper = rho.upper)
    
    beta2.bc <- ests$beta2
    beta2.SE <- sqrt((((1-bar$a4)/2)-(2*bar$a2)-(bar$b2^2)+
                          (2*bar$a2/bar$r)*(bar$a3+(bar$a2*(1-bar$a2)/bar$r)))/n)
        beta2.upper <- beta2.bc + qval*beta2.SE
        beta2.lower <- beta2.bc - qval*beta2.SE
    beta2 <- c(estimate = beta2.bc, lower = beta2.lower, upper = beta2.upper)
    
    mu.bc <- ests$mu
    mu.SE <- sqrt((1-bar$a2)/(2*n*r2))
        mu.upper <- (mu.bc + qval * mu.SE) %% (2*pi)
        mu.lower <- (mu.bc - qval * mu.SE) %% (2*pi)
    mu <- c(estimate = mu.bc, lower = mu.lower, upper = mu.upper)
    
    alpha2.bc <- ests$alpha2
    alpha2.SE <- sqrt((((1+bar$a4)/2)-(bar$a2*bar$a2)+
                           (2*bar$b2/bar$r)*(bar$b3+(bar$b2*(1-bar$a2)/bar$r)))/n)
        alpha2.upper <- alpha2.bc + qval*alpha2.SE
        alpha2.lower <- alpha2.bc - qval*alpha2.SE
    alpha2 <- c(estimate = alpha2.bc, lower = alpha2.lower, upper = alpha2.upper)
    
    list(alpha = alpha, mu = mu, rho = rho, beta2 = beta2, alpha2 = alpha2)
}

#-------------------------------------------------------------------------------------------------
# maximum likelihood estimation of Jones-Pewsey parameters
JP.mle <- function(data) {
    n <- length(data)
    s <- sum(sin(data))
    c <- sum(cos(data))
    
    # starting values: von Mises ML estimates
    muvM <- atan2(s,c) %% (2*pi)
    kapvM <- A1inv(sqrt(s*s+c*c)/n)
    
    # specify Jones-Pewsey log-likelihood function to optimize
    JPnll <- function(p){
        mu <- p[1] ; kappa <- p[2] ; psi <- p[3] ; parlim <- abs(kappa*psi)
        if (parlim > 10) {
            y <- 9999.0
            return(y)
        } else {
            ncon <- JP.NCon(kappa, psi)
            y <- -sum(log(JP.pdf(data, mu, kappa, psi, ncon)))
            return (y)
        }
    }
    
    # optimize the specified function - try-catch added to capture optimization errors
    out <- tryCatch({
        optim(par=c(muvM, kapvM, 0.001), fn=JPnll, gr = NULL, method = "L-BFGS-B",
              lower = c(muvM-pi, 0, -Inf), upper = c(muvM+pi, Inf, Inf), hessian = TRUE)
    }, 
    error=function(err) {list(value = NA, par = rep(NA, 3), HessMat = NA)} )
    
    list(maxll = -out$value,
         mu = out$par[1] %% (2*pi),
         kappa = out$par[2],
         psi = out$par[3],
         HessMat = out$hessian)
}

# approximate normal-theory confidence intervals for Jones-Pewsey MLE, 
# obtained by solving the Hessian matrix to obtain the squared standard error
JP.ci.nt <- function(jp.ests, alpha = 0.05) {
    quant <- qnorm(1-alpha/2)
    infmat <- solve(jp.ests$HessMat)
    standerr <- sqrt(diag(infmat))
    
    list(alpha = alpha,
         mu = c(est = jp.ests$mu, lower = jp.ests$mu-(quant*standerr[1]), 
                                  upper = jp.ests$mu+(quant*standerr[1])),
         kappa = c(est = jp.ests$kappa, lower = jp.ests$kappa-(quant*standerr[2]), 
                                        upper = jp.ests$kappa+(quant*standerr[2])),
         psi = c(est = jp.ests$psi, lower = jp.ests$psi-(quant*standerr[3]), 
                                        upper = jp.ests$psi+(quant*standerr[3]))    )
}

#=================================================================================================
# ASSESSING GOODNESS OF FIT
#=================================================================================================
# probability integral transform test for goodness of fit of von Mises distribution
# (parametric bootstrap version)
vM.GoF.boot <- function(data, mu.0, kappa.0, B = 9999, show.progress = T) {
    
    # get 2*pi*F(theta)
    n <- length(data)
    tdf <- pvonmises(data, mu.0, kappa.0, from=circular(0), tol = 1e-06)
    cunif <- circular(2*pi*tdf)
    
    unif.test.0 <- rep(0,2)
    nxtrm <- rep(1,2)
    
    unif.test.0[1] <- kuiper.test(cunif)$statistic
    unif.test.0[2] <- watson.test(cunif)$statistic
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 2:(B+1)) {
        # create bootstrap sample, obtain parameters
        bootstrap.sample <- rvonmises(n, mu.0, kappa.0)
        vM.mle <- mle.vonmises(bootstrap.sample, bias = T)
        mu.1 <- vM.mle$mu
        kappa.1 <- vM.mle$kappa
        
        tdf <- pvonmises(bootstrap.sample, mu.1, kappa.1, from=circular(0), tol = 1e-06)
        cunif <- circular(2*pi*tdf)
        
        # count number of tests where bootstrap test statistic is greater than observed
        nxtrm[1] <- nxtrm[1] + (kuiper.test(cunif)$statistic >= unif.test.0[1])
        nxtrm[2] <- nxtrm[2] + (watson.test(cunif)$statistic >= unif.test.0[2])
        
        if (show.progress) {setTxtProgressBar(pb,b)}      
    }
    p.val <- nxtrm/(B+1)
    names(p.val) <- c("kuiper", "watson")
    if (show.progress) {close(pb)}      
    p.val
}

# probability integral transform test for goodness of fit of Jones-Pewsey distribution
# (parametric bootstrap version)
JP.GoF.boot <- function(data, mu, kappa, psi, B = 9999, show.progress = T) {
    n <- length(data)
    muhat0 <- mu
    kaphat0 <- kappa
    psihat0 <- psi
    ncon0 <- JP.NCon(kaphat0, psihat0)
    
    # 2*pi*F(theta)
    tdf <- 0 
    for (j in 1:n) {
        tdf[j] <- JP.df(data[j], muhat0, kaphat0, psihat0, ncon0)
    }
    
    cunif <- circular(2*pi*tdf)
    nxtrm <- rep(1,2)
    
    unif.test.0 <- c(kuiper.test(cunif)$statistic, 
                     watson.test(cunif)$statistic)
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = B, style = 3)}
    for (b in 2:(B+1)) {
        # generate bootstrap samples and parameter estimates
        bootstrap.sample <- JP.sim(n, muhat0, kaphat0, psihat0, ncon0)
        JPmleRes <- JP.mle(bootstrap.sample) 
        
        # if J-P MLE can't be found, treat as rejected test
        if (!is.na(JPmleRes$mu)) {
            muhat1 <- JPmleRes$mu
            kaphat1 <- JPmleRes$kappa
            psihat1 <- JPmleRes$psi
            ncon1 <- JP.NCon(kaphat1, psihat1)
            tdf <- 0
            for (j in 1:n) {
                tdf[j] <- JP.df(bootstrap.sample[j], muhat1, kaphat1, psihat1, ncon1)
            }
            cunif <- circular(2*pi*tdf)
            
            nxtrm[1] <- nxtrm[1] + (kuiper.test(cunif)$statistic >= unif.test.0[1])
            nxtrm[2] <- nxtrm[2] + (watson.test(cunif)$statistic >= unif.test.0[2])
            if (show.progress) {setTxtProgressBar(pb, b)}
        }}
    pval <- nxtrm/(B+1)
    names(pval) <- c("kuiper", "watson")
    if (show.progress) {close(pb)}
    return(pval)
}

#-------------------------------------------------------------------------------------------------
# von Mises P-P plot and residuals
vM.PP <- function(data, mu, kappa)  {
    edf <- ecdf(data)
    tdf <- pvonmises(data, mu, kappa, from=circular(0), tol = 1e-06)
    
    plot.default(tdf, edf(data), pch=20, xlim=c(0,1), ylim=c(0,1),
        xlab = "von Mises distribution function", ylab = "Empirical distribution function")
    lines(c(0,1), c(0,1), lwd=2, col = "lightseagreen")
    # return residuals
    edf(data) - tdf
}

# von Mises Q-Q plot and residuals
vM.QQ <- function(data, mu, kappa)  {
    edf <- ecdf(data)
    tqf <- qvonmises(edf(data), mu, kappa, from=circular(0), tol = 1e-06)
    
    plot.default(tqf, data, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), 
                 xlab = "von Mises quantile function", ylab = "Empirical quantile function")
    lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "lightseagreen")
    # return residuals
    data - tqf
}

#-------------------------------------------------------------------------------------------------
# Jones-Pewsey P-P plot and residuals
JP.PP <- function(data, mu, kappa, psi) {
    ncon <- JP.NCon(kappa, psi)     # normalising constant
    n <- length(data)
    edf <- ecdf(data) 
    tdf <- 0 
    for (j in 1:n) {tdf[j] <- JP.df(data[j], mu, kappa, psi, ncon)}
    
    plot.default(tdf, edf(data), pch=20, xlim=c(0,1), ylim=c(0,1),
         xlab = "Jones-Pewsey distribution function", ylab = "Empirical distribution function")
    lines(c(0,1), c(0,1), lwd=2, col = "lightseagreen")
    # return residuals
    edf(data) - tdf
}

# Jones-Pewsey Q-Q plot and residuals
JP.QQ <- function(data, mu, kappa, psi) {
    ncon <- JP.NCon(kappa, psi)
    n <- length(data)
    edf <- ecdf(data) 
    tqf <- 0 
    for (j in 1:n) {tqf[j] <- JP.qf(edf(data)[j], mu, kappa, psi, ncon)}
    
    plot.default(tqf, data, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi),
                 xlab = "Jones-Pewsey quantile function", ylab = "Empirical quantile function") 
    lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "lightseagreen")
    # return residuals
    data - tqf
}


#=================================================================================================
# COMPARISON OF MULTIPLE SAMPLES
#=================================================================================================
# Watson's test of common mean direction
watson.common.mean.test <- function(samples) {
    
    data <- unlist(samples)
    N <- length(data)
    g <- length(samples)
    sample.sizes <- 0
    for (i in 1:g) {sample.sizes[i] <- length(samples[[i]])}
    
    size.csum <- cumsum(sample.sizes) 
    
    delhat <- 0 
    tbar <- 0
    
    for (k in 1:g) {
        sample <- samples[[k]]
        
        tm1 <- trigonometric.moment(sample, p=1)
        tm2 <- trigonometric.moment(sample, p=2)
        
        Rbar1 <- tm1$rho
        Rbar2 <- tm2$rho 
        tbar[k] <- tm1$mu %% (2*pi)
        
        delhat[k] <- (1-Rbar2)/(2*Rbar1*Rbar1)
    }
    
    dhatmax <- max(delhat) 
    dhatmin <- min(delhat)
    
    if (dhatmax/dhatmin <= 4) { # use P procedure
        
        CP <- 0
        SP <- 0
        dhat0 <- 0
        
        for (k in 1:g) {
            CP <- CP + sample.sizes[k]*cos(tbar[k])
            SP <- SP + sample.sizes[k]*sin(tbar[k])
            dhat0 <- dhat0 + sample.sizes[k]*delhat[k] 
        }
        
        dhat0 <- dhat0/N
        RP <- sqrt(CP*CP+SP*SP)
        
        Yg <- 2*(N-RP)/dhat0
    } else {
        
        CM <- 0 
        SM <- 0
        Yg <- 0
        
        for (k in 1:g) {
            CM <- CM + (sample.sizes[k]*cos(tbar[k])/delhat[k])
            SM <- SM + (sample.sizes[k]*sin(tbar[k])/delhat[k])
            Yg <- Yg + (sample.sizes[k]/delhat[k]) 
        }
        RM <- sqrt(CM*CM+SM*SM)
        Yg <- 2*(Yg-RM)
    }
    
    pval = pchisq(Yg, g-1, lower.tail = F)
    list(Y.g = Yg, p.val = pval, disp.ratio = dhatmax/dhatmin)
}

# Wallraff test of common concentration
wallraff.concentration.test <- function(samples) {
    
    data <- unlist(samples)
    g <- length(samples)
    g.id <- c()
    for (i in 1:g) {
        g.id <- c(g.id, rep(i, length(samples[[i]])))
    }
    
    N <- length(data)
    #    sample.sizescsum <- cumsum(sample.sizes) 
    
    tbar <- circular(0) 
    distdat <- c()
    for (k in 1:g) {
        
        #        dist <- 0 
        sample <- samples[[k]]
        
        tm1 <- trigonometric.moment(sample, p=1) 
        tbar[k] <- tm1$mu %% (2*pi)
        
        dist <- pi-abs(pi-abs(sample-tbar[k]))
        distdat <- c(distdat, dist)
    }
    
    TestRes <- kruskal.test(distdat, g = g.id)
    
    list(p.val = TestRes$p.value, result = TestRes)
} 

# support function to calculate sine and consine rank scores
cs.unif.scores <- function(samples) {
    
    data <- unlist(samples)
    N <- length(data)
    ranks <- rank(data, ties.method="random")
    cos.u.scores <- cos(ranks*2*pi/N)
    sin.u.scores <- sin(ranks*2*pi/N)
    
    list(cos.scores = cos.u.scores, sin.scores = sin.u.scores)    
}

# Mardia-Wheeler-Watson test of commono distribution
mww.common.dist.LS <- function(cs.scores, sample.sizes) {
    
    N <- sum(sample.sizes)
    
    g <- length(sample.sizes)
    g.id <- c()
    for (i in 1:g) {
        g.id <- c(g.id, rep(i, sample.sizes[i]))
    }
    
    Wg <- 0
    
    for (k in 1:g) {
        cos.u.scores.k <- cs.scores$cos.scores[g.id == k] 
        sin.u.scores.k <- cs.scores$sin.scores[g.id == k] 
        
        sum.cos.k.sq <- (sum(cos.u.scores.k))^2
        sum.sin.k.sq <- (sum(sin.u.scores.k))^2
        
        Wg <- Wg+(sum.cos.k.sq+sum.sin.k.sq)/length(g.id[g.id == k])
    }
    Wg <- 2*Wg 
    p.val <- pchisq(Wg, 2*(g-1), lower.tail = F)
    
    list(statistic = Wg, p.val = p.val)    
}

# randomisation version of Watson's two-sample test
watson.two.test.rand <- function(data1, data2, NR = 9999, show.progress = T){
    
    wats.obs <- watson.two.test(data1, data2)$statistic 
    nxtrm <- 1
    
    n1 <- length(data1)
    n2 <- length(data2)
    N <- n1+n2
    
    combined <- c(data1, data2)
    
    if (show.progress) {pb <- txtProgressBar(min = 0, max = NR, style = 3)}
    
    for (r in 1:NR) {
        rs <- sample(combined)
        rs1 <- rs[1:n1]
        rs2 <- rs[(n1+1):N]
        
        wats.rand <- watson.two.test(rs1, rs2)$statistic
        
        if (wats.rand >= wats.obs) {nxtrm <- nxtrm+1}
        if (show.progress) {setTxtProgressBar(pb, r)}
    }
    
    pval <- nxtrm/(NR+1) 
    if (show.progress) {close(pb)}
    list(observed = wats.obs, p.val = pval)
}