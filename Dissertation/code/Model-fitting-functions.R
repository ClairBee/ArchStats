require(circular)

#===========================================================================================
# TESTS OF REFLECTIVE SYMMETRY
#===========================================================================================
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

#===========================================================================================
# PARAMETER ESTIMATION
#===========================================================================================
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
    beta2.SE <- sqrt((((1-bar$a4)/2)-(2*bar$a2)-(bar$b2^2)+(2*bar$a2/bar$r)*(bar$a3+(bar$a2*(1-bar$a2)/bar$r)))/n)
        beta2.upper <- beta2.bc + qval*beta2.SE
        beta2.lower <- beta2.bc - qval*beta2.SE
    beta2 <- c(estimate = beta2.bc, lower = beta2.lower, upper = beta2.upper)
    
    mu.bc <- ests$mu
    mu.SE <- sqrt((1-bar$a2)/(2*n*r2))
        mu.upper <- (mu.bc + qval * mu.SE) %% (2*pi)
        mu.lower <- (mu.bc - qval * mu.SE) %% (2*pi)
    mu <- c(estimate = mu.bc, lower = mu.lower, upper = mu.upper)
    
    alpha2.bc <- ests$alpha2
    alpha2.SE <- sqrt((((1+bar$a4)/2)-(bar$a2*bar$a2)+(2*bar$b2/bar$r)*(bar$b3+(bar$b2*(1-bar$a2)/bar$r)))/n)
        alpha2.upper <- alpha2.bc + qval*alpha2.SE
        alpha2.lower <- alpha2.bc - qval*alpha2.SE
    alpha2 <- c(estimate = alpha2.bc, lower = alpha2.lower, upper = alpha2.upper)
    
    list(alpha = alpha, mu = mu, rho = rho, beta2 = beta2, alpha2 = alpha2)
}

#----------------------------------------------------------------------------------------
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
    
    # optimize the specified function
    out <- optim(par=c(muvM, kapvM, 0.001), fn=JPnll, gr = NULL, method = "L-BFGS-B", lower = c(muvM-pi, 0, -Inf), upper = c(muvM+pi, Inf, Inf), hessian = TRUE)
    
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
         mu = c(est = jp.ests$mu, lower = jp.ests$mu-(quant*standerr[1]), upper = jp.ests$mu+(quant*standerr[1])),
         kappa = c(est = jp.ests$kappa, lower = jp.ests$kappa-(quant*standerr[2]), upper = jp.ests$kappa+(quant*standerr[2])),
         psi = c(est = jp.ests$psi, lower = jp.ests$psi-(quant*standerr[3]), upper = jp.ests$psi+(quant*standerr[3]))    )
}

#===========================================================================================
# ASSESSING GOODNESS OF FIT
#===========================================================================================
# probability integral transform test for goodness of fit of von Mises distribution
# includes direct distributional Watson test of goodness of fit
vM.GoF <- function(data, mu, kappa, display = T) {
    tdf <- pvonmises(data, circular(mu), kappa, from=circular(0), tol = 1e-06)  
    cunif <- circular(2*pi*tdf)         # 2*pi*F(theta)                                         
    
    # only print test results if specified (otherwise, only return test statistics)
    if (display) {
        print(watson.test(data, dist = "vonmises")) # distributional test
        print(kuiper.test(cunif))
        print(watson.test(cunif))
    } else {
        list(watson.vM.stat = watson.test(data, dist = "vonmises")$statistic,
             kuiper.unif.stat = kuiper.test(cunif)$statistic,
             watson.unif.stat = watson.test(cunif)$statistic,
    }      
}

# probability integral transform test for goodness of fit of Jones-Pewsey distribution
JP.GoF <- function(data, mu, kappa, psi, display = T) {
    n <- length(data)
    tdf <- 0
    ncon <- JP.NCon(kappa, psi)
    
    for (j in 1:n) {tdf[j] <- JP.df(data[j], mu, kappa, psi, ncon)}     # F(theta)
    cunif <- circular(2*pi*tdf)                                         # 2*pi*F(theta)
    
    if (display) {
        print(kuiper.test(cunif))
        print(watson.test(cunif))
    } else {
        list(kuiper.statistic = kuiper.test(cunif)$statistic,
             watson.statistic = watson.test(cunif)$statistic)
    }
}

#----------------------------------------------------------------------------------------
# von Mises P-P plot and residuals
vM.PP <- function(data, mu, kappa)  {
    edf <- ecdf(data)
    tdf <- pvonmises(data, mu, kappa, from=circular(0), tol = 1e-06)
    
    plot.default(tdf, edf(data), pch=20, xlim=c(0,1), ylim=c(0,1), xlab = "von Mises distribution function", ylab = "Empirical distribution function")
    lines(c(0,1), c(0,1), lwd=2, col = "lightseagreen")
    # return residuals
    edf(data) - tdf
}

# von Mises Q-Q plot and residuals
vM.QQ <- function(data, mu, kappa)  {
    edf <- ecdf(data)
    tqf <- qvonmises(edf(data), mu, kappa, from=circular(0), tol = 1e-06)
    
    plot.default(tqf, data, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), xlab = "von Mises quantile function", ylab = "Empirical quantile function")
    lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "lightseagreen")
    # return residuals
    data - tqf
}

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
    
    plot.default(tqf, data, pch=20, xlim=c(0,2*pi), ylim=c(0,2*pi), xlab = "Jones-Pewsey quantile function", ylab = "Empirical quantile function") 
    lines(c(0,2*pi), c(0,2*pi), lwd=2, col = "lightseagreen")
    # return residuals
    data - tqf
}


#===========================================================================================
# EXPECTATION-MAXIMISATION ALGORITHM
#===========================================================================================
# E-M algorithm to fit a mixture of k von Mises distributions
EM.vonmises <- function(x, k, max.runs = 1000, conv = 0.00001) {
    
    x <- circular(x)
    # provide starting values for mu, kappa, alpha
    mu <- circular(runif(k, 0, max(x)))
    kappa <- runif(k,0,1)
    alpha <- runif(k,0,1)
    alpha <- alpha/sum(alpha) # normalise to sum to 1
    
    # Support function - calculate log-likelihood
    log.likelihood <- function(x, mu, kappa, alpha, k) {
        l <- matrix(nrow = k, ncol = length(x))
        for (i in 1:k) {
            l[i,] <- alpha[i] * dvonmises(x, mu[i], kappa[i])
        }
        sum(log(colSums(l)))
    }
    
    log.lh <- log.likelihood(x, mu, kappa, alpha, k)
    new.log.lh <- abs(log.lh) + 100
    n = 0
    
    while ((abs(log.lh - new.log.lh) > conv) && (n < max.runs)) {
        
        # Estimation - calculate z_ij
        z <- matrix(0, ncol = length(x), nrow = k)
        for (i in 1:k){
            z[i,] <- alpha[i] * dvonmises(x, mu[i], kappa[i])
        }
        all.z <- colSums(z)
        for (i in 1:k){
            z[i,] <- z[i,] / all.z
        }
        
        # Maximisation - update parameters
        for (i in 1:k) {
            alpha[i] <- sum(z[i,]) / length (x)
            mu[i] <- atan2(sum(z[i,] * sin(x)), sum(z[i,] * cos(x)))
            kappa[i] <- A1inv(sum(z[i,] * (cos(x - mu[i]))) / sum(z[i,]))
            
            # correct for negative kappa if necessary
            if (kappa[i] < 0) {
                kappa[i] <- abs(kappa[i])
                mu[i] <- mu[i] + pi
            }
        }
        
        # calculate log-likelihoods for comparison
        log.lh <- new.log.lh
        new.log.lh <- log.likelihood(x, mu, kappa, alpha, k)
        n <- n + 1
    }
    
    # Output: if model hasn't converged, show error message 
    if ((abs(log.lh - new.log.lh) > conv)) {
        cat ("Data hasn't converged after", n, "iterations; \n",
             "Difference in log-likelihoods is", 
             round(abs(log.lh - new.log.lh),6))
    } else {
        list(k = k, mu = mu %% (2*pi), kappa = kappa, alpha = alpha, log.lh = new.log.lh)
    }
}

# modified algorithm to fit a mixture of 1 uniform and (k-1) von Mises
EM.u.vonmises <- function(x, k, max.runs = 1000, conv = 0.00001) {
    
    # E-M algorithm with one component fixed as uniform
    x <- circular(x)
    # provide starting values for mu, kappa, alpha
    mu <- circular(runif(k, 0, max(x)))
    kappa <- c(0, runif(k-1,0,1))
    alpha <- runif(k,0,1)
    alpha <- alpha/sum(alpha) # normalise to sum to 1
    
    # Support function - calculate log-likelihood
    
    log.likelihood <- function(x, mu, kappa, alpha, k) {
        l <- matrix(nrow = k, ncol = length(x))
        for (i in 1:k) {
            l[i,] <- alpha[i] * dvonmises(x, mu[i], kappa[i])
        }
        sum(log(colSums(l)))
    }
    
    log.lh <- log.likelihood(x, mu, kappa, alpha, k)
    new.log.lh <- abs(log.lh) + 100
    n = 0
    
    while ((abs(log.lh - new.log.lh) > conv) && (n < max.runs)) {
        
        # Estimation - calculate z_ij
        z <- matrix(0, ncol = length(x), nrow = k)
        for (i in 1:k){
            z[i,] <- alpha[i] * dvonmises(x, mu[i], kappa[i])
        }
        all.z <- colSums(z)
        for (i in 1:k){
            z[i,] <- z[i,] / all.z
        }
        
        # Maximisation - update parameters
        for (i in 1:k) {
            alpha[i] <- sum(z[i,]) / length (x)
            mu[i] <- atan2(sum(z[i,] * sin(x)), sum(z[i,] * cos(x)))
            kappa[i] <- A1inv(sum(z[i,] * (cos(x - mu[i]))) / sum(z[i,]))   
            
            # correct for negative kappa if necessary
            if (kappa[i] < 0) {
                kappa[i] <- abs(kappa[i])
                mu[i] <- mu[i] + pi
            }
        }
        
        # Sort parameters by kappa for identifiability
        alpha <- alpha[order(kappa)]
        mu <- mu[order(kappa)]
        z <- z[order(kappa),]
        kappa <- kappa[order(kappa)]
        # fix smallest kappa as 0 (uniform)
        kappa[1] <- 0
        
        # calculate log-likelihoods for comparison
        log.lh <- new.log.lh
        new.log.lh <- log.likelihood(x, mu, kappa, alpha, k)
        n <- n + 1
    }
    
    # Output: if model hasn't converged, show error message
    if ((abs(log.lh - new.log.lh) > conv)) {
        cat ("Data hasn't converged after", n, "iterations; \n",
             "Difference in log-likelihoods is", 
             round(abs(log.lh - new.log.lh),6))
    } else {
        list(k = k, mu = mu %% (2*pi), kappa = kappa, alpha = alpha,
             log.lh = new.log.lh)
    }
}

# plot von Mises mixture model and components
plot.EM.vonmises <- function(varToPlot, modelToPlot, h.breaks = 20) {
    
    varToPlot = matrix(varToPlot)       # convert from circular data
    
    # get density of each component
    x <- circular(seq(0, 2*pi, 0.01))
    components <- matrix(nrow = modelToPlot$k, ncol = length(x))
    for (i in 1:modelToPlot$k) {
        components[i,] <- dvonmises(x, circular(modelToPlot$mu[i]), modelToPlot$kappa[i]) * modelToPlot$alpha[i]
    }
    mixt <- colSums(components)
    
    # plot original data
    y.max <- max(c(mixt, hist(varToPlot, plot = F, breaks = h.breaks)$density)) * 1.1 # rescale y axis
    labl <- paste("Mixture of", modelToPlot$k,"von Mises")
    hist(varToPlot, freq = F, ylim = c(0, y.max), main = "", col = "lightgrey", xlab = labl, xlim = c(0, 2*pi), breaks = h.breaks)
    
    # plot vM components
    for (i in 1:modelToPlot$k) {
        lines(matrix(x), components[i,], col = i+1, lwd = 2)
    }
    # add mixture model
    lines(matrix(x), mixt, lwd = 2)
}

# obtain cluster numbers using winner-takes-all assignment over mixture model
mvM.clusters <- function(data, model) {
    components <- matrix(nrow = model$k, ncol = length(data))
    for (i in 1:model$k) {
        components[i, ] <- dvonmises(data, circular(model$mu[i]), model$kappa[i])
    }
    apply(components, 2, which.max)
}

# mixture von Mises P-P plot and residuals
mvM.PP <- function(data, mu, kappa, alpha) {
    edf <- ecdf(data)
    
    l <- matrix(nrow = length(mu), ncol = length(data))
    for (i in 1:length(mu)) {
        l[i,] <- alpha[i] * pvonmises(data, mu[i], kappa[i], from = circular(0), tol = 1e-06)
    }
    tdf <- colSums(l)
    
    plot.default(tdf, edf(data), pch = 20, xlim = c(0, 1), ylim = c(0, 1), xlab = "mixture von Mises distribution function", ylab = "Empirical distribution function")
    lines(c(0, 1), c(0, 1), lwd = 2, col = "lightseagreen")
    edf(data) - tdf
}
