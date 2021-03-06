
# fn mledist() from {fitdistrplus}

mledist <- function (data,  # format: df with 2 cols [left] and [right]; if right-censored, the [right] col will be NA 
                     distr,  # character string, e.g. "weibull" 
                     start = NULL,  # initial guess 
                     fix.arg = NULL,  # params to hold fixed 
                     optim.method = "default", 
                     lower = -Inf,  # relevant for L-BFGS-B
                     upper = Inf,   # relevant for L-BFGS-B
                     custom.optim = NULL, 
                     weights = NULL, 
                     silent = TRUE, 
                     gradient = NULL,  # fn to return gradients for BFGS, CG, L-BFGS-B methods. If null, will use a finite-difference approach  
                     checkstartfix = FALSE, ...) {
    
    if (!is.character(distr)) 
        stop("distr must be a character string naming a distribution")
    else distname <- distr
    
    # get dweibull(), the weibull pdf ----
    ddistname <- paste("d", distname, sep = "")
    argddistname <- names(formals(ddistname))
    
    if (!exists(ddistname, mode = "function")) 
        stop(paste("The ", ddistname, " function must be defined"))
    
    if (is.null(custom.optim)) 
        optim.method <- match.arg(optim.method, c("default", 
                                                  "Nelder-Mead", 
                                                  "BFGS", "CG", "L-BFGS-B", 
                                                  "SANN", "Brent"))
    
    start.arg <- start  # starting guess for params 
    if (is.vector(start.arg)) 
        start.arg <- as.list(start.arg)
    
    txt1 <- "data must be a numeric vector of length greater than 1 for non censored data"
    txt2 <- "or a dataframe with two columns named left and right and more than one line for censored data"
    
    # dealing with weights: ----
    if (!is.null(weights)) {
        if (any(weights < 0)) 
            stop("weights should be a vector of integers greater than 0")
        if (!is.allint.w(weights)) 
            stop("weights should be a vector of (strictly) positive integers")
        if (length(weights) != NROW(data)) 
            stop("weights should be a vector with a length equal to the observation number")
        warning("weights are not taken into account in the default initial values")
    }
    
    
    if (is.vector(data)) {
        cens <- FALSE
        if (!(is.numeric(data) & length(data) > 1)) 
            stop(paste(txt1, txt2))
    }
    else {
        cens <- TRUE
        censdata <- data
        if (!(is.vector(censdata$left) & is.vector(censdata$right) & 
              length(censdata[, 1]) > 1)) 
            stop(paste(txt1, txt2))
        
        # get pweibull(), the weibull CDF ----  
        pdistname <- paste("p", distname, sep = "")
        if (!exists(pdistname, mode = "function")) 
            stop(paste("The ", pdistname, " function must be defined to apply maximum likelihood to censored data"))
    }
    
    
    if (cens) {
        dataformat <- cens2pseudo(censdata)  # todo: what is this? 
        data <- dataformat$pseudo
        rcens <- dataformat$rcens  # vector of right-censored obs? 
        lcens <- dataformat$lcens
        icens <- dataformat$icens
        ncens <- dataformat$ncens  # vector of uncensored obs? 
        irow <- cens2idxrow(censdata)
        irow.rcens <- irow$rcens
        irow.lcens <- irow$lcens
        irow.icens <- irow$icens
        irow.ncens <- irow$ncens
    }
    
    
    if (!checkstartfix) {
        arg_startfix <- manageparam(start.arg = start, 
                                    fix.arg = fix.arg, 
                                    obs = data,
                                    distname = distname)
        hasnodefaultval <- sapply(formals(ddistname), is.name)
        arg_startfix <- checkparamlist(arg_startfix$start.arg, 
                                       arg_startfix$fix.arg,
                                       argddistname,
                                       hasnodefaultval)
        if (is.function(fix.arg)) 
            fix.arg.fun <- fix.arg
        else fix.arg.fun <- NULL
    }
    else {
        arg_startfix <- list(start.arg = start, fix.arg = fix.arg)
        fix.arg.fun <- NULL
    }
    
    
    vstart <- unlist(arg_startfix$start.arg)
    
    if (is.null(vstart)) 
        stop("Starting values could not be NULL with checkstartfix=TRUE")
    
    fix.arg <- arg_startfix$fix.arg
    
    # special case: uniform dist: ---------
    if (distname == "unif") {
        par <- c(min = min(data), max = max(data))
        res <- list(estimate = par[!names(par) %in% names(fix.arg)], 
                    convergence = 0, 
                    loglik = NA, 
                    hessian = NA, 
                    optim.function = NA, 
                    fix.arg = fix.arg)
        return(res)
    }
    
    # case 1: uncensored data with no weights: ----
    if (!cens && is.null(weights)) {
        if ("log" %in% argddistname) {
            fnobj <- function(par, fix.arg, obs, ddistnam) {
                -sum(do.call(ddistnam, c(list(obs), as.list(par), 
                                         as.list(fix.arg), log = TRUE)))
            }
        }
        else {
            # > function to optimize: -----
            
            # this is just the standard log-likelihood function, based on
            # product of pdf values at each data point 
            fnobj <- function(par, 
                              fix.arg, 
                              obs,
                              ddistnam) {
                
                -sum(log(do.call(ddistnam, 
                                 c(list(obs), 
                                   as.list(par),
                                   as.list(fix.arg)
                                   )
                                 )
                         )
                     )
            }
        }
    }
    
    # case 2: censored data with no weights: ----
    else if (cens && is.null(weights)) {
        argpdistname <- names(formals(pdistname))
        
        if (("log" %in% argddistname) & ("log.p" %in% 
                                         argpdistname)) 
            
            # > function to optimize: -------
            fnobjcens <- function(par, 
                                  fix.arg, 
                                  rcens, 
                                  lcens, 
                                  icens, 
                                  ncens, 
                                  ddistnam, 
                                  pdistnam) {
                
                loglik_1 <- -sum(do.call(ddistnam,
                                         c(list(ncens),  # only uncensored 
                                           as.list(par),
                                           as.list(fix.arg),
                                           list(log = TRUE)
                                           )
                                         )
                                 ) 
                
                
                loglik_2 <- -sum(do.call(pdistnam,
                                          c(list(lcens),  # only left-censored 
                                            as.list(par),
                                            as.list(fix.arg),
                                            list(log = TRUE)
                                            )
                                          )
                                  ) 
                
                loglik_3 <- -sum(do.call(pdistnam, 
                                          c(list(rcens),  # only right-censored 
                                            as.list(par),
                                            as.list(fix.arg),
                                            list(lower.tail = FALSE),
                                            list(log = TRUE)
                                            )
                                          )
                                  ) 
                
                loglik_4 <- -sum(log(do.call(pdistnam, 
                                             c(list(icens$right),  # only interval-censored? 
                                               as.list(par),
                                               as.list(fix.arg)
                                               )
                                             ) - do.call(pdistnam, 
                                                         c(list(icens$left), 
                                                           as.list(par),
                                                           as.list(fix.arg)
                                                           )
                                                         )
                                     )
                                 )
                
                return(loglik_1 + loglik_2 + loglik_3 + loglik_4)
            }
                
        else {
            
            #___________________________----
            # main function to optimize: -----
            #___________________________----
            
            # in the case of right-censoring only, this is the same as Harrell, p426 
            
            fnobjcens <- function(par, 
                                  fix.arg, 
                                  rcens, 
                                  lcens, 
                                  icens, 
                                  ncens, 
                                  ddistnam, 
                                  pdistnam) {
                
                loglike_1 <- -sum(log(do.call(ddistnam,  # pdf 
                                             c(list(ncens),  # uncensored 
                                               as.list(par),
                                               as.list(fix.arg)
                                               )
                                             )
                                     )
                                 )
                
                loglike_2 <- -sum(log(do.call(pdistnam,  # cdf 
                                              c(list(lcens),  # left-censored 
                                                as.list(par),
                                                as.list(fix.arg)
                                                )
                                              )
                                      )
                                  )
                
                loglike_3 <- -sum(log(1 - do.call(pdistnam,  # survival function 
                                                  c(list(rcens),  # right-censored 
                                                    as.list(par),
                                                    as.list(fix.arg)
                                                    )
                                                  )
                                      )
                                  ) 
                
                loglike_4 <- -sum(log(do.call(pdistnam, 
                                              c(list(icens$right),  # interval censored 
                                                as.list(par), 
                                                as.list(fix.arg))) - do.call(pdistnam, 
                                                                             c(list(icens$left),  # interval censored 
                                                                               as.list(par),
                                                                               as.list(fix.arg)
                                                                               )
                                                                             )
                                      )
                                  )
                
                return(loglike_1 + loglike_2 + loglike_3 + loglike_4)
            }
        }
    }
    
    # case 3: uncensored with weights: ----
    else if (!cens && !is.null(weights)) {
        fnobj <- function(par, fix.arg, obs, ddistnam) {
            -sum(weights * log(do.call(ddistnam, c(list(obs), 
                                                   as.list(par), as.list(fix.arg)))))
        }
    }
    
    # case 4: censored with weights: ----
    else if (cens && !is.null(weights)) {
        fnobjcens <- function(par, fix.arg, rcens, lcens, icens, 
                              ncens, ddistnam, pdistnam) {
            p1 <- log(do.call(ddistnam, c(list(ncens), as.list(par), 
                                          as.list(fix.arg))))
            p2 <- log(do.call(pdistnam, c(list(lcens), as.list(par), 
                                          as.list(fix.arg))))
            p3 <- log(1 - do.call(pdistnam, c(list(rcens), as.list(par), 
                                              as.list(fix.arg))))
            p4 <- log(do.call(pdistnam, c(list(icens$right), 
                                          as.list(par), as.list(fix.arg))) - do.call(pdistnam, 
                                                                                     c(list(icens$left), as.list(par), as.list(fix.arg))))
            -sum(weights[irow.ncens] * p1) - sum(weights[irow.lcens] * 
                                                     p2) - sum(weights[irow.rcens] * p3) - sum(weights[irow.icens] * 
                                                                                                   p4)
        }
    }
    
    owarn <- getOption("warn")
    
    
    
    #**************************************************************************
    # start optimization steps ---- 
    if (is.null(custom.optim)) {
        
        hasbound <- any(is.finite(lower) | is.finite(upper))
        
        if (optim.method == "default") {
            meth <- ifelse(length(vstart) > 1, 
                           "Nelder-Mead",  
                           "BFGS")
        }
        else meth <- optim.method
        
        if (meth == "BFGS" && hasbound && is.null(gradient)) {
            meth <- "L-BFGS-B"
            txt1 <- "The BFGS method cannot be used with bounds without provided the gradient."
            txt2 <- "The method is changed to L-BFGS-B."
            warning(paste(txt1, txt2))
        }
        
        options(warn = ifelse(silent, -1, 0))
        if (hasbound) {
            if (!is.null(gradient)) {
                opt.fun <- "constrOptim"
            }
            else {
                if (meth == "Nelder-Mead") 
                    opt.fun <- "constrOptim"
                else if (meth %in% c("L-BFGS-B", "Brent")) 
                    opt.fun <- "optim"
                else {
                    txt1 <- paste("The method", meth, "cannot be used by constrOptim() nor optim() without gradient and bounds.")
                    txt2 <- "Only optimization methods L-BFGS-B, Brent and Nelder-Mead can be used in such case."
                    stop(paste(txt1, txt2))
                }
            }
            if (opt.fun == "constrOptim") {
                npar <- length(vstart)
                lower <- as.double(rep_len(lower, npar))
                upper <- as.double(rep_len(upper, npar))
                haslow <- is.finite(lower)
                Mat <- diag(npar)[haslow, ]
                hasupp <- is.finite(upper)
                Mat <- rbind(Mat, -diag(npar)[hasupp, ])
                colnames(Mat) <- names(vstart)
                rownames(Mat) <- paste0("constr", 1:NROW(Mat))
                Bnd <- c(lower[is.finite(lower)], -upper[is.finite(upper)])
                names(Bnd) <- paste0("constr", 1:length(Bnd))
                initconstr <- Mat %*% vstart - Bnd
                if (any(initconstr < 0)) 
                    stop("Starting values must be in the feasible region.")
                if (!cens) {
                    opttryerror <- try(opt <- constrOptim(theta = vstart, 
                                                          f = fnobj,
                                                          ui = Mat,
                                                          ci = Bnd,
                                                          grad = gradient, 
                                                          fix.arg = fix.arg, 
                                                          obs = data, 
                                                          ddistnam = ddistname, 
                                                          hessian = !is.null(gradient),
                                                          method = meth, 
                                                          ...),
                                       silent = TRUE)
                }
                else opttryerror <- try(opt <- constrOptim(theta = vstart, 
                                                           f = fnobjcens, ui = Mat, ci = Bnd, grad = gradient, 
                                                           ddistnam = ddistname, rcens = rcens, lcens = lcens, 
                                                           icens = icens, ncens = ncens, pdistnam = pdistname, 
                                                           fix.arg = fix.arg, hessian = !is.null(gradient), 
                                                           method = meth, ...), silent = TRUE)
                if (!inherits(opttryerror, "try-error")) 
                    if (length(opt$counts) == 1) 
                        opt$counts <- c(opt$counts, NA)
            }
            else {
                if (!cens) 
                    opttryerror <- try(opt <- optim(par = vstart, 
                                                    fn = fnobj, fix.arg = fix.arg, obs = data, 
                                                    gr = gradient, ddistnam = ddistname, hessian = TRUE, 
                                                    method = meth, lower = lower, upper = upper, 
                                                    ...), silent = TRUE)
                else opttryerror <- try(opt <- optim(par = vstart, 
                                                     fn = fnobjcens, fix.arg = fix.arg, gr = gradient, 
                                                     rcens = rcens, lcens = lcens, icens = icens, 
                                                     ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, 
                                                     hessian = TRUE, method = meth, lower = lower, 
                                                     upper = upper, ...), silent = TRUE)
            }
        }
        else {
            opt.fun <- "optim"
            if (!cens) 
                opttryerror <- try(opt <- optim(par = vstart, 
                                                fn = fnobj, fix.arg = fix.arg, obs = data, 
                                                gr = gradient, ddistnam = ddistname, hessian = TRUE, 
                                                method = meth, lower = lower, upper = upper, 
                                                ...), silent = TRUE)
            else opttryerror <- try(opt <- optim(par = vstart, 
                                                 fn = fnobjcens, fix.arg = fix.arg, gr = gradient, 
                                                 rcens = rcens, lcens = lcens, icens = icens, 
                                                 ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, 
                                                 hessian = TRUE, method = meth, lower = lower, 
                                                 upper = upper, ...), silent = TRUE)
        }
        options(warn = owarn)
        if (inherits(opttryerror, "try-error")) {
            warnings("The function optim encountered an error and stopped.")
            if (getOption("show.error.messages")) 
                print(attr(opttryerror, "condition"))
            return(list(estimate = rep(NA, length(vstart)), convergence = 100, 
                        loglik = NA, hessian = NA, optim.function = opt.fun, 
                        fix.arg = fix.arg, optim.method = meth, fix.arg.fun = fix.arg.fun, 
                        counts = c(NA, NA)))
        }
        if (opt$convergence > 0) {
            warnings("The function optim failed to converge, with the error code ", 
                     opt$convergence)
        }
        if (is.null(names(opt$par))) 
            names(opt$par) <- names(vstart)
        res <- list(estimate = opt$par, convergence = opt$convergence, 
                    value = opt$value, hessian = opt$hessian, optim.function = opt.fun, 
                    optim.method = meth, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                    weights = weights, counts = opt$counts, optim.message = opt$message, 
                    loglik = -opt$value)
    } else {
        
        options(warn = ifelse(silent, -1, 0))
        if (!cens) 
            opttryerror <- try(opt <- custom.optim(fn = fnobj, 
                                                   fix.arg = fix.arg, obs = data, ddistnam = ddistname, 
                                                   par = vstart, ...), silent = TRUE)
        else opttryerror <- try(opt <- custom.optim(fn = fnobjcens, 
                                                    fix.arg = fix.arg, rcens = rcens, lcens = lcens, 
                                                    icens = icens, ncens = ncens, ddistnam = ddistname, 
                                                    pdistnam = pdistname, par = vstart, ...), silent = TRUE)
        options(warn = owarn)
        if (inherits(opttryerror, "try-error")) {
            warnings("The customized optimization function encountered an error and stopped.")
            if (getOption("show.error.messages")) 
                print(attr(opttryerror, "condition"))
            return(list(estimate = rep(NA, length(vstart)), convergence = 100, 
                        loglik = NA, hessian = NA, optim.function = custom.optim, 
                        fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                        counts = c(NA, NA)))
        }
        if (opt$convergence > 0) {
            warnings("The customized optimization function failed to converge, with the error code ", 
                     opt$convergence)
        }
        if (is.null(names(opt$par))) 
            names(opt$par) <- names(vstart)
        argdot <- list(...)
        method.cust <- argdot$method
        res <- list(estimate = opt$par, convergence = opt$convergence, 
                    value = opt$value, hessian = opt$hessian, optim.function = custom.optim, 
                    optim.method = method.cust, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                    weights = weights, counts = opt$counts, optim.message = opt$message, 
                    loglik = -opt$value)
    }
    return(res)
}