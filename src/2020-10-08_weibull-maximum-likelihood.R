

#******************************************************************
# Maximum likelihood estimation of Weibull parameters 
# 2020-10-08
# Nayef 

#******************************************************************

library(stringr)

# Create functions: ---- 
wbl_2p_density <- function(x, shape, scale){
    y <- shape * scale * x^(shape-1)
    return(y)
}

wbl_2p_survival <- function(x, shape, scale){
    y <- scale * x^shape 
}

wbl_2p_loglik <- function(data_uncensored, data_censored, shape1, scale1){
    
    # first get result for all uncensored: 
    uncensored_loglik <- 0
    for (i in 1:length(data_uncensored)) {
        x <- data_uncensored[i]
        loglik <- log(wbl_2p_density(x, shape1, scale1))
        uncensored_loglik <- uncensored_loglik + loglik 
    }
    
    # next get result for all data: 
    all_loglik <- 0 
    combined_data <- c(data_uncensored, data_censored)
    for (i in 1:length(combined_data)) {
        x <- combined_data[i]
        loglik <- log(wbl_2p_survival(x, shape1, scale1))
        all_loglik <- all_loglik + loglik
    }
    
    # final logli: 
    final_loglik <- uncensored_loglik - all_loglik
    
}



# generate data without censoring: ----
data <- rweibull(50, 2, 1)

# try out lots of param guesses: 
param_guesses <- list(guess1 = c(shape = 1, scale = 1), 
                      guess2 = c(shape = 2, scale = 1))

for (i in 1:length(param_guesses)){
    shape_1 <- param_guesses[[i]][1]
    scale_1 <- param_guesses[[i]][2]
    print(str_glue("Shape: {shape_1}, scale: {scale_1}"))
    
    loglik <- wbl_2p_loglik(data_uncensored = data, 
                            data_censored = NULL,
                            shape1 = shape_1,
                            scale1 = scale_1)
    param_guesses[[i]][3] <- loglik
    names(param_guesses[[i]][3]) <- "loglik"

}

param_guesses
