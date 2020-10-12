

#******************************************************************
# Maximum likelihood estimation of Weibull parameters 
# 2020-10-08
# Nayef 

#******************************************************************

library(tidyverse)
library(fitdistrplus)

# Create functions: ---- 
wbl_2p_neg_loglik <- function(df_data, shape, scale){
    par <- c(shape = shape, scale = scale)
    ncens <- df_data %>% filter(is_censored == 0) %>% pull(value)
    rcens <- df_data %>% filter(is_censored == 1) %>% pull(value)
    
    loglike_1 <- -sum(log(do.call(dweibull,  # pdf 
                                  c(list(ncens),  # uncensored 
                                    as.list(par)
                                    )
                                  )
                          )
                      )
    
    loglike_2 <- -sum(log(1 - do.call(pweibull,  # survival function 
                                      c(list(rcens), # right-censored
                                        as.list(par)
                                        )   
                                      )
                          )
                      )
    
    return(loglike_1 + loglike_2)
}



# generate data without censoring: ----
param_shape <- 2
param_scale <- 3
sample_n <- 100 

data <- rweibull(sample_n, param_shape, param_scale)
density(data) %>% plot

df0_data <- tibble(value = data, is_censored = 0)

# evaluate loglik with several param guesses: ------- 
# set up guesses: 
shape_grid <- seq(1, 2.5, length.out = sample_n)
scale_grid <- seq(1, 5, length.out = sample_n)

df1_loglik <- 
    expand.grid(shape = shape_grid,
                scale = scale_grid) %>% 
    as.tibble() 

# calculations: 
df1_loglik <- 
    df1_loglik %>% 
    mutate(neg_loglik_value = purrr::map2_dbl(shape, scale, 
                                          wbl_2p_neg_loglik, 
                                          df_data = df0_data))

# results: 
# df1_loglik
df1_loglik %>% arrange(neg_loglik_value)
wbl_2p_neg_loglik(df0_data, shape = param_shape, scale = param_scale)

# contour plot
df1_loglik %>% 
    ggplot(aes(x = shape,
               y = scale, 
               z = neg_loglik_value)) + 
    geom_contour(breaks=seq(0, 200, by=1)) + 
    
    geom_vline(xintercept = param_shape, col="red") +
    geom_hline(yintercept = param_scale, col="red") +
    
    scale_fill_gradient(low="blue", high="red") + 
    labs(title = "Minimizing the negative of the log-likelihood", 
         subtitle = "Red lines show true values")


# check: fitting with fitdistrplus: 
fit <- fitdist(df0_data$value, "weibull")
summary(fit)
