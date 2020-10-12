

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
sample_n <- 500 

data <- rweibull(sample_n, param_shape, param_scale)
density(data) %>% plot

df0_data <- tibble(value = data, is_censored = 0)

# evaluate loglik with several param guesses: ------- 
# set up guesses: 
df1_loglik <- 
    tibble(shape = seq(0.1, 10, length.out = sample_n),
           scale = runif(sample_n, 0, 100)) %>%
    
    # add the correct values: 
    bind_rows(data.frame(shape = param_shape,
                         scale = param_scale))
    
# calculations: 
df1_loglik <- 
    df1_loglik %>% 
    mutate(loglik_value = purrr::map2_dbl(shape, scale, 
                                          wbl_2p_loglik, 
                                          df_data = df0_data))

# results: 
# df1_loglik
df1_loglik %>% arrange(loglik_value)
wbl_2p_loglik(df0_data, shape = param_shape, scale = param_scale)

df1_loglik %>% 
    ggplot(aes(x = shape,
               y = scale, 
               z = loglik_value)) + 
    geom_contour() + 
    # facet_wrap(~as.factor(scale)) + 
    # geom_vline(xintercept = param_shape, col="red") + 
    labs(title = "Minimizing the negative of the log-likelihood", 
         subtitle = "Red line shows actual shape param value")


# check: fitting with fitdistrplus: 
fit <- fitdist(df0_data$value, "weibull")
summary(fit)
