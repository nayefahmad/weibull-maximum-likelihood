

#******************************************************************
# Maximum likelihood estimation of Weibull parameters 
# 2020-10-08
# Nayef 

#******************************************************************

library(tidyverse)
library(fitdistrplus)

# Create functions: ---- 
wbl_2p_hazard <- function(x, shape, scale){
    y <- shape * scale * x^(shape-1)
    return(y)
}

wbl_2p_chf <- function(x, shape, scale){
    y <- scale * x^shape 
    return(y)
}

wbl_2p_loglik <- function(df_data, shape, scale){
    data_uncensored <- df_data %>% filter(is_censored == 0) %>% pull(value)
    all_data <- df_data %>% pull(value)
    
    # first get result for all uncensored: 
    log_hazard <- 0
    for (i in 1:length(data_uncensored)) {
        x <- data_uncensored[i]
        log_hazard_new <- log(wbl_2p_hazard(x, shape, scale))
        log_hazard <- log_hazard + log_hazard_new 
    }
    
    # next get result for all data: 
    chf <- 0 
    # combined_data <- c(data_uncensored, data_censored)
    for (i in 1:length(all_data)) {
        x <- all_data[i]
        chf_new <- wbl_2p_chf(x, shape, scale)
        chf <- chf + chf_new
    }
    
    # final loglik: 
    final_loglik <- log_hazard - chf 
    return(final_loglik)
}



# generate data without censoring: ----
param_shape <- 2
param_scale <- 3

data <- rweibull(500, param_shape, param_scale)
density(data) %>% plot

df0_data <- tibble(value = data, is_censored = 0)

# evaluate loglik with several param guesses: ------- 
# set up guesses: 
df1_loglik <- 
    tibble(shape = seq(0.1, 3, length.out = 500),
           scale = rep(seq(2.99, 3.02, length.out = 5), 100)) %>%
    
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
df1_loglik %>% arrange(desc(loglik_value))
wbl_2p_loglik(df0_data, shape = param_shape, scale = param_scale)

df1_loglik %>% 
    ggplot(aes(x = shape, y = loglik_value)) + 
    geom_point() + 
    facet_wrap(~as.factor(scale))


# check: fitting with fitdistrplus: 
fit <- fitdist(df0_data$value, "weibull")
summary(fit)
