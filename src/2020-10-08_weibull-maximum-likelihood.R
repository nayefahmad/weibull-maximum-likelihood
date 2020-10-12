

#******************************************************************
# Maximum likelihood estimation of Weibull parameters 
# 2020-10-08
# Nayef 

#******************************************************************

library(tidyverse)
library(fitdistrplus)

# 0) Create functions: ---- 
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

# reference: https://statwonk.com/weibull.html
rweibull_cens <- function(n, shape, scale) {
    a_random_death_time <- rweibull(n, shape = shape, scale = scale) 
    a_random_censor_time <- rweibull(n, shape = shape, scale = scale)
    observed_time <- pmin(a_random_censor_time, a_random_death_time)
    censor <- as.integer(observed_time == a_random_death_time)
    tibble(value = observed_time, is_censored = censor)
}


#********************************************************************
# 1) data without censoring: ----
param_shape <- 2
param_scale <- 3
sample_n <- 100 
censoring <- "none"

data <- rweibull(sample_n, param_shape, param_scale)
density(data) %>% plot

df0_data <- tibble(value = data, is_censored = 0)

# 1.1) evaluate loglik with several param guesses: ------- 
# set up guesses: 
shape_grid <- seq(1, 5, length.out = sample_n)
scale_grid <- seq(1, 5, length.out = sample_n)

# note that we must pass a grid to geom_contour( ) below. Reference: https://vincenzocoia.com/post/contour_plots/
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

fitted_shape <- df1_loglik %>% arrange(neg_loglik_value) %>% slice(1) %>% pull(shape)
fitted_scale <- df1_loglik %>% arrange(neg_loglik_value) %>% slice(1) %>% pull(scale)


# 1.2) contour plot: ----
df1_loglik %>% 
    ggplot(aes(x = shape,
               y = scale, 
               z = neg_loglik_value)) + 
    geom_contour(breaks=seq(0, 200, by=1)) + 
    
    geom_vline(xintercept = param_shape, col="red") +
    geom_hline(yintercept = param_scale, col="red") +
    
    scale_fill_gradient(low="blue", high="red") + 
    labs(title = "Minimizing the negative log-likelihood of the two-parameter Weibull distribution", 
         subtitle = str_glue("Red lines show true values  \nCensoring: {censoring}  \nSample size: {sample_n}"))


# check: fitting with fitdistrplus: 
fit <- fitdist(df0_data$value, "weibull")
summary(fit)

# 1.3) data with fit: ---- 
df0_data %>% 
    ggplot(aes(x = value)) + 
    geom_density() + 
    stat_function(fun = dweibull, 
                  args = list(shape = fitted_shape, 
                              scale = fitted_scale), 
                  col = "dodgerblue2") + 
    stat_function(fun = dweibull, 
                  args = list(shape = param_shape, 
                              scale = param_scale), 
                  col = "red") + 
    labs(title = "Data and fitted distribution", 
         subtitle = "Black: data  \nRed: true distribution  \nBlue: fitted distribution")



#******************************************************************************
# 2) data with censoring: ----
param_shape <- 2
param_scale <- 3
sample_n <- 100 
censoring <- "right-censored"

df2_censdata <- rweibull_cens(sample_n, param_shape, param_scale)
df2_censdata$value %>% density() %>% plot()


# 2.1) evaluate loglik with several param guesses: ------- 
# set up guesses: 
shape_grid <- seq(1, 5, length.out = sample_n)
scale_grid <- seq(1, 5, length.out = sample_n)

# note that we must pass a grid to geom_contour( ) below. Reference: https://vincenzocoia.com/post/contour_plots/
df3_loglik_censdata <- 
    expand.grid(shape = shape_grid,
                scale = scale_grid) %>% 
    as.tibble() 

# calculations: 
df3_loglik_censdata <- 
    df3_loglik_censdata %>% 
    mutate(neg_loglik_value = purrr::map2_dbl(shape, scale, 
                                              wbl_2p_neg_loglik, 
                                              df_data = df2_censdata))

# results: 
# df1_loglik
df3_loglik_censdata %>% arrange(neg_loglik_value)
wbl_2p_neg_loglik(df2_censdata, shape = param_shape, scale = param_scale)

fitted_shape <- df3_loglik_censdata %>% arrange(neg_loglik_value) %>% slice(1) %>% pull(shape)
fitted_scale <- df3_loglik_censdata %>% arrange(neg_loglik_value) %>% slice(1) %>% pull(scale)


# 2.2) contour plot: -----
df3_loglik_censdata %>% 
    ggplot(aes(x = shape,
               y = scale, 
               z = neg_loglik_value)) + 
    geom_contour(breaks=seq(0, 200, by=1)) + 
    
    geom_vline(xintercept = param_shape, col="red") +
    geom_hline(yintercept = param_scale, col="red") +
    
    scale_fill_gradient(low="blue", high="red") + 
    labs(title = "Minimizing the negative log-likelihood of the two-parameter Weibull distribution", 
         subtitle = str_glue("Red lines show true values  \nCensoring: {censoring}  \nSample size: {sample_n}"))


# check: fitting with fitdistrplus: 
fit <- fitdist(df2_censdata$value, "weibull")
summary(fit)


# 2.3) data with fit: ---- 
df2_censdata %>% 
    ggplot(aes(x = value)) + 
    geom_density() + 
    stat_function(fun = dweibull, 
                  args = list(shape = fitted_shape, 
                              scale = fitted_scale), 
                  col = "dodgerblue2") +
    stat_function(fun = dweibull, 
                  args = list(shape = param_shape, 
                              scale = param_scale), 
                  col = "red") + 
    labs(title = "Data and fitted distribution", 
         subtitle = "Black: data  \nRed: true distribution  \nBlue: fitted distribution")
    
