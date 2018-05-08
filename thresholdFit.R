library(rstan)
library(dplyr)
library(ggplot2)
library(reshape2)
library(matrixStats)
library(scales)
library(knitr)
library(tikzDevice)
opts_chunk$set(dev = 'pdf')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores(), tikzDefaultEngine = 'xetex')
theme_set(theme_light(base_size = 8))

## Set working directory according to where your scripts are located.
setwd('../DATA_TRAFFIC_STOPS')
stops <- read.csv('https://raw.githubusercontent.com/jkastelan/DAA2018_jlk635/master/Data/cleanTXstops.csv')


########################################################

plot_benchmark_test = function(obs) {
  races = as.character(levels(obs$driver_race))
  mx = max(obs$search_rate)
  df = obs %>% 
    filter(driver_race == 'White') %>% 
    right_join(obs %>% filter(driver_race != 'White'), by = 'police_department')
  
  ggplot(df) + geom_point(aes(x=search_rate.x, y=search_rate.y, size=num_stops.y), shape =1 ) +
    geom_abline(slope=1, intercept=0, linetype='dashed') +
    scale_y_continuous('Minority search rate\n', limits=c(0, mx), labels=percent, expand = c(0, 0)) +
    scale_x_continuous('\nWhite search rate', limits=c(0, mx), labels=percent, expand = c(0, 0)) +
    scale_size_area(max_size=15) +
    theme(legend.title = element_blank(),
          legend.background = element_rect(fill = 'transparent'),
          panel.margin.x=unit(1.5, "cm")) +
    guides(size=FALSE) + facet_grid( .~driver_race.y)
}

plot_benchmark_test(stops)


########################################################

plot_outcome_test = function(obs) {
  races = as.character(levels(obs$driver_race))
  mx = max(obs$hit_rate)
  df = obs %>% filter(driver_race == 'White') %>% 
    right_join(obs %>% filter(driver_race != 'White'), by='police_department')
  
  ggplot(df) + geom_point(aes(x=hit_rate.x, y=hit_rate.y, size=num_stops.y), shape = 1) +
    geom_abline(slope=1, intercept=0, linetype='dashed') +
    scale_y_continuous('Minority hit rate\n', limits=c(0, mx), labels=percent, expand = c(0, 0)) +
    scale_x_continuous('\nWhite hit rate', limits=c(0, mx), labels=percent, expand = c(0, 0)) +
    scale_size_area(max_size=15) +
    theme(legend.title = element_blank(),
          legend.background = element_rect(fill = 'transparent'),
          panel.margin.x=unit(1.5, "cm")) +
    guides(size=FALSE, color = FALSE) + facet_grid( .~driver_race.y)
}

plot_outcome_test(stops)

########################################################

stops %>%
  group_by(driver_race) %>%
  summarise(total_stops = sum(num_stops),
            total_searches = sum(num_searches),
            total_hits = sum(num_hits),
            search_rate = round(total_searches / total_stops * 100, 1),
            hit_rate = round(total_hits/ total_searches * 100, 1)) %>%
  setNames(c('Driver Race', 'Total Stops', 'Total Searches',
             'Total Hits', 'Search Rate', 'Hit Rate'))

########################################################


stan_data = with(stops, list(
   N = nrow(stops),
   D = length(unique(police_department)),
   R = length(unique(driver_race)),
   d = as.integer(police_department),
   r = as.integer(driver_race),
   n = num_stops,
   s = num_searches,
   h = num_hits))

model <- stan_model(file = 'threshold_test.stan')

fit <- sampling(
  model, data = stan_data, iter=8000,
  init = 'random', chains=5,
  cores=5, refresh=50, warmup = 2500,
  control = list(adapt_delta = 0.95,
                 max_treedepth = 12,
                 adapt_engaged = TRUE))

post = rstan::extract(fit)

########################################################


signal_to_p = function(x, phi, delta){
  #Checked. Converts x -> p. 
  p = phi * dnorm(x, delta, 1) / (phi * dnorm(x, delta, 1) + (1 - phi) * dnorm(x, 0, 1));
  return(p)
}

plot_department_thresholds = function(obs, post) {
  colors = c('blue', 'black', 'red')
  races = as.character(levels(obs$driver_race))
  obs$thresholds = colMeans(signal_to_p(post$t_i, post$phi, post$delta))
  mx = max(obs$thresholds)
  df = obs %>% filter(driver_race == 'White') %>%
    right_join(obs %>% filter(driver_race != 'White'), by = 'police_department')
  
  ggplot(df) + 
    geom_point(aes(x=thresholds.x, y=thresholds.y, size = num_stops.y), alpha=0.8, shape = 1) +
    geom_abline(slope=1, intercept=0, linetype='dashed') +
    scale_y_continuous('Minority threshold\n', limits=c(0,mx), labels=percent, expand=c(0, 0)) +
    scale_x_continuous('\nWhite threshold', limits=c(0,mx), labels=percent, expand=c(0, 0)) +
    scale_size_area(max_size=15) +
    theme(legend.position=c(0.0,1.0),
          legend.justification=c(0,1),
          legend.title = element_blank(),
          legend.background = element_rect(fill = 'transparent'),
          panel.margin.x=unit(1.5, "cm")) +
    scale_color_manual(values = colors[-1], labels=races[-1]) +
    guides(size=FALSE) + facet_grid(.~driver_race.y)
}



plot_department_thresholds(stops, post)

########################################################


stops = stops %>% 
    mutate(thresholds = colMeans(signal_to_p(post$t_i, post$phi, post$delta))) %>%
    group_by(police_department) %>%
    mutate(total_stops = sum(num_stops)) %>%
    ungroup()
  
  na_replace = function(x, r) ifelse(is.finite(x), x, r)
  
  accumrowMeans = function(M, i, w = rep(1, nrow(M)), imax = max(i)) {
    t(sapply(1:imax, function(j) (i == j)*na_replace(w/sum(w[i == j]),0))) %*% M
  }
  
  avg_thresh = accumrowMeans(t(signal_to_p(post$t_i, post$phi, post$delta)),
                             as.integer(stops$driver_race), stops$total_stops)
  
  data.frame(levels(stops$driver_race),
             sprintf('%.3f', rowMeans(avg_thresh)),
             apply(rowQuantiles(avg_thresh, probs = c(0.025, 0.975)), 1,
                   function(x) paste0('(', paste0(sprintf('%.3f',x), collapse = ', '), ')'))
             ) %>%
    setNames(c('Driver Race', 'Average Threshold', '95% Credible Interval'))


########################################################


search_rate_ppc <- function(obs, post, ylim = 0.03) {
  obs$pred_search_rate = colMeans(post$search_rate)
  ggplot(data=obs, aes(x=pred_search_rate, y=pred_search_rate-search_rate)) +
    geom_point(aes(size=num_stops, color=driver_race), alpha = 0.8) + 
    scale_size_area(max_size=10) +
    scale_x_continuous('\nPredicted search rate', labels=percent)+
    scale_y_continuous('Search rate prediction error\n', labels=percent, limits=c(-ylim, ylim)) +
    geom_abline(slope=0, intercept=0, linetype='dashed') +
    theme(legend.position=c(1.0,0),
          legend.justification=c(1,0),
          legend.title = element_blank(),
          legend.background = element_rect(fill = 'transparent')) +
    scale_color_manual(values=c('blue','black','red', 'green4', 'purple')) +
    guides(size=FALSE)
}

hit_rate_ppc <- function(obs, post, ylim = 0.2) {
  obs$pred_hit_rate = colMeans(post$hit_rate)
  ggplot(data=obs, aes(x=pred_hit_rate, y=hit_rate-pred_hit_rate)) +
    geom_point(aes(size=num_stops, color=driver_race), alpha=0.8) + 
    scale_size_area(max_size=10) +
    scale_x_continuous('\nPredicted hit rate', labels=percent) +
    scale_y_continuous('Hit rate prediction error\n', labels=percent, limits = c(-ylim, ylim)) +
    geom_abline(slope=0, intercept=0, linetype='dashed') +
    theme(legend.position=c(1.0,0),
          legend.justification=c(1,0), 
          legend.title = element_blank(),
          legend.background = element_rect(fill = 'transparent'))+
    scale_color_manual(values=c('blue','black','red', 'green4', 'purple')) +
    guides(size=FALSE)
}

search_rate_ppc(stops, post)

hit_rate_ppc(stops, post)



########################################################

s <- summary(fit)
Rhat <- s$summary[,'Rhat']
n_eff <- s$summary[,'n_eff']
summary(Rhat, na.rm=T)
summary(n_eff)

#########################################################

#Report the maximum Rhat value:

max(summary(fit)$summary[,'Rhat'], na.rm=T)   #where fit is stan fit object

#And the minimum number of effective samples:

min(summary(fit)$summary[,'n_eff'], na.rm=T)

