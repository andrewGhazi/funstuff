library(tidyverse)
library(magrittr)
library(data.table)
library(rstanarm)
covid = data.table(s = c(5,90), # moderna I think?
                   h = c(14995,14910),
                   vac = c(1,0))

covid_glm = stan_glm(cbind(s,h) ~ vac, family = binomial(),
         data = covid)
inv_logit = function(x){ exp(x) / (1 + exp(x))}
1 - inv_logit(covid_glm$stan_summary['vac', c('2.5%', 'mean', '97.5%')])

pfizer = data.table(s = c(8,162),
                    h = c(22e4-8,22e4-162),
                    vac = c(1,0))
pfizer_glm = stan_glm(cbind(s,h) ~ vac, family = binomial(),
                     data = pfizer, 
                     cores = 4, chains = 4, iter = 4e5+1000, warmup = 1000)
inv_logit = function(x){ exp(x) / (1 + exp(x))}
1 - inv_logit(pfizer_glm$stan_summary['vac', c('2.5%', 'mean', '97.5%')])
