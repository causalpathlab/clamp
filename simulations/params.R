library(tidyverse)


####################################################################
## Parameters settings for synthetic data

nn <- 500
# pp <- c(1000, 3000, 5000)
pp <- c(200, 1000)
# h2 <- c(0.01, 0.1, 0.2, 0.4)
h2 <- c(1, 0.5, 0.2)
n_effect_vars <- c(5, 10)


param_comb <- as.data.frame(expand.grid(n = nn, p = pp, h2 = h2,
                                        n_effect_vars = n_effect_vars)) %>%
  mutate(ID = 1 : n()) %>%
  slice(rep(1 : n(), each = 100)) %>% # (repeat each simulation 100 times)
  mutate(seed = 1 : n()) %>%
  select(ID, seed, n, p, h2, n_effect_vars)


write.table(param_comb, file = "./simulations/params-synthetic.txt",
            row.names = F,col.names = F)
# save as `params-synthetic.txt` in the remote host.
