## understanding merge

df.pracMeta <- data.frame(sim_gen = c(200, 200, 200, 200), inv_id = c(1880, 2111, 11240, 9597), 
                          inv_length = c(40427, 1079, 33333, 34797), seed = c(rep(3383298, 2), rep(3383299, 2)))

df.pracTime <- data.frame(sim_gen = c(800, 1000, 1200, 1400, 800, 1000, 1200, 1400), 
                          inv_id = c(1880, 2111, 1880, 2111, 11240, 9597, 11240, 9597), 
                          freq = c(0.07, 0.2, 0.09, 0.4, 0.07, 0.2, 0.09, 0.4), 
                          seed = c(rep(3383298, 4), rep(3383299, 4)))

merge(df.pracMeta, df.pracTime, by.x = c("inv_id", "seed"), by.y = c("inv_id", "seed"), all.y = TRUE)
