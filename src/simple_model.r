library(tidyverse)

singletons <- readRDS("/net/snowwhite/home/beckandy/research/singleton_case_control/data/singleton_tables/chr22.Rds")
controls <- readRDS("/net/snowwhite/home/beckandy/research/singleton_case_control/data/control_tables/chr22.Rds")

mut_type <- "AT_GC"

final <- bind_rows((singletons[[mut_type]] %>% select(Motif, n) %>% mutate(mut = 1)),
                  (controls[[mut_type]] %>% select(Motif, n) %>% mutate(mut = 0)))

for(x in c(1:4,6:9)){
    final[[paste0("c", x)]] <- substr(final$Motif,x,x)
}

mod1 <- glm(n ~ (c1 + c2 + c3 + c4 + c6 + c7 + c8 + c9 +mut)^2, data = final, family = "poisson")