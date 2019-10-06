library(tidyverse)

control_dir <- "/net/snowwhite/home/beckandy/research/singleton_case_control/data/control_tables/"
singleton_dir <- "/net/snowwhite/home/beckandy/research/singleton_case_control/data/singleton_tables/"

gw_control <- readRDS(paste0(control_dir,"chr1.RData"))
gw_single <- readRDS(paste0(singleton_dir,"chr1.RData"))
mut_types <- names(gw_control)

merge_counts <- function(x, y){
    x <- x %>% select(Motif, n)
    y <- y %>% select(Motif, n)
    z <- full_join(x, y, by = "Motif") %>%
        replace_na(list(n.x = 0, n.y = 0)) %>%
        mutate(n = n.x + n.y) %>%
        select(Motif, n)
    return(z)
}

for(i in 2:22){
    print(paste0("Chromosome ", i, "..."))
    control <- readRDS(paste0(control_dir,"chr",i,".RData"))
    singleton <- readRDS(paste0(singleton_dir,"chr",i,".RData"))
    for(mut_type in mut_types){
        gw_control[[mut_type]] <- merge_counts(control[[mut_type]], gw_control[[mut_type]])
        gw_single[[mut_type]] <- merge_counts(singleton[[mut_type]], gw_single[[mut_type]])
    }
}

saveRDS(gw_control, paste0(control_dir,"all_chrom.RData"))
saveRDS(gw_single, paste0(singleton_dir,"all_chrom.RData"))