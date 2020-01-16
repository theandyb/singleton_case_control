library(tidyverse)

print("loading the data...")
control <- readRDS("/net/snowwhite/home/beckandy/research/singleton_case_control/data/control_tables/all_chrom.RData")
singleton <- readRDS("/net/snowwhite/home/beckandy/research/singleton_case_control/data/singleton_tables/all_chrom.RData")

get_df <- function(singleton, control, mut_type){
    if(stringr::str_sub(mut_type,1,2) == "GC"){
        df <- full_join(control[[paste0("cpg_", mut_type)]], control[[mut_type]], by = "Motif") %>%
            replace_na(list(n.x = 0, n.y = 0)) %>%
            mutate(n = n.x + n.y) %>%
            select(Motif, n) %>%
            mutate(mut = 0)
        df2 <- bind_rows(singleton[[paste0("cpg_", mut_type)]], singleton[[mut_type]]) %>%
            mutate(mut = 1)
        df <- bind_rows(df, df2)
        df$mut <- factor(df$mut, levels=c(0,1))
    } else {
        df <- control[[mut_type]] %>%
            mutate(mut = 0)
        df2 <- singleton[[mut_type]] %>%
            mutate(mut = 1)
        df <- bind_rows(df, df2)
    }
    
    for(i in c(1:10,12:21)){
        colName = paste0("c", i)
        df <- df %>%
            mutate(!!rlang::sym(colName) := stringr::str_sub(Motif, i, i))
    }
    
    return(df)
}

categories <- names(control)[!grepl("^cpg",names(control))]
for(mut_cat in categories){
    outFile <- paste0("/net/snowwhite/home/beckandy/research/singleton_case_control/data/by_cat/", mut_cat, ".Rdata")
    print(paste0("Getting data for ", mut_cat, "..."))
    df <- get_df(singleton, control, mut_cat)
    print(paste0("Writing to ", outFile, "..."))
    saveRDS(df, file = outFile)
    print("Done!")
}