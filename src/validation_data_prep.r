library(tidyverse)
library(magrittr)

#' Function for splitting motif into characters for the validation data
split_motif <- function(df){
    for(i in c(1,2,3,4,6,7,8,9)){
            colName = paste0("c", i)
            df %<>% mutate(!!rlang::sym(colName) := stringr::str_sub(SEQ, i, i))
    }
    
    df$c1 <- as.factor(df$c1)
    df$c2 <- as.factor(df$c2)
    df$c3 <- as.factor(df$c3)
    df$c4 <- as.factor(df$c4)
    df$c6 <- as.factor(df$c6)
    df$c7 <- as.factor(df$c7)
    df$c8 <- as.factor(df$c8)
    df$c9 <- as.factor(df$c9)
    
    return(df)
}

validation_file <- "/net/snowwhite/home/beckandy/research/smaug-redux/output/predicted/validation_sites_9mer.txt"
input_sites <- read_tsv(validation_file, 
                        col_names = c("CHR", "POS", "MU", "OBS", "Category", "SEQ", "ID"),
                        col_types = cols(
                            CHR = col_integer(),
                            POS = col_integer(),
                            MU = col_double(),
                            OBS = col_factor(levels = 0:1),
                            Category = col_factor(levels = c("AT_CG", "AT_GC", "AT_TA", "GC_AT", "GC_CG", "GC_TA")),
                            SEQ = col_character(),
                            ID = col_character()))
input_sites %<>% filter(MU > 0) %>%
    split_motif %>%
    select(-c(CHR,POS,SEQ))

# Do the subsampling
input_dnms <- input_sites %>% filter(ID != "all")
input_sites %<>% filter(ID == "all")

# Subsample the validation sites
set.seed(1337)
input_sites %<>% sample_n(1000000)

print("Saving...")
save(input_sites, input_dnms, file = "/net/snowwhite/home/beckandy/research/singleton_case_control/data/validation_100k.RData")