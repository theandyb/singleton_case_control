library(tidyverse)
library(magrittr)

print("Loading the data")

##### PHASE UNO - Load the data
# Get genomewide counts
inputDir <- "/net/snowwhite/home/beckandy/research/smaug-redux/motif_counts/9-mers/full/"
gwCount <- read_tsv(paste0(inputDir, "chr", 1, ".9-mer_motifs_full.txt"), col_name= T, col_types = cols(
              CHR = col_double(),
              Motif = col_character(),
              nMotifs = col_double()
            )) %>%
    select(Motif, nMotifs)

for(chrom in 2:22){
    fileName <- paste0(inputDir, "chr", chrom, ".9-mer_motifs_full.txt")
    df <- read_tsv(fileName, col_name= T, col_types = cols(
              CHR = col_double(),
              Motif = col_character(),
              nMotifs = col_double()
            )) %>% 
        select(Motif, nMotifs)
    gwCount <- full_join(gwCount, df, by = "Motif") %>%
        replace_na(list(nMotifs.x = 0, nMotifs.y = 0)) %>%
        mutate(nMotifs = nMotifs.x + nMotifs.y) %>%
        select(Motif, nMotifs)
    rm(df)
}

gwCount <- gwCount %>%
    rename(n = nMotifs)

singleton <- read_tsv("/net/snowwhite/home/beckandy/research/smaug-redux/summaries/filtered/chr1.singletons",
                     col_names = T, 
                     col_types = cols(
                          CHR = col_double(),
                          POS = col_double(),
                          REF = col_character(),
                          ALT = col_character(),
                          AA = col_character(),
                          AN = col_double(),
                          Motif = col_character(),
                          Category = col_character()
                        )) %>% 
    group_by(Motif, Category) %>% summarize(nMut = n())

for(i in 2:22){
    singleton %<>% bind_rows(
        read_tsv(paste0("/net/snowwhite/home/beckandy/research/smaug-redux/summaries/filtered/chr", i, ".singletons"),
                  col_names = T, 
                     col_types = cols(
                          CHR = col_double(),
                          POS = col_double(),
                          REF = col_character(),
                          ALT = col_character(),
                          AA = col_character(),
                          AN = col_double(),
                          Motif = col_character(),
                          Category = col_character()
                        )) %>% 
    group_by(Motif, Category) %>% summarize(nMut = n()))
}

final <- full_join(singleton, gwCount, by = "Motif") %>%
    drop_na
rm(gwCount, singleton)

##### PHASE DOS - Subset and save the data

get_global_df <- function(final, mut_cat){
    if(grepl("GC_", mut_cat)){
        sub_global <- final %>% filter(grepl(mut_cat, Category, fixed=T)) %>% select(Motif, nMut, n)
        }
    else {
        sub_global <- final %>% filter(Category == mut_cat) %>% select(Motif, nMut, n)
    }

    for(i in c(1,2,3,4,6,7,8,9)){
            colName = paste0("c", i)
            sub_global <- sub_global %>%
                mutate(!!rlang::sym(colName) := stringr::str_sub(Motif, i, i))
    }
    
    sub_global$c1 <- as.factor(sub_global$c1)
    sub_global$c2 <- as.factor(sub_global$c2)
    sub_global$c3 <- as.factor(sub_global$c3)
    sub_global$c4 <- as.factor(sub_global$c4)
    sub_global$c6 <- as.factor(sub_global$c6)
    sub_global$c7 <- as.factor(sub_global$c7)
    sub_global$c8 <- as.factor(sub_global$c8)
    sub_global$c9 <- as.factor(sub_global$c9)
    
    return(sub_global %>% ungroup %>% select(-Motif))
}

mut_cats <- c("AT_TA", "AT_CG", "AT_GC", "GC_CG", "GC_AT", "GC_TA")
for(mc in mut_cats){
    print( paste0( "Category ", mc, "..." ))
    fOut <- paste0("/net/snowwhite/home/beckandy/research/singleton_case_control/data/global_train/", mc, ".csv")
    df <- get_global_df(final, mc)
    write_csv(df, fOut)
    print("Done!")
}