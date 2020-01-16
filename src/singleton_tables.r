library(tidyverse)
library(janitor)

args = commandArgs(trailingOnly=TRUE)

in_file <- args[1]
out_file <- args[2]

df <- read_tsv(in_file) %>% 
    select(Motif, Category) %>%
    rowwise %>% 
    mutate(Motif = substr(Motif, 1,21))

cats <- unique(df$Category)

tables <- list()
for(mt in cats){
    tables[[mt]] <- df %>%
        filter(Category == mt) %>%
        tabyl(Motif)
}
saveRDS(tables, file = out_file)