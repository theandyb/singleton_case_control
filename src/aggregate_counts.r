library(tidyverse)
library(janitor)

args = commandArgs(trailingOnly=TRUE)

file_location <- args[1] #"/net/snowwhite/home/beckandy/research/singleton_case_control/data/chr22.csv"
out_file <- args[2]
# read file in
df <- read_csv(file_location, col_names = c("Cat", "CHR", "Motif", "Pos")) %>%
    mutate(Center = substr(Motif, 5, 5))

rc <- function(z){
  rc1 <- function(zz){
    s <- strsplit(zz, split = "")[[1]]
    s <- rev(s)
    dchars <- strsplit("ACGTMRWSYKVHDBNI", split = "")[[1]]
    comps <- strsplit("TGCAKYWSRMBDHVNI", split = "")[[1]]
    s <- s[s %in% dchars] # remove spaces etc
    s <- dchars[match(s, comps)]
    s <- paste0(s, collapse = "")
    return(s)
  }
  z <- toupper(z)
  tmpnames <- names(z)
  res <- unname(sapply(z, rc1))
  if(!is.null(attr(z, "quality"))){
    strev <- function(x) sapply(lapply(lapply(unname(x), charToRaw), rev), rawToChar)
    attr(res, "quality") <- unname(sapply(attr(z, "quality"), strev))
  }
  names(res) <- tmpnames
  return(res)
}

df1 <- df %>%
    filter(Center == "T" | Center == "G")

df2 <- df %>%
    filter(Center == "A" | Center == "C")

# Get reverse complement when center is T or G
df1 <- df1 %>%
    mutate(Motif = rc(Motif))

# Put it all back together
df <- bind_rows(df1, df2) %>% select(Cat, Motif)
rm(df1, df2)

# Tabulate some stuff!
## First: dummy vars for positions
#for(i in c(1:4, 6:9)){
#    df[[paste0("c",i)]] <- with(df,substr(Motif,i,i))
#}
#df <- df %>% select(Cat, c1:c9)

cats <- unique(df$Cat)

## Table for each category
tables <- list()
for(mt in cats){
    tables[[mt]] <- df %>%
        filter(Cat == mt) %>%
        tabyl(Motif)
}
saveRDS(tables, file = out_file)