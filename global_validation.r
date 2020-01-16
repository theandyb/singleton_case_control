#' global_validation.r
#' Extract rate estimates from poisson models
#' Fit likelihood for validation data set

library(tidyverse)
library(magrittr)
library(gtools)
library(rlang)

# Some function definitions for yinz

#' Function for loading the global training data (mutation and genomewide counts)
get_global_train <- function(mType){
    nucl <- c("A", "T", "C", "G")
    fName <- paste0("/net/snowwhite/home/beckandy/research/singleton_case_control/data/global_train/", mType, ".csv")
    df <- read_csv(fName, col_types = cols(nMut = col_integer(),
                                          n = col_integer(),
                                          c1 = col_factor(levels = nucl),
                                          c2 = col_factor(levels = nucl),
                                          c3 = col_factor(levels = nucl),
                                          c4 = col_factor(levels = nucl),
                                          c6 = col_factor(levels = nucl),
                                          c7 = col_factor(levels = nucl),
                                          c8 = col_factor(levels = nucl),
                                          c9 = col_factor(levels = nucl)))
    return(df)
}

#' Function for splitting motif into characters for the validation data
# split_motif <- function(df){
#     for(i in c(1,2,3,4,6,7,8,9)){
#             colName = paste0("c", i)
#             df %<>% mutate(!!rlang::sym(colName) := stringr::str_sub(SEQ, i, i))
#     }
#     return(df)
# }

#' Generate df to extract rate estimates from models
#'
#'
rate_df <- function(nSites, siteNames){
    df <- data.frame(permutations(n = 4, r = nSites, v = c("A","T","C","G"), repeats.allowed=T))
    names(df) <- siteNames
    df$n <- 1
    return(df)
}

# Also load the validation sites
# data.frames: input_sites and input_dnms
load("/net/snowwhite/home/beckandy/research/singleton_case_control/data/validation_100k.RData")


categories <- unique(input_sites$Category)
########################################
#           Fitting models             #
########################################

# First: three-mer with and without interaction

final_val_df <- data.frame()

for(mType in categories){
    print(mType)
    df <- get_global_train(mType)
    val_df <- input_sites %>% filter(Category == mType)
    val_dnm_df <- input_dnms %>% filter(Category == mType)
    val_df <- bind_rows(val_df, val_dnm_df)
    
    mod1 <- glm(nMut ~ c4 + c6 + offset(log(n)), data = df, family = "poisson")
    mod2 <- glm(nMut ~ c4 + c6 + c4:c6 + offset(log(n)), data = df, family = "poisson")
    rates <- rate_df(2, c("c4","c6"))
    rates %<>% mutate(r1 = exp(predict(mod1, newdata = .)))
    rates %<>% mutate(r2 = exp(predict(mod2, newdata = .)), resid1 = r2 - r1)
    
    final_val_df %<>% bind_rows(left_join(val_df, rates, by = c("c4", "c6")))
}

mod_val1 <- glm(OBS ~ r1, data = final_val_df, family=binomial())
mod_val2 <- glm(OBS ~ r1 + resid1, data = final_val_df, family=binomial())

lmtest::lrtest(mod_val1, mod_val2)

model_comparison <- function(form1, form2, input_sites, input_dnms, final, nSites, siteNames){
    final_val_df <- data.frame()
    categories <- unique(input_sites$Category)

    print("Fitting rate models...")
    for(mType in categories){
        print(mType)
        df <- get_global_train(mType)
        val_df <- input_sites %>% filter(Category == mType)
        val_dnm_df <- input_dnms %>% filter(Category == mType)
        val_df <- bind_rows(val_df, val_dnm_df)

        mod1 <- glm(as.formula(form1), data = df, family = "poisson")
        mod2 <- glm(as.formula(form2), data = df, family = "poisson")
        rates <- rate_df(nSites, siteNames)
        rates %<>% mutate(r1 = exp(predict(mod1, newdata = .)))
        rates %<>% mutate(r2 = exp(predict(mod2, newdata = .)), resid1 = r2 - r1)

        final_val_df %<>% bind_rows(left_join(val_df, rates, by = siteNames))
    }
    
    print("Fitting validation models...")
    mod_val1 <- glm(OBS ~ r1, data = final_val_df, family=binomial())
    mod_val2 <- glm(OBS ~ r1 + resid1, data = final_val_df, family=binomial())
    print(lmtest::lrtest(mod_val1, mod_val2))
}

# Example: c4+c6+c4:c6 versus c3 + c4 + c6 + c7 + c4:c6
model_comparison("nMut ~ c4 + c6 + c4:c6 + offset(log(n))", 
                 "nMut ~ c3 + c4 + c6 + c7 + c4:c6 + offset(log(n))", 
                 input_sites, input_dnms, final, 4, c("c3","c4","c6","c7"))

# And the logical next step...
model_comparison("nMut ~ c3 + c4 + c6 + c7 + c4:c6 + offset(log(n))", 
                 "nMut ~ (c3 + c4 + c6 + c7)^2 + offset(log(n))", 
                 input_sites, input_dnms, final, 4, c("c3","c4","c6","c7"))

model_comparison("nMut ~ (c3 + c4 + c6 + c7)^2 + c4:c6 + offset(log(n))", 
                 "nMut ~ (c3 + c4 + c6 + c7)^3 + offset(log(n))", 
                 input_sites, input_dnms, final, 4, c("c3","c4","c6","c7"))

model_comparison("nMut ~ (c3 + c4 + c6 + c7)^3 + c4:c6 + offset(log(n))", 
                 "nMut ~ (c3 + c4 + c6 + c7)^4 + offset(log(n))", 
                 input_sites, input_dnms, final, 4, c("c3","c4","c6","c7"))

# DON'T ACTUALLY RUN THIS ONE!
# Example: 7-mer with or without +/- 4
#model_comparison("nMut ~ (c2+c3+c4+c6+c7+c8)^6 + offset(log(n))", 
#                 "nMut ~  (c2+c3+c4+c6+c7+c8)^6 + c1 + c9 + offset(log(n))", 
#                 input_sites, input_dnms, final, 8, c("c1","c2","c3","c4","c6","c7","c8","c9"))

# Timing experiment
library(biglm)

mType <- "GC_TA"
df <- get_global_train(mType)
val_df <- input_sites %>% filter(Category == mType)
val_dnm_df <- input_dnms %>% filter(Category == mType)
val_df <- bind_rows(val_df, val_dnm_df)

system.time(
    mod1 <- glm(nMut ~ c4 + c6 + offset(log(n)), data = df, family = "poisson"))

system.time(
    mod1 <- bigglm(nMut ~ c4 + c6 + offset(log(n)), data = df, family = poisson()))

# Attempt to exploit sparsity
library(Matrix)
library(speedglm)
x <- sparse.model.matrix(~ (c3 + c4 + c6 + c7)^2 + offset(log(n)), df)
system.time(mod1 <- speedglm.wfit(df$nMut, x, intercept = FALSE, family = poisson()))