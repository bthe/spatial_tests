library(spatioTemporalIndices)
# source("prep_ldists_for_setupData.R")
# out <- prep_ldists_for_setupData(
#   input_path = "ldists.csv",
#   output_path = "ldists_setupData.rds",
#   species = NA_integer_,
#   length_interval = 5,
#   min_length = 5,
#   max_length = 120
# )

conf <- defConf(years = 1985:2025, lengthGroups = seq(5, 100, by = 5))

source(
  "spatioTemporalIndices/R/utils.R"
)
source('tidy_prep.R')

#out <- readRDS("ldists_setupData.rds")

dat <- setupData(data = df, conf = conf)


par <- setPar(dat, conf)
map <- setMap(par, conf)


start_time <- Sys.time()
run = fitModel(data = dat, par = par, conf = conf, map = map)
end_time <- Sys.time()
timeUsed = end_time - start_time
timeUsed
