source('../R/sureLDA.R')

L <- readRDS("../data/simdata.rds")

set.seed(123)
surelda_run <- sureLDA(L$X, L$weight, L$ICD, L$NLP, L$HU, L$filter)
surelda_scores <- surelda_run$scores
surelda_probs <- surelda_run$probs
str(surelda_run)
