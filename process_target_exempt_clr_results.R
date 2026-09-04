library(data.table)
library(mlr3)
library(mlr3resampling)
library(mlr3batchmark)
library(batchtools)

dataname <- "MovingPictures"
reg.dir <- paste0(dataname, "_target_exempt_clr")
reg <- batchtools::loadRegistry(reg.dir)

## Process results after the jobs finish.
batchtools::getStatus(reg = reg)
jobs.after <- batchtools::getJobTable(reg = reg)
table(jobs.after$error)
jobs.after[!is.na(error)]

ids <- jobs.after[is.na(error), job.id]
bmr <- mlr3batchmark::reduceResultsBatchmark(ids, reg = reg)
save(
  bmr,
  file = paste0(
    "/projects/genomic-ml/da2343/necromass/",
    dataname,
    "_target_exempt_clr-bmr.RData"
  )
)

score.dt <- mlr3resampling::score(bmr)
score.dt[, model := data.table::fcase(
  learner_id == "cv_glmnet_comp" & train.subsets == "same", "glmnet_same_comp",
  learner_id == "cv_glmnet_comp" & train.subsets == "all", "glmnet_all_comp",
  learner_id == "fuser_all_comp" & train.subsets == "all", "fuser_all_comp",
  learner_id == "featureless_comp" & train.subsets == "same", "featureless_same_comp"
)]
score.dt <- score.dt[!is.na(model)]
score.dt[, log_regr.mse := log10(regr.mse)]

fwrite(
  score.dt,
  paste0(
    "/projects/genomic-ml/da2343/necromass/",
    dataname,
    "_target_exempt_clr.score.tall.csv"
  )
)
