library(data.table)
library(mlr3)
library(mlr3learners)
library(mlr3misc)
library(mlr3resampling)
library(R6)
library(paradox)
library(checkmate)
library(glmnet)
library(Matrix)

set.seed(42)

source("/projects/genomic-ml/da2343/necromass/l2_fusion.R")

dataname <- "TwinsUK"
task.dt <- data.table::fread(
  paste0("/projects/genomic-ml/da2343/ml_project_1/data/microbe_ds/",
         dataname, "_11_15.csv")
)
taxa.columns <- setdiff(names(task.dt), "Group_ID")
## These columns already contain log10(count + 1), as used in the paper.
log.count.mat <- as.matrix(task.dt[, ..taxa.columns])
group.id <- task.dt$Group_ID
task.list <- list()

for (out.i in seq_along(taxa.columns)) {
  target.name <- taxa.columns[[out.i]]
  predictor.index <- setdiff(seq_along(taxa.columns), out.i)

  ## The target is excluded from the center used for this regression.
  target.exempt.center <- rowMeans(log.count.mat[, predictor.index, drop = FALSE])
  target.response <- log.count.mat[, out.i] - target.exempt.center
  target.predictors <- sweep(
    log.count.mat[, predictor.index, drop = FALSE],
    MARGIN = 1,
    STATS = target.exempt.center,
    FUN = "-"
  )

  colnames(target.predictors) <- taxa.columns[predictor.index]
  task.dt <- data.table(target.predictors)
  task.dt[, (target.name) := target.response]
  task.dt[, Group_ID := group.id]

  reg.task <- mlr3::TaskRegr$new(
    id = target.name,
    backend = task.dt,
    target = target.name
  )
  reg.task$col_roles$subset <- "Group_ID"
  reg.task$col_roles$stratum <- "Group_ID"
  reg.task$col_roles$feature <- taxa.columns[predictor.index]
  task.list[[target.name]] <- reg.task
}

## Fused L2 learner used in the current analysis.
LearnerRegrFuserL2Comp <- R6Class(
  "LearnerRegrFuserL2Comp",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      ps <- ps(
        lambda = p_dbl(0, default = 0.01, tags = "train"),
        gamma = p_dbl(0, default = 0.1, tags = "train")
      )
      ps$values <- list(lambda = 0.01, gamma = 0.1)
      super$initialize(
        id = "regr.fuser_l2_comp",
        param_set = ps,
        feature_types = c("integer", "numeric"),
        label = "Fused L2 target exempt CLR",
        packages = c("glmnet", "Matrix")
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv <- self$param_set$get_values(tags = "train")
      subset.id <- task$col_roles$subset
      X.train <- as.matrix(task$data(cols = task$feature_names))
      y.train <- as.numeric(task$data(cols = task$target_names)[[1]])
      groups <- as.vector(task$data(cols = subset.id)[[1]])

      beta.estimate <- fusedL2DescentGLMNet(
        X = X.train,
        y = y.train,
        groups = groups,
        lambda = pv$lambda,
        gamma = pv$gamma,
        scaling = FALSE
      )

      self$model <- list(
        beta = beta.estimate,
        groups = groups,
        formula = task$formula(),
        data = task$data(),
        pv = pv
      )
      self$model
    },
    .predict = function(task) {
      subset.id <- task$col_roles$subset
      X.test <- as.matrix(task$data(cols = task$feature_names))
      new.groups <- as.vector(task$data(cols = subset.id)[[1]])
      y.predict <- predictFusedL2(
        beta.mat = self$model$beta,
        newX = X.test,
        newGroups = new.groups,
        groups = self$model$groups
      )
      list(response = as.numeric(y.predict))
    }
  )
)

## Learners and nested tuning.
grid.search <- mlr3tuning::TunerBatchGridSearch$new()
grid.search$param_set$values$resolution <- 5

subtrain.valid.cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
subtrain.valid.cv$param_set$values$folds <- 2
subtrain.valid.cv$param_set$values$subsets <- "A"

fuser.learner <- LearnerRegrFuserL2Comp$new()
fuser.learner$param_set$values$lambda <- paradox::to_tune(
  levels = c(0.001, 0.01, 0.1, 1)
)
fuser.learner$param_set$values$gamma <- paradox::to_tune(
  levels = c(0.001, 0.01, 0.1, 1)
)
fuser.learner.tuned <- mlr3tuning::auto_tuner(
  tuner = grid.search,
  learner = fuser.learner,
  resampling = subtrain.valid.cv,
  measure = mlr3::msr("regr.mse")
)
fuser.learner.tuned$id <- "fuser_all_comp"
fuser.learner.tuned$encapsulate(
  method = "evaluate",
  fallback = mlr3::LearnerRegrFeatureless$new()
)

glmnet.learner <- mlr3learners::LearnerRegrCVGlmnet$new()
glmnet.learner$id <- "cv_glmnet_comp"
glmnet.learner$param_set$values$alpha <- 1
glmnet.learner$param_set$values$nfolds <- 3
glmnet.learner$param_set$values$grouped <- TRUE
glmnet.learner$encapsulate(
  method = "evaluate",
  fallback = mlr3::LearnerRegrFeatureless$new()
)

featureless.learner <- mlr3::LearnerRegrFeatureless$new()
featureless.learner$id <- "featureless_comp"

## Same and All use the existing SAC resampling implementation.
same.cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
same.cv$param_set$values$folds <- 5
same.cv$param_set$values$subsets <- "S"

all.cv <- mlr3resampling::ResamplingSameOtherSizesCV$new()
all.cv$param_set$values$folds <- 5
all.cv$param_set$values$subsets <- "A"

same.grid <- mlr3::benchmark_grid(
  tasks = task.list,
  learners = list(glmnet.learner, featureless.learner),
  resamplings = same.cv
)
all.grid <- mlr3::benchmark_grid(
  tasks = task.list,
  learners = list(glmnet.learner, fuser.learner.tuned),
  resamplings = all.cv
)
reg.bench.grid <- data.table::rbindlist(
  list(same.grid, all.grid),
  use.names = TRUE,
  fill = TRUE
)

## Submit the benchmark jobs.
reg.dir <- paste0(dataname, "_target_exempt_clr")
unlink(reg.dir, recursive = TRUE)
reg <- batchtools::makeExperimentRegistry(
  file.dir = reg.dir,
  seed = 1,
  source = "/projects/genomic-ml/da2343/necromass/l2_fusion.R",
  packages = c(
    "data.table", "mlr3", "mlr3learners", "mlr3misc",
    "mlr3resampling", "mlr3tuning", "mlr3batchmark",
    "R6", "paradox", "checkmate", "glmnet", "Matrix"
  )
)

mlr3batchmark::batchmark(
  reg.bench.grid,
  store_models = TRUE,
  reg = reg
)
job.table <- batchtools::getJobTable(reg = reg)
chunks <- data.frame(job.table, chunk = 1)
batchtools::submitJobs(
  chunks,
  resources = list(
    walltime = 60 * 480,
    memory = 1024,
    ncpus = 1,
    ntasks = 1,
    chunks.as.arrayjobs = TRUE
  ),
  reg = reg
)
