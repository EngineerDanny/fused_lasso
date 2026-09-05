library(data.table)
library(mlr3)
library(mlr3resampling)
library(mlr3batchmark)
library(batchtools)

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
dataname <- "TwinsUK"
reg.dir <- paste0(dataname, "_target_exempt_clr")
parts.dir <- paste0(reg.dir, "_scores")
manifest.file <- file.path(parts.dir, "manifest.rds")

if (mode == "prepare") {
  reg <- loadRegistry(reg.dir)
  jobs <- getJobTable(reg = reg)
  stopifnot(nrow(jobs) > 0, nrow(findNotDone(reg = reg)) == 0)
  # A job.name contains a complete task/learner/resampling configuration.
  manifest <- split(jobs$job.id, jobs$job.name)
  dir.create(parts.dir, showWarnings = FALSE)
  if (file.exists(manifest.file)) {
    stopifnot(identical(readRDS(manifest.file), manifest))
  } else {
    saveRDS(manifest, manifest.file)
  }
  cat(length(manifest), "collection tasks. Submit with:\n")
  cat(sprintf("sbatch --array=1-%d collect_clr_array.sbatch\n", length(manifest)))
} else if (mode == "collect") {
  manifest <- readRDS(manifest.file)
  index <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", args[2]))
  stopifnot(length(index) == 1, !is.na(index), index >= 1, index <= length(manifest))
  output <- file.path(parts.dir, sprintf("part-%04d.rds", index))
  if (file.exists(output)) {
    saved <- readRDS(output)
    stopifnot(identical(saved$ids, manifest[[index]]))
    cat("Already collected", index, "\n")
  } else {
    started <- proc.time()[[3]]
    reg <- loadRegistry(reg.dir)
    cat("Collecting", index, "of", length(manifest), ":", length(manifest[[index]]), "jobs\n")
    bmr <- reduceResultsBatchmark(manifest[[index]], reg = reg,
                                 store_backends = FALSE, unmarshal = FALSE)
    scores <- mlr3resampling::score(bmr)
    scores[, model := fcase(
      learner_id == "cv_glmnet_comp" & train.subsets == "same", "glmnet_same_comp",
      learner_id == "cv_glmnet_comp" & train.subsets == "all", "glmnet_all_comp",
      learner_id == "fuser_all_comp" & train.subsets == "all", "fuser_all_comp",
      learner_id == "featureless_comp" & train.subsets == "same", "featureless_same_comp"
    )]
    stopifnot(nrow(scores) == length(manifest[[index]]), !anyNA(scores$model))
    scores <- scores[, .(task_id, learner_id, train.subsets, test.fold,
                         test.subset, regr.mse, model, log_regr.mse = log10(regr.mse))]
    seconds <- proc.time()[[3]] - started
    # Rename only after a complete write; an interrupted task can be retried.
    temporary <- paste0(output, ".", Sys.getpid(), ".tmp")
    saveRDS(list(ids = manifest[[index]], scores = scores, seconds = seconds), temporary)
    stopifnot(file.rename(temporary, output))
    cat("Saved", output, "in", round(seconds, 1), "seconds\n")
  }
} else if (mode == "combine") {
  manifest <- readRDS(manifest.file)
  files <- file.path(parts.dir, sprintf("part-%04d.rds", seq_along(manifest)))
  if (!all(file.exists(files))) stop("Some collection tasks are unfinished. Rerun the array to resume.")
  parts <- lapply(files, readRDS)
  stopifnot(all(vapply(seq_along(parts), function(i) identical(parts[[i]]$ids, manifest[[i]]), logical(1))))
  scores <- rbindlist(lapply(parts, `[[`, "scores"))
  stopifnot(nrow(scores) == sum(lengths(manifest)))
  output <- paste0(dataname, "_target_exempt_clr.score.tall.csv")
  temporary <- paste0(output, ".tmp")
  fwrite(scores, temporary)
  stopifnot(file.rename(temporary, output))
  cat("Saved", nrow(scores), "scores to", output, "\n")
  cat("Median collection task time:", median(vapply(parts, `[[`, numeric(1), "seconds")), "seconds\n")
} else {
  stop("Use: Rscript collect_clr_array.R prepare | collect INDEX | combine")
}
