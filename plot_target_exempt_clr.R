library(data.table)
library(ggplot2)

datanames <- c("MovingPictures", "TwinsUK", "necromass")
score.files <- paste0(
  "/projects/genomic-ml/da2343/necromass/",
  datanames,
  "_target_exempt_clr.score.tall.csv"
)
names(score.files) <- datanames
score.files <- score.files[file.exists(score.files)]

score.tall <- rbindlist(lapply(names(score.files), function(dataset.name) {
  score.dt <- fread(score.files[[dataset.name]])
  score.dt[, dataset := dataset.name]
  score.dt
}))

## Put the four models on the same row for each test split and taxon.
score.wide <- dcast(
  score.tall,
  dataset + test.fold + test.subset + task_id ~ model,
  value.var = "regr.mse",
  fun.aggregate = mean
)

## Negative differences mean that fuser has lower prediction error.
score.wide[, `:=`(
  glmnet_all = fuser_all_comp - glmnet_all_comp,
  glmnet_same = fuser_all_comp - glmnet_same_comp
)]

taxon.differences <- score.wide[, .(
  `fuser - glmnet all` = mean(glmnet_all),
  `fuser - glmnet same` = mean(glmnet_same)
), by = .(dataset, task_id)]

plot.dt <- melt(
  taxon.differences,
  id.vars = c("dataset", "task_id"),
  variable.name = "comparison",
  value.name = "mse_difference"
)

plot.dt[, comparison := factor(
  comparison,
  levels = c(
    "fuser - glmnet all",
    "fuser - glmnet same"
  ),
  labels = c(
    "glmnet all",
    "glmnet same"
  )
)]

gg <- ggplot(plot.dt, aes(comparison, mse_difference)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.10,
    height = 0
  ) +
  facet_wrap(~dataset, scales = "free_y") +
  scale_y_continuous(
    "MSE(fuser) - MSE(algorithm)",
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  scale_x_discrete(NULL) +
  labs(
    title = "Prediction performance after target exempt CLR",
    caption = "Values below zero favor fuser."
  )

ggsave(
  "target_exempt_clr_mse_difference.png",
  gg,
  width = 10,
  height = 4,
  dpi = 500
)
