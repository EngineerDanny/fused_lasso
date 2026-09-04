library(data.table)
library(ggplot2)

dataname <- "TwinsUK"

score.tall <- fread(
  paste0(
    "/projects/genomic-ml/da2343/necromass/",
    dataname,
    "_target_exempt_clr.score.tall.csv"
  )
)

## Put the four models on the same row for each test split and taxon.
score.wide <- dcast(
  score.tall,
  test.fold + test.subset + task_id ~ model,
  value.var = "regr.mse",
  fun.aggregate = mean
)

## Positive differences mean that fuser has lower prediction error.
score.wide[, `:=`(
  glmnet_all = glmnet_all_comp - fuser_all_comp,
  glmnet_same = glmnet_same_comp - fuser_all_comp
)]

taxon.differences <- score.wide[, .(
  `glmnet all - fuser` = mean(glmnet_all),
  `glmnet same - fuser` = mean(glmnet_same)
), by = task_id]

plot.dt <- melt(
  taxon.differences,
  id.vars = "task_id",
  variable.name = "comparison",
  value.name = "mse_difference"
)

plot.dt[, comparison := factor(
  comparison,
  levels = c(
    "glmnet all - fuser",
    "glmnet same - fuser"
  ),
  labels = c(
    "glmnet all minus fuser",
    "glmnet same minus fuser"
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
  scale_y_continuous(
    "MSE difference",
    expand = expansion(mult = c(0.08, 0.12))
  ) +
  scale_x_discrete(NULL) +
  labs(
    title = "Prediction performance after target exempt CLR",
    subtitle = paste0(
      dataname,
      ": each point is the mean test MSE difference for one taxon"
    ),
    caption = "Values above zero favor fuser."
  )

print(gg)

ggsave(
  paste0(dataname, "_target_exempt_clr_mse_difference.png"),
  gg,
  width = 7,
  height = 5,
  dpi = 500
)
