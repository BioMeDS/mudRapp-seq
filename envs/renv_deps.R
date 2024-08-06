# contains calls to dependencies that are otherwise not automatically inferred
# e.g. svglite as dependency of ggplot2, only if ggsave is used with svg output
# see https://github.com/rstudio/renv/issues/1930
library(svglite)
library(Hmisc)
