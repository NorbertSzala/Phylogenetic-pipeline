suppressPackageStartupMessages({
  library(ape)
  library(TreeDist)
  library(optparse)
})


option_list <- list(
  make_option("--tree1", type="character"),
  make_option("--tree2", type="character"),
  make_option("--out", type="character"),
  make_option("--metric", type="character", default="RF")
)
opt <- parse_args(OptionParser(option_list = option_list))


t1 <- unroot(read.tree(opt$tree1))
t2 <- unroot(read.tree(opt$tree2))

common <- intersect(t1$tip.label, t2$tip.label)
t1 <- drop.tip(t1, setdiff(t1$tip.label, common))
t2 <- drop.tip(t2, setdiff(t2$tip.label, common))

stopifnot(length(t1$tip.label) >= 3)

tip_map <- data.frame(
  letter = LETTERS[seq_along(t1$tip.label)],
  species = t1$tip.label
)


t1$tip.label <- seq_along(t1$tip.label)
t2$tip.label <- seq_along(t2$tip.label)


metric_fun <- switch(opt$metric,
  RF     = RobinsonFouldsMatching,
  InfoRF = InfoRobinsonFoulds,
  SPI    = SharedPhylogeneticInfo,
  JRF    = JaccardRobinsonFoulds,
  stop("Unknown metric")
)


png(opt$out, width = 2000, height = 900, res = 150)

VisualizeMatching(
  metric_fun,
  t1,
  t2,
  Plot = TreeDistPlot,
  matchZeros = FALSE
)
legend(
  "topleft",
  legend = paste(tip_map$letter, tip_map$species),
  cex = 0.6,
  bty = "n"
)

dev.off()
