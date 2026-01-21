#!/usr/bin/env Rscript

library(ape)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_rf_differences.R ref_tree.nwk query_tree.nwk output.png")
}

ref_file   <- args[1]
query_file <- args[2]
out_png    <- args[3]

# ----------------------------
# Load trees
# ----------------------------
ref   <- unroot(read.tree(ref_file))
query <- unroot(read.tree(query_file))

# remove gene IDs if present
ref$tip.label   <- sub("\\|.*$", "", ref$tip.label)
query$tip.label <- sub("\\|.*$", "", query$tip.label)

# keep only common taxa
common <- intersect(ref$tip.label, query$tip.label)
ref   <- drop.tip(ref, setdiff(ref$tip.label, common))
query <- drop.tip(query, setdiff(query$tip.label, common))

# ----------------------------
# Get splits (RF definition)
# ----------------------------
# prop.part() = bipartitions
ref_splits   <- prop.part(ref)
query_splits <- prop.part(query)

# convert splits to comparable strings
split_to_string <- function(split, n) {
  v <- rep(0, n)
  v[split] <- 1
  paste(v, collapse = "")
}

n <- length(query$tip.label)

ref_strings <- sapply(ref_splits, split_to_string, n = n)
qry_strings <- sapply(query_splits, split_to_string, n = n)

# splits unique to query tree = RF differences
diff_strings <- setdiff(qry_strings, ref_strings)

# ----------------------------
# Color edges in query tree
# ----------------------------
edge_cols <- rep("black", nrow(query$edge))

for (i in seq_along(query_splits)) {
  split_str <- split_to_string(query_splits[[i]], n)
  if (split_str %in% diff_strings) {
    node <- attr(query_splits, "node")[i]
    edge_cols[which(query$edge[,2] == node)] <- "red"
  }
}

# ----------------------------
# Plot
# ----------------------------
png(out_png, width = 1000, height = 900)

plot(
  query,
  edge.color = edge_cols,
  cex = 0.8,
  main = "Query tree\n(red = RF-unique splits)"
)

legend(
  "topleft",
  legend = c("shared split", "RF-unique split"),
  col = c("black", "red"),
  lwd = 2,
  bty = "n"
)

dev.off()
