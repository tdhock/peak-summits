source("packages.R")

## Compute several models and confront with labels.
data(Mono27ac)
data.dir <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  "chr11:60000-580000")
dir.create(data.dir, recursive=TRUE, showWarnings=FALSE)
write.table(
  Mono27ac$coverage, file.path(data.dir, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")

all.rev.dt <- Mono27ac$coverage[.N:1, .(
  chrom, count,
  bases=chromEnd-chromStart,
  chromEnd=cumsum(bases),
  orig.chromStart=chromStart,
  orig.chromEnd=chromEnd)]
out.rev.dt <- all.rev.dt[, .(
  chrom,
  chromStart=as.integer(chromEnd-bases),
  chromEnd=as.integer(chromEnd),
  count=as.integer(count))]
rev.dir <- file.path(
  tempfile(),
  "H3K27ac-H3K4me3_TDHAM_BP",
  "samples",
  "Mono1_H3K27ac",
  "S001YW_NCMLS",
  "problems",
  out.rev.dt[, sprintf(
    "%s:%d-%d",
    chrom[1],
    chromStart[1],
    chromEnd[.N])])
dir.create(rev.dir, recursive=TRUE, showWarnings=FALSE)
write.table(
  out.rev.dt, file.path(rev.dir, "coverage.bedGraph"),
  col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t")



pen.str <- "3000"
direction.vec <- list(
  forward=data.dir,
  reverse=rev.dir)
loss.dt.list <- list()
for(direction in names(direction.vec)){
  direction.dir <- direction.vec[[direction]]
  for(allow.free.changes in c(TRUE, FALSE)){
    fit <- problem.PeakSegFPOP(
      direction.dir, pen.str, allow.free.changes=allow.free.changes)
    loss.dt.list[[paste(direction, allow.free.changes)]] <- data.table(
      direction, allow.free.changes, fit$loss)
  }
}
loss.dt <- do.call(rbind, loss.dt.list)

loss.dt[order(allow.free.changes), .(
  allow.free.changes, direction, total.cost, penalty, peaks)]


