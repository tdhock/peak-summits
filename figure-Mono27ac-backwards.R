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
  chromEnd=as.integer(cumsum(bases)),
  orig.chromStart=chromStart,
  orig.chromEnd=chromEnd)]
all.rev.dt[, chromStart := as.integer(chromEnd-bases)]
out.rev.dt <- all.rev.dt[, .(
  chrom,
  chromStart,
  chromEnd,
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

start.join <- all.rev.dt[, .(chromStart, orig.chromEnd)]
end.join <- all.rev.dt[, .(chromEnd, orig.chromStart)]

pen.str <- "3000"
direction.vec <- list(
  forward=data.dir,
  reverse=rev.dir)
loss.dt.list <- list()
segs.dt.list <- list()
for(direction in names(direction.vec)){
  direction.dir <- direction.vec[[direction]]
  for(allow.free.changes in c(TRUE, FALSE)){
    fit <- problem.PeakSegFPOP(
      direction.dir, pen.str, allow.free.changes=allow.free.changes)
    fwd.segs <- if(direction=="reverse"){
      start.dt <- start.join[fit$segments, on="chromStart"]
      end.dt <- end.join[start.dt, on="chromEnd"]
      end.dt[.N:1, .(
        chrom,
        chromStart=orig.chromStart,
        chromEnd=orig.chromEnd,
        mean)]
    }else{
      fit$segments[, .(chrom, chromStart, chromEnd, mean)]
    }
    segs.dt.list[[paste(direction, allow.free.changes)]] <- data.table(
      direction, allow.free.changes, fwd.segs)
    loss.dt.list[[paste(direction, allow.free.changes)]] <- data.table(
      direction, allow.free.changes, fit$loss)
  }
}
segs.dt <- do.call(rbind, segs.dt.list)
loss.dt <- do.call(rbind, loss.dt.list)

loss.dt[order(allow.free.changes), .(
  allow.free.changes, direction, total.cost, penalty, peaks)]

all.diff.dt <- segs.dt[, {
  get1 <- function(d){
    dt <- .SD[direction==d]
    dt[, chromStart1 := chromStart+1L]
    setkey(dt, chromStart1, chromEnd)
    dt
  }
  fwd.dt <- get1("forward")
  rev.dt <- get1("reverse")
  over.dt <- foverlaps(fwd.dt, rev.dt)
  over.dt[, max.chromStart := ifelse(
    chromStart < i.chromStart, i.chromStart, chromStart)]
  over.dt[, min.chromEnd := ifelse(
    chromEnd < i.chromEnd, chromEnd, i.chromEnd)]
  diff.dt <- over.dt[order(max.chromStart)]
  diff.dt[, all.equal(min.chromEnd[-.N], max.chromStart[-1])]
  diff.dt[, max.mean := ifelse(mean < i.mean, i.mean, mean)]
  diff.dt[, min.mean := ifelse(mean < i.mean, mean, i.mean)]
  diff.dt
}, by=allow.free.changes]

segs.dt[, edges := ifelse(allow.free.changes, 4, 2)]
some.diff.dt <- all.diff.dt[min.mean<max.mean]
gg.zoomout <- ggplot()+
  ggtitle("Differences between directions highlighted with black circles")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(edges ~ ., labeller=label_both)+
  geom_step(aes(
    chromStart, count),
    color="grey50",
    data=Mono27ac$coverage)+
  geom_step(aes(
    chromStart, mean, color=direction, size=direction, linetype=direction),
    data=segs.dt)+
  scale_size_manual(values=c(forward=2, reverse=1))+
  geom_point(aes(
    max.chromStart, min.mean),
    size=5,
    shape=21,
    data=some.diff.dt)
png("figure-Mono27ac-backwards-zoomout.png", 12, 6, units="in", res=100)
print(gg.zoomout)
dev.off()


gg.middle <- ggplot()+
  ggtitle("Differences highlighted with black rectangles")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(edges ~ ., labeller=label_both)+
  geom_step(aes(
    chromStart, count),
    color="grey50",
    data=Mono27ac$coverage)+
  geom_step(aes(
    chromStart, mean, color=direction, size=direction, linetype=direction),
    data=segs.dt[edges==4])+
  scale_size_manual(values=c(forward=2, reverse=1))+
  geom_rect(aes(
    xmin=max.chromStart, xmax=min.chromEnd,
    ymin=min.mean, ymax=max.mean),
    size=1,
    color="black",
    data=some.diff.dt)+
  coord_cartesian(xlim=c(207000, 209000), ylim=c(25, 35))
png("figure-Mono27ac-backwards-middle.png", 12, 6, units="in", res=100)
print(gg.middle)
dev.off()

gg.end <- ggplot()+
  ggtitle("Differences between directions highlighted with grey rectangles")+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(edges ~ ., labeller=label_both)+
  geom_step(aes(
    chromStart, count),
    data=Mono27ac$coverage)+
  geom_rect(aes(
    xmin=max.chromStart, xmax=min.chromEnd,
    ymin=min.mean, ymax=max.mean),
    alpha=0.5,
    data=some.diff.dt)+
  geom_segment(aes(
    chromStart, mean,
    xend=chromEnd, yend=mean,
    color=direction, size=direction, linetype=direction),
    data=segs.dt[edges==4])+
  scale_size_manual(values=c(forward=1, reverse=0.75))+
  coord_cartesian(xlim=c(578500, 580000), ylim=c(0, 2))
png("figure-Mono27ac-backwards.png", 12, 6, units="in", res=100)
print(gg.end)
dev.off()
