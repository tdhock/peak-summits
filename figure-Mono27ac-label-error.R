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

pen.vec <- list(
  ##"way too many peaks"="50",
  "Too many peaks"="500",
  "Zero label errors"="3000",
  "Too few peaks"="50000")
Mono27ac$labels[, chromStart1 := chromStart + 1L]
state.labels <- Mono27ac$labels[annotation %in% c("peakStart", "peakEnd")]
summit.labels <- Mono27ac$labels[annotation %in% c("noPeaks", "peaks")]
lab.err.list <- list()
fit.segs.list <- list()
fit.state.list <- list()
l <- function(pname){
  sprintf("Penalty=%s\n%s", pen.vec[pname], pname)
}
m <- function(pen.name){
  factor(l(pen.name), l(names(pen.vec)))
}
for(pen.name in names(pen.vec)){
  pen.str <- pen.vec[[pen.name]]
  fit <- problem.PeakSegFPOP(data.dir, pen.str, allow.free.changes=TRUE)
  sorted.segs <- fit$segments[order(chromStart)]
  sorted.segs[, diff.status := c(0, diff(status=="peak"))]
  sorted.segs[, peak.i := cumsum(diff.status>0)]
  sorted.segs[, diff.after := c(diff(mean), NA)]
  sorted.segs[, annotation := ifelse(status=="peak", "peakStart", "peakEnd")]
  state.dt <- sorted.segs[, {
    extreme <-.SD[which.max(chromEnd)]
    list(
      stateStart=min(chromStart),
      stateEnd=max(chromEnd),
      extremeMid=extreme[, as.integer((chromEnd+chromStart)/2)],
      extremeMean=extreme$mean
    )}, by=list(peak.i, annotation)]
  state.dt[, stateStart1 := stateStart + 1L]
  setkey(state.dt, annotation, stateStart1, stateEnd)
  setkey(state.labels, annotation, chromStart1, chromEnd)
  state.over <- foverlaps(state.labels, state.dt)
  state.err <- state.over[, {
    states <- sum(!is.na(peak.i))
    list(
      fp=as.integer(1 < states),
      fn=as.integer(0 == states)
    )}, by=list(chrom, chromStart, chromEnd, annotation)]
  summit.dt <- state.dt[annotation=="peakStart", list(extremeMid)]
  summit.dt[, extremeMid0 := extremeMid]
  setkey(summit.dt, extremeMid, extremeMid0)
  setkey(summit.labels, chromStart1, chromEnd)
  summit.over <- foverlaps(summit.labels, summit.dt)
  summit.err <- summit.over[, {
    summits <- sum(!is.na(extremeMid))
    if(annotation=="noPeaks"){
      list(fp=as.integer(0<summits), fn=0L)
    }else{#peaks
      list(fp=0L, fn=as.integer(0==summits))
    }
  }, by=list(chrom, chromStart, chromEnd, annotation)]
  pen.fac <- m(pen.name)
  lab.err.list[[pen.name]] <- data.table(
    pen.fac, rbind(state.err, summit.err))
  fit.state.list[[pen.name]] <- data.table(
    pen.fac, state.dt)
  fit.segs.list[[pen.name]] <- data.table(
    pen.fac, sorted.segs)
}
lab.err <- do.call(rbind, lab.err.list)
lab.err[, status := ifelse(
  0 < fp, "false positive", ifelse(
    0 < fn, "false negative", "correct"))]
fit.state <- do.call(rbind, fit.state.list)
fit.segs <- do.call(rbind, fit.segs.list)

ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
lab.min <- Mono27ac$labels[1, chromStart]
lab.max <- Mono27ac$labels[.N, chromEnd]
lab.min <- 3e5
lab.max <- 4e5
type.colors <- c(
  data="grey50",
  model="blue")
fit.lines <- fit.segs[, data.table(
  pos=as.integer(t(cbind(chromStart, chromEnd))),
  mean=as.numeric(t(cbind(mean, mean)))
), by=list(pen.fac)]
peak.y <- -1.3
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(pen.fac ~ .)+
  penaltyLearning::geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd,
    fill=annotation),
    color="grey",
    size=0.25,
    data=Mono27ac$labels)+
  penaltyLearning::geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd,
    linetype=status),
    fill=NA,
    data=lab.err)+
  scale_fill_manual("label", values=ann.colors)+
  geom_step(aes(
    chromStart, count, color=type),
    data=data.table(type="data", Mono27ac$coverage))+
  geom_line(aes(
    pos, mean, color=type),
    alpha=0.7,
    size=0.7,
    data=data.table(type="model", fit.lines))+
  ## geom_segment(aes(
  ##   chromStart, mean,
  ##   xend=chromEnd, yend=mean),
  ##   color="green",
  ##   size=1,
  ##   data=fit.segs)+
  ## geom_vline(aes(
  ##   xintercept=extremeMid),
  ##   data=fit.state[annotation=="peakStart"],
  ##   color="green",
  ##   linetype="dashed")+
  geom_point(aes(
    extremeMid, peak.y, color=type),
    data=data.table(type="model", fit.state[annotation=="peakStart"]),
    shape=1)+
  scale_linetype_manual(
    "error type",
    limits=c("correct", 
             "false negative",
             "false positive"),
    values=c(correct=0,
             "false negative"=3,
             "false positive"=1))+
  xlab("Position on chromosome 11")+
  scale_y_continuous(
    "Number of aligned DNA sequence reads
(measure of H3K27ac histone modification)",
breaks=seq(0, 10, by=2))+
  guides(color="none")+
  scale_color_manual(values=type.colors)+
  geom_text(aes(
    x,y,label=label, color=type, hjust=hjust),
    data=data.table(pen.fac=m("Too many peaks"), rbind(
      data.table(hjust=1, x=4e5, y=3, label="Noisy coverage data", type="data"),
      data.table(hjust=0.5, x=31e4, y=9, label="Mean model", type="model"),
      data.table(hjust=0, x=375000, y=peak.y, label="Predicted peaks", type="model"))))
gg.zoom <- gg+
  coord_cartesian(xlim=c(lab.min, lab.max), ylim=c(-2, 10))
##print(gg.zoom)
png("figure-Mono27ac-label-error-zoom.png", 10, 4, res=300, units="in")
print(gg.zoom)
dev.off()
## gg.zoom <- gg+
##   coord_cartesian(xlim=c(50e4, 51e4), ylim=c(-2, 45))+
##   guides(linetype="none", fill="none")
## png("figure-Mono27ac-label-error-zoom2.png", 10, 4, res=300, units="in")
## print(gg.zoom)
## dev.off()
gg.out <- gg+
  coord_cartesian(xlim=c(2e5, 5.8e5), ylim=c(-2, 45), expand=FALSE)
png("figure-Mono27ac-label-error.png", 10, 4, res=300, units="in")
print(gg.out)
dev.off()
##system("display figure-Mono27ac-label-error.png")
