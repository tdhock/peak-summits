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

pen.vec <- c("5000", "3500")
Mono27ac$labels[, chromStart1 := chromStart + 1L]
lab <- function(chromStart, chromEnd, annotation){
  data.table(chrom="chr11", chromStart, chromEnd, annotation)
}
my.labels <- rbind(
  lab(175000, 2e5, "noPeaks"),
  lab(3e5, 32e4, "peaks"))
my.labels[, chromStart1 := chromStart + 1L]
state.labels <- my.labels[annotation %in% c("peakStart", "peakEnd")]
summit.labels <- my.labels[annotation %in% c("noPeaks", "peaks")]
lab.err.list <- list()
fit.segs.list <- list()
fit.state.list <- list()
for(pen.i in seq_along(pen.vec)){
  pen.str <- pen.vec[[pen.i]]
  fpop <- problem.PeakSegFPOP(
    data.dir, pen.str,
    allow.free.changes=TRUE)
  state.dt <- fpop$states
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
  l <- function(pname){
    sprintf("penalty=%s\n%s", pen.vec[pname], pname)
  }
  ##pen.fac <- factor(l(pen.name), l(names(pen.vec)))
  peaks.int <- fpop$loss$peaks
  lab.err.list[[paste(peaks.int)]] <- data.table(
    peaks.int, rbind(state.err, summit.err))
  fit.state.list[[paste(peaks.int)]] <- data.table(
    peaks.int, state.dt)
  fit.segs.list[[paste(peaks.int)]] <- data.table(
    peaks.int, fpop$segments)
}
lab.err <- do.call(rbind, lab.err.list)
lab.err[, status := ifelse(
  0 < fp, "false positive", ifelse(
    0 < fn, "false negative", "correct"))]
fit.state <- do.call(rbind, fit.state.list)
fit.segs <- do.call(rbind, fit.segs.list)

gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(peaks.int ~ .)+
  geom_step(aes(
    chromStart, count),
    color="grey50",
    data=Mono27ac$coverage)+
  geom_segment(aes(
    chromStart, mean,
    xend=chromEnd, yend=mean),
    color="green",
    size=1,
    data=fit.segs)+
  geom_point(aes(
    extremeMid, extremeMean),
    data=fit.state[annotation=="peakStart"],
    shape=21,
    color="black",
    fill="green")+
  xlab("Position on chromosome")+
  ylab("Number of aligned reads")
print(gg)

if(file.exists("Mono27acPeaks.RData")){
  load("Mono27acPeaks.RData")
}else{
  pdpa <- PeakSegOptimal::PeakSegPDPAchrom(Mono27ac$coverage, 13L)
  save(pdpa, file="Mono27acPeaks.RData")
}

show.segs <- data.table(pdpa$segments)[peaks %in% 8:12]
show.segs <- data.table(pdpa$segments)[peaks %in% fit.segs$peaks.int]
type.colors <- c(
  data="grey50",
  model="blue")
peak.y <- -4
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(peaks ~ .)+
  geom_step(aes(
    chromStart, count, color=type),
    data=data.table(type="data", Mono27ac$coverage))+
  geom_segment(aes(
    chromStart, mean,
    color=type,
    xend=chromEnd, yend=mean),
    size=1,
    data=data.table(type="model", show.segs))+
  geom_point(aes(
    chromStart, peak.y, color=type),
    shape=1,
    data=data.table(type="model", show.segs[status=="peak"]))+
  xlab("Position on chromosome")+
  ylab("Number of aligned reads")+
  scale_color_manual(values=type.colors)
print(gg)

ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
P <- 10
pdpa.segs <- show.segs[peaks==P]
pdpa.peaks <- pdpa.segs[status=="peak"]
getLines <- function(dt){
  dt[, data.table(
    pos=as.integer(t(cbind(chromStart, chromEnd))),
    mean=as.numeric(t(cbind(mean, mean))))]
}
pdpa.lines <- getLines(pdpa.segs)
fpop.segs <- fit.segs[peaks.int==P]
fpop.lines <- getLines(fpop.segs)
m <- function(val){
  factor(val, c("Two edges", "Four edges"), c(
    "Previous model:
one segment
per peak",
    "Proposed model:
several segments
per peak"))
}
both.lines <- rbind(
  data.table(model=m("Four edges"), fpop.lines),
  data.table(model=m("Two edges"), pdpa.lines))
pdpa.err <- PeakError::PeakErrorChrom(pdpa.peaks, my.labels)
common <- c(
  "chromStart", "chromEnd", "annotation", 
  "fp", "fn", "status")
both.err <- rbind(
  data.table(model=m("Four edges"), lab.err[peaks.int==P, ..common]),
  data.table(model=m("Two edges"), pdpa.err[, common]))
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(model ~ .)+
  penaltyLearning::geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd,
    fill=annotation),
    color="grey",
    size=0.25,
    data=my.labels)+
  penaltyLearning::geom_tallrect(aes(
    xmin=chromStart, xmax=chromEnd,
    linetype=status),
    fill=NA,
    data=both.err)+
  scale_fill_manual("label", values=ann.colors)+
  scale_linetype_manual(
    "error type",
    limits=c("correct", 
             "false negative",
             "false positive"),
    values=c(correct=0,
             "false negative"=3,
             "false positive"=1))+
  geom_step(aes(
    chromStart, count, color=type),
    data=data.table(type="data", Mono27ac$coverage))+
  geom_line(aes(
    pos, mean, color=type),
    size=0.7,
    alpha=0.7,
    data=data.table(type="model", both.lines))+
  geom_point(aes(
    (chromStart+chromEnd)/2, peak.y, color=type),
    shape=1,
    data=data.table(type="model", model=m("Two edges"), pdpa.peaks))+
  geom_point(aes(
    extremeMid, peak.y, color=type),
    shape=1,
    data=data.table(
      type="model",
      model=m("Four edges"), fit.state[annotation=="peakStart"& peaks.int==P]))+
  xlab("Position on chromosome 11")+
  scale_y_continuous(
    "Aligned DNA sequences
(H3K27ac histone modification)"
  )+
  coord_cartesian(xlim=c(15e4, 58e4), ylim=c(-6, 44))+
  guides(color="none")+
  scale_color_manual(values=type.colors)+
  geom_text(aes(
    x,y,label=label, color=type),
    hjust=0,
    size=2.9,
    data=data.table(model=m("Two edges"), rbind(
      data.table(x=325000, y=20, label="Noisy DNA sequence count data", type="data"),
      data.table(x=51e4, y=20, label="Model mean", type="model"),
      data.table(x=335000, y=peak.y, label="Model peaks", type="model"))))
png("figure-Mono27ac-summits-vs-peaks.png", 7, 2.5, res=300, units="in")
print(gg)
dev.off()
