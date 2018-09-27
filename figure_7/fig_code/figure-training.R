works_with_R(
  "3.5.1",
  penaltyLearning="2018.9.4",
  data.table="1.11.6",
  ggplot2="3.0.0",
  tidyverse="1.2.1",
  latex2exp="0.4.0",
  "jewellsean/FastLZeroSpikeInference@e4fcf32333bac8d01122b5352032b1aa28c6f714")

source("utils.R")

## data parameters
fps <- 100 # frames per second 
data_set <- 7 # spike finder dataset number 
ind <- 13 # cell number 
subset_is <- 1:20000 # subset of the data used to estimate the solution 
gam <- 0.98 # decay parameter 

## tuning parameters 
lambda_p2 <- 3.6
lambda_p3 <- 3.76

## file i/o paramters 
spikefinder_base_dir <- "../spike_finder_data/"
fig_save_dir <- "../figures/"

## 1. Read and process spike finder data
data_dir <- paste0(spikefinder_base_dir, "spikefinder.train/", data_set, ".train.")
calcium_dat <- readr::read_csv(paste0(data_dir, "calcium.csv")) %>% rename_all(col_renamer)
spike_dat <- readr::read_csv(paste0(data_dir, "spikes.csv")) %>% rename_all(col_renamer)
data <- subset_data(calcium_dat, spike_dat, subset_is)

times <- subset_is * (1 / fps)
dt <- data.table(seconds=times, calcium=data$calcium_d)
fwrite(dt[, list(calcium)], "../../Cpattern1D/data.txt", col.names=FALSE)

pen.vec <- c(
  "Too few spikes"="500",
  "Zero label errors"="0.5",
  "Too many spikes"="0.05")
m <- function(val){
  factor(val, names(pen.vec), sprintf("Penalty=%s\n%s", pen.vec[names(pen.vec)], names(pen.vec)))
}
mean.dt.list <- list()
for(pen.name in names(pen.vec)){
  penalty <- pen.vec[[pen.name]]
  res.dt <- fread(paste("cd ../../Cpattern1D && bin/Debug/multimodal", penalty))
  res.dt[, orig.mean := mean]
  seconds.between.data <- 1/fps
  sec.w <- seconds.between.data/3
  res.dt[, start.i := start +1L]
  res.dt[, end.i := c(nrow(dt), start.i[-.N]-1L)]
  res.dt[, start.seconds := dt$seconds[start.i]-sec.w ]
  res.dt[, end.seconds := dt$seconds[end.i]+sec.w ]
  res.dt[, trunc.mean := ifelse(orig.mean < 0.5, 0, orig.mean) ]
  res.dt[, zero.i := cumsum(trunc.mean != 0)]
  zero.means <- res.dt[trunc.mean==0, list(
    seg.mean=mean(orig.mean)
  ), by=list(zero.i)]
  res.dt[zero.means, mean := seg.mean, on=list(zero.i)]
  over.dt <- res.dt[dt, on=list(start.seconds < seconds, end.seconds > seconds)]
  over.dt[, change.after := c(diff(mean), NA)]
  over.dt[, seconds.after := c(start.seconds[-1], NA)]
  over.dt[, spike.i := cumsum(change.after < 0)]
  mean.dt.list[[pen.name]] <- data.table(pen.fac=m(pen.name), over.dt)
}
mean.dt <- do.call(rbind, mean.dt.list)
spike.dt <- mean.dt[0 < change.after & 0 < mean, list(
  start.seconds=min(start.seconds),
  end.seconds=max(start.seconds)
), by=list(spike.i, pen.fac)]
spike.dt[, mid.seconds := (start.seconds+end.seconds)/2]

ann.map <- c(
  oneSpike="1change",
  noSpikes="0changes")
lab <- function(xmin, xmax, label){
  data.table(xmin, xmax, label, annotation=ann.map[label], problem=1)
}
label.dt <- rbind(
  ##lab(125, 127, "noSpikes")
  lab(147, 148, "oneSpike"),
  lab(145, 146, "oneSpike"),
  lab(130, 132, "oneSpike"))
ann.colors <- c(
  noPeaks="#f6f4bf",
  noSpikes="#f6f4bf",
  peakStart="#ffafaf",
  oneSpike="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
xmin <- 120
xmax <- 150
spike.dt[, problem := 1]
models.dt <- spike.dt[, list(
  models=.N
  ), by=list(pen.fac, problem)]
err.list <- penaltyLearning::labelError(
  models.dt, label.dt, spike.dt,
  change.var="mid.seconds",
  problem.var="problem",
  label.vars=c("xmin", "xmax"),
  model.vars="pen.fac")
type.colors <- c(
  data="grey50",
  model="blue")
show.spikes <- spike.dt[xmin < mid.seconds & mid.seconds < xmax]
show.data <- dt[xmin < seconds & seconds < xmax]
spike.y <- -1.5
gg.out <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(pen.fac ~ .)+
  penaltyLearning::geom_tallrect(aes(
    xmin=xmin, xmax=xmax, fill=label),
    color="grey",
    size=0.5,
    data=label.dt)+
  penaltyLearning::geom_tallrect(aes(
    xmin=xmin, xmax=xmax, linetype=status),
    fill=NA,
    color="black",
    data=err.list$label.errors)+
  scale_fill_manual(values=ann.colors)+
  scale_linetype_manual(
    "error type",
    limits=c("correct", 
             "false negative",
             "false positive"),
    values=c(correct=0,
             "false negative"=3,
             "false positive"=1))+
  geom_point(aes(
    seconds, calcium, color=type),
    shape=1,
    data=data.table(type="data", show.data))+
  ## geom_line(aes(
  ##   seconds, calcium),
  ##   color="grey50",
  ##   data=show.data)+
  geom_line(aes(
    start.seconds, mean, color=type),
    data=data.table(
      type="model",
      mean.dt[xmin < start.seconds & start.seconds < xmax]),
    size=0.5)+
  geom_point(aes(
    mid.seconds, spike.y, color=type),
    shape=1,
    data=data.table(type="model", show.spikes))+
  ## geom_vline(aes(
  ##   xintercept=mid.seconds),
  ##   color=model.color,
  ##   size=0.5,
  ##   linetype="dashed",
  ##   data=show.spikes)+
  scale_y_continuous(
    "Calcium concentration (measure of neural activity)",
    breaks=seq(0, 10, by=2),
    limits=c(-2, 10)
  )+
  scale_x_continuous("Time (seconds)")+
  guides(color="none")+
  scale_color_manual(values=type.colors)+
  geom_text(aes(
    x,y,label=label, color=type, hjust=hjust),
    data=data.table(pen.fac=m("Too few spikes"), rbind(
      data.table(
        hjust=0, x=132.5, y=7, label="Noisy activity data", type="data"),
      data.table(
        hjust=0.5, x=141, y=8, label="Mean model", type="model"),
      data.table(
        hjust=1, x=127, y=spike.y, label="Predicted spikes", type="model"))))
##print(gg.out)
png("../../figure-neuro-training.png", 10, 4, units="in", res=300)
print(gg.out)
dev.off()
##system("display ../../figure-AR1-multimodal.png")
