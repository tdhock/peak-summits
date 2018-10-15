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

## 2. Estimate spikes and calcium concentrations corresponding to eqn. (2) and (3)
fit_2 <- estimate_spikes(dat = data$calcium_d, gam = gam, lambda = lambda_p2, 
                         constraint = T, estimate_calcium = T)


times <- subset_is * (1 / fps)
dt <- data.table(seconds=times, calcium=data$calcium_d, AR1=fit_2$estimated_calcium)
fwrite(dt[, list(calcium)], "../../Cpattern1D/data.txt", col.names=FALSE)
penalty <- 5
res.dt <- fread(paste("cd ../../Cpattern1D && bin/Debug/multimodal", penalty))
seconds.between.data <- 1/fps
sec.w <- seconds.between.data/3
res.dt[, start.i := start +1L]
res.dt[, end.i := c(nrow(dt), start.i[-.N]-1L)]
res.dt[, start.seconds := dt$seconds[start.i]-sec.w ]
res.dt[, end.seconds := dt$seconds[end.i]+sec.w ]
res.dt[, Multimodal := ifelse(mean<0.1, 0, mean) ]
over.dt <- res.dt[dt, on=list(start.seconds < seconds, end.seconds > seconds)]
tall.dt <- melt(
  over.dt,
  measure.vars=c("Multimodal", "AR1"),
  variable.name="model")
tall.dt[, change.after := c(diff(value), NA)]
tall.dt[, seconds.after := c(start.seconds[-1], NA)]
tall.dt[, spike.i := cumsum(change.after < 0)]
tall.dt[, thresh := ifelse(model=="Multimodal", 0.5, 0)]
tall.dt[, thresh := ifelse(model=="Multimodal", 0, 0)]
m <- function(model){
  factor(
    model,
    c("AR1", "Multimodal"),
    c("Previous model:
AR1 changepoint
Jewell et al 2017", "Proposed model:
changepoint with
graph constraints"))
}
tall.dt[, model.fac := m(model)]
spike.dt <- tall.dt[0 < change.after & thresh < value, list(
  start.seconds=min(start.seconds),
  end.seconds=max(start.seconds)
), by=list(spike.i, model.fac)]
spike.dt[, mid.seconds := (start.seconds+end.seconds)/2]
spike.dt[169.75 < start.seconds & start.seconds < 169.85]

gg <- ggplot()+
  theme_bw()+
  ##theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.fac)+
  geom_point(aes(
    seconds, calcium),
    color="grey50",
    data=dt)+
  geom_line(aes(
    start.seconds, value),
    data=tall.dt,
    color="green")+
  geom_vline(aes(
    xintercept=mid.seconds),
    color="green",
    linetype="dotted",
    data=spike.dt)

if(FALSE){
  
  gg+
    geom_vline(xintercept=c(168, 172))

  gg+coord_cartesian(xlim=c(0, 10))

  gg+coord_cartesian(xlim=c(168, 172))

  gg+coord_cartesian(xlim=c(160, 180))

  gg+coord_cartesian(xlim=c(169.5, 170))
  
}

lab <- function(xmin, xmax){
  data.table(xmin, xmax, label="oneSpike", annotation="1change", problem=1)
}
label.dt <- rbind(
  lab(166.5, 166.75),
  lab(169.7, 169.9))
ann.colors <- c(
  noPeaks="#f6f4bf",
  peakStart="#ffafaf",
  oneSpike="#ffafaf",
  peakEnd="#ff4c4c",
  peaks="#a445ee")
xmin <- 166
xmax <- 171
spike.dt[, problem := 1]
models.dt <- spike.dt[, list(
  models=.N
  ), by=list(model.fac, problem)]
err.list <- penaltyLearning::labelError(
  models.dt, label.dt, spike.dt,
  change.var="mid.seconds",
  problem.var="problem",
  label.vars=c("xmin", "xmax"),
  model.vars="model.fac")
type.colors <- c(
  data="grey50",
  model="blue")
show.spikes <- spike.dt[xmin < mid.seconds & mid.seconds < xmax]
show.data <- dt[xmin < seconds & seconds < xmax]
spike.y <- -1.5
gg.out <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(model.fac ~ .)+
  penaltyLearning::geom_tallrect(aes(
    xmin=xmin, xmax=xmax, fill=label),
    color=NA,
    data=label.dt)+
  penaltyLearning::geom_tallrect(aes(
    xmin=xmin, xmax=xmax, linetype=status),
    fill=NA,
    color="black",
    data=err.list$label.errors)+
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
    start.seconds, value, color=type),
    data=data.table(
      type="model",
      tall.dt[xmin < start.seconds & start.seconds < xmax]),
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
    "Fluorescence intensity
(Measure of neural activity)",
    breaks=seq(0, 10, by=2),
    limits=c(-2, 10)
  )+
  scale_x_continuous("Time (seconds)")+
  guides(color="none")+
  scale_color_manual(values=type.colors)+
  geom_text(aes(
    x,y,label=label, color=type, hjust=hjust),
    size=3,
    data=data.table(model.fac=m("AR1"), rbind(
      data.table(
        hjust=0, x=167, y=2, label="Noisy activity data", type="data"),
      data.table(
        hjust=1, x=169.5, y=5, label="Mean model", type="model"),
      data.table(
        hjust=1, x=169.4, y=spike.y, label="Predicted spikes", type="model"))))
##print(gg.out)
png("../../figure-AR1-multimodal.png", 7, 2.5, units="in", res=300)
print(gg.out)
dev.off()
##system("display ../../figure-AR1-multimodal.png")
