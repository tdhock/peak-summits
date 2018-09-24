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

## plot parameters 
wd <- 0.1 # rectangle box width 

## 1. Read and process spike finder data
data_dir <- paste0(spikefinder_base_dir, "spikefinder.train/", data_set, ".train.")
calcium_dat <- readr::read_csv(paste0(data_dir, "calcium.csv")) %>% rename_all(col_renamer)
spike_dat <- readr::read_csv(paste0(data_dir, "spikes.csv")) %>% rename_all(col_renamer)
data <- subset_data(calcium_dat, spike_dat, subset_is)


## 2. Estimate spikes and calcium concentrations corresponding to eqn. (2) and (3)
fit_2 <- estimate_spikes(dat = data$calcium_d, gam = gam, lambda = lambda_p2, 
                         constraint = T, estimate_calcium = T)

fit_3 <- estimate_spikes(dat = data$calcium_d, gam = gam, lambda = lambda_p3, 
                         constraint = F, estimate_calcium = T)


## 3. Plot results from (2) and (3) around a region of interest and 
## 4. Save plot as pdf in specified location 
neg_spikes <- which(fit_3$estimated_calcium[fit_3$spikes] - gam * fit_3$estimated_calcium[fit_3$spikes - 1] < 0 )
neg_spikes_magnitude <- fit_3$estimated_calcium[fit_3$spikes[neg_spikes]] -gam * fit_3$estimated_calcium[fit_3$spikes[neg_spikes - 1]]
large_neg_spike_idx <- sort(neg_spikes_magnitude, index.return = T)$ix[1]

## plot around large negative spike
neg_spike_time <- fit_3$spikes[neg_spikes[large_neg_spike_idx]] * (1 / fps)
ts <- neg_spike_time + 2 * c(-1, 1)

## setup plot
times <- subset_is * (1 / fps)
true_spike <- which(data$spikes > 0) * (1 / fps)

##pdf(paste0(fig_save_dir, "fig7.pdf"), height = 8, width = 16)
##par(mfrow = c(1, 2), mar = c(5, 3, 4, 2) + 0.1)
plot_estimates(times, data$calcium_d, true_spike, fit_3$spikes, xlim = ts, xlab = "Time (s)")
lines(times, fit_3$estimated_calcium, col = "blue")
rect(neg_spike_time - wd, -1.25, neg_spike_time + wd, 10, border = "darkred", lwd = 2)

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
tall.dt[, model.fac := factor(model, c("AR1", "Multimodal"), c("Autoregressive order 1 model\nJewell et al 2017", "Multimodal regression model\nProposal based on graph constraints"))]
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

gg+
  geom_vline(xintercept=c(168, 172))

gg+coord_cartesian(xlim=c(0, 10))

gg+coord_cartesian(xlim=c(168, 172))

gg+coord_cartesian(xlim=c(160, 180))

gg+coord_cartesian(xlim=c(169.5, 170))

lab <- function(xmin, xmax){
  data.table(xmin, xmax, label="oneSpike", annotation="1change", problem=1)
}
label.dt <- rbind(
  lab(166.5, 166.75),
  lab(169.6, 169.9))
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
gg.out <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ model.fac)+
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
    seconds, calcium),
    color="grey50",
    shape=1,
    data=dt[xmin < seconds & seconds < xmax])+
  geom_line(aes(
    start.seconds, value),
    data=tall.dt[xmin < start.seconds & start.seconds < xmax],
    size=0.5,
    color="green")+
  geom_vline(aes(
    xintercept=mid.seconds),
    color="green",
    size=0.2,
    linetype="dashed",
    data=spike.dt[xmin < mid.seconds & mid.seconds < xmax])+
  scale_y_continuous("Calcium concentration")+
  scale_x_continuous("Time (seconds)")
##print(gg.out)
png("../../figure-AR1-multimodal.png", 7, 3, units="in", res=200)
print(gg.out)
dev.off()
##system("display ../../figure-AR1-multimodal.png")

## plot_estimates(times, data$calcium_d, true_spike, fit_2$spikes, xlim = ts, xlab = "Time (s)")
## lines(times, fit_2$estimated_calcium, col = "blue")
## rect(neg_spike_time - wd, -1.25, neg_spike_time + wd, 10, border = "darkred", lwd = 2)
##dev.off()


