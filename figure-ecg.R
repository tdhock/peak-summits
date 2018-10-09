works_with_R(
  "3.5.1",
  ggplot2="3.0.0",
  data.table="1.11.6")

csv.vec <- Sys.glob("ecg-data/output/*.csv")
ecg <- fread(csv.vec[6], fill=TRUE)
ecg.mat <- t(as.matrix(ecg))
ecg.dt <- data.table(Y=ecg.mat[,1], t=ecg.mat[,2])
model.dt <- data.table(ecg.mat[,-(1:2)])[!is.na(V1)]
setnames(model.dt, c(
  "Y.P", "t.P",
  "Y.Q", "t.Q",
  "Y.R", "t.R",
  "Y.S", "t.S",
  "Y.T", "t.T"))
molt <- melt(model.dt, measure=patterns("^Y.", "^t."))
setnames(molt, c("var.int", "Y", "t"))
molt[, letter := c("P", "Q", "R", "S", "T")[var.int] ]
## row 1  is all of the Y values for the ECGII signal
## row 2  is all of the t values for the ECGII signal
## row 3  is all of the Y values for the P peaks
## row 4  is all of the t values for the P peaks
## row 5  is all of the Y values for the Q peaks
## row 6  is all of the t values for the Q peaks
## row 7  is all of the Y values for the R peaks
## row 8  is all of the t values for the R peaks
## row 9  is all of the Y values for the S peaks
## row 10 is all of the t values for the S peaks
## row 11 is all of the Y values for the T peaks
## row 12 is all of the t values for the T peaks

some <- function(dt){
  dt[t %between% c(52e3, 53e3)]
}
some.ecg <- some(ecg.dt)
some.model <- some(molt)

gg <- ggplot()+
  theme_bw()+
  geom_line(aes(
    t, Y),
    color="grey50",
    data=some.ecg)+
  ## geom_point(aes(
  ##   t, Y),
  ##   color="grey50",
  ##   data=ecg.dt)+
  geom_text(aes(
    t, Y, label=letter),
    data=some.model[letter %in% c("Q", "R", "S")])
print(gg)

fwrite(some.ecg[, list(Y)], "Cpattern1D/data.txt", col.names=FALSE)


## graph << Edge(beta, 0, 1, "down",   0);
## graph << Edge(0,    1, 2, "up", 2000);
## graph << Edge(0,    2, 3, "down", 5000);
## graph << Edge(0,    3, 4, "up", 2000);
## graph << Edge(0,    4, 5, "up", 1000);
## graph << Edge(0,    5, 6, "up", 0);
## graph << Edge(0,    6, 7, "down", 0);
## graph << Edge(0,    7, 8, "down", 0);
## graph << Edge(0,    8, 0, "up", 0);
penalty <- 80000000
res.dt <- fread(paste("cd Cpattern1D && bin/Debug/9states2", penalty))
res.dt[, start.i := start +some.ecg$t[1]]
res.dt[, end.i := c(some.ecg[.N, t], start.i[-.N]-1L)]
res.dt[, cum.0 := cumsum(state==0)]
res.dt[, letter := c("beforeQ", "Q", "R", "S", "S1", "S2", "peak", "afterPeak", "foo")[state+1] ]
mean.dt <- res.dt[order(start.i), data.table(
  t=as.numeric(rbind(start.i-0.5, end.i+0.5)),
  Y=as.numeric(rbind(mean, mean)))]
gg+
  geom_line(aes(
    t, Y),
    data=mean.dt,
    color="green")+
  geom_point(aes(
  (start.i+end.i)/2, mean),
  color="black",
  fill="green",
  shape=21,
  data=res.dt[letter=="R"])+
  coord_cartesian(xlim=c(52000, 52900), expand=FALSE)

m <- function(x){
  factor(x, c("previous", "changepoint"), c(
    "Previous model", "Graph-constrained\nchangepoint model"))
}
model.dt <- rbind(res.dt[, data.table(
  model=m("changepoint"),
  Y=mean,
  t=ifelse(letter=="Q", end.i, (start.i+end.i)/2),
  letter)], some.model[, data.table(
    model=m("previous"),
    Y, t, letter)])[letter %in% c("Q", "R", "S")]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(model ~ .)+
  geom_line(aes(
    t, Y),
    color="grey50",
    data=some.ecg)+
  geom_line(aes(
    t, Y),
    data=data.table(model=m("changepoint"), mean.dt),
    color="blue")+
  geom_label(aes(
    t, Y, label=letter),
    color="blue",
    size=3,
    label.padding=grid::unit(0.1, "lines"),
    alpha=0.6,
    data=model.dt)+
  coord_cartesian(xlim=c(52000, 52900), expand=FALSE)
png("figure-ecg.png", 10, 4, res=300, units="in")
print(gg)
dev.off()
