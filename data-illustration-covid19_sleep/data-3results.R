rm(list=ls())
library("data.table")

subfolder <- "data-tvem/"
myfiles <- list.files(subfolder)
myfiles <- myfiles[grep(pattern=".Rdata",myfiles)]
simres <- NULL
for (ll in myfiles) {
  load(paste0(subfolder,ll))
  simres <- c(simres,list(res.dt.long))
  rm(res.dt.long)
}
obs.est <- simres[[1]] # observed estimate
setkey(obs.est)
simres <- simres[-1]
# remove bootstrap samples with extreme estimates
simres <- simres[unlist(lapply(simres, function(x) x[,max(abs(est))]))<10]
length(simres)
sim_res <- rbindlist(simres)
rm(simres)

res.boot <- sim_res[, as.list(
  c("se"=sd(est), quantile(est, probs=c(.025,.975)))), by=c("d","t","meth")]
setkey(res.boot)
res <- merge(obs.est,res.boot, by=c("d","t","meth"))
setkey(res)

pdf("res-plot.pdf", width=8,height=6)
par(mfrow=c(2,2))
for (dd in 1:max(res[,d])) {
  gest.allsep <- res[d==dd & meth=="all_separate"]
  gest.allsep.cols <- rep("grey50",nrow(gest.allsep))
  gest.allsep.cols[unlist(gest.allsep[,6,with=FALSE]>0) | 
                     unlist(gest.allsep[,7,with=FALSE]<0)] <- "red"
  plot(gest.allsep[,est], 
       ylim=range(res[,6:7,with=FALSE]), pch=20,col=gest.allsep.cols,
       ylab="effect estimate",xlab="study week",xaxt="n",yaxt="n",
       main=paste("lag", eval(dd)))
  axis(side=1,at=gest.allsep[,t],cex.axis=.8)
  axis(side=2,
       at=floor(min(res[,6:7,with=FALSE])):ceiling(max(res[,6:7,with=FALSE])),
       cex.axis=.75)
  lines(res[d==dd & meth=="gest_tvem",est],lwd=1.25)
  abline(h=0,lty=3)
  lines(unlist(res[d==dd & meth=="gest_tvem",6,with=FALSE]),lty=2)
  lines(unlist(res[d==dd & meth=="gest_tvem",7,with=FALSE]),lty=2)
  for (i in 1:nrow(gest.allsep)) {
    lines(rep(i,2),unlist(gest.allsep[i,6:7,with=FALSE]),lwd=.75,
          col=gest.allsep.cols[i])
  }
}
dev.off()
