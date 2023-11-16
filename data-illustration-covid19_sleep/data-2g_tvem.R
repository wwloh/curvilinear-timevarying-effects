# load organized data for analysis
source("data_1prep.R")
library("tvem")
# load helper functions
if (grepl("apple",sessionInfo()[[1]]$platform)) {
  source("../time_varying-helper_funs.R")
  source("../time_varying-lagged_effect-wrapper.R")
} else {
  source("time_varying-helper_funs.R")
  source("time_varying-lagged_effect-wrapper.R")
}

Data.OBS <- data.table(Data)
setkey(Data.OBS)
rm(Data)

# parallel on desktop
library("parallel")
n.cores <- detectCores(logical = FALSE)
(n.cores <- round(n.cores/2)) # use half the number of available cores

OneDataEst <- function(Data) {

meths <- c("all_separate","gest_tvem")  
res <- NULL
for (j in 1:length(meths)) {
  # apply each method in turn to the same dataset
  meth.j <- meths[j]
  
  # set up data for g-estimation
  Data.ys <- data.table(Data)
  Data.ys[,t:=t-1] # g-estimation function indexes first wave as t=0
  n.lags <- 4L # lags of effects
  
  if (meth.j=="all_separate") {
    # each outcome is handled as a separate end-of-study outcome
    res.j <- lapply(1:T_, function(ys) {
      res.ys <- OneEst(Data.ys,ys,Lnames=c(Cnames,Lnames),
                       Treat_name=Xname,Outcome_name=Yname,use.DR=TRUE)
      delta.est <- tail(unlist(res.ys$est),n=n.lags)
      if (length(delta.est)<n.lags) {
        delta.est <- c(rep(NA,n.lags-length(delta.est)),delta.est)
      }
      return(cbind("d"=n.lags:1L,"t"=ys-((n.lags:1L)-1L),"s"=ys+1L,delta.est))
    })
    res.j <- data.table(do.call(rbind,res.j))
    res.j[, s := NULL]
    res.j <- res.j[t>0]
  } else if (meth.j=="gest_tvem") {
    # g-estimation + TVEM
    res.all_lags <- OneEstLagged(
      Data=Data.ys,T_,Lnames=c(Cnames,Lnames),
      Treat_name=Xname,Outcome_name=Yname, use.DR=TRUE,
      time.coef.fun="tvem")
    
    res.j <- data.table(do.call(rbind,res.all_lags$est[1:n.lags]))
  }
  setnames(res.j,c("d","t",meth.j))
  setkey(res.j)
  res[[j]] <- res.j
}

# merge all the estimates into a single data.table
res.dt <- Reduce(function(...)
  merge(..., by=c("d","t"),all=TRUE), 
  res[!unlist(lapply(res,is.null))])
setkey(res.dt)

# long format with method type as a column
res.dt.long <- melt(
  res.dt,id.vars=colnames(res.dt)[1:2],
  measure.vars=,value.name="est",
  variable.name="meth",variable.factor=FALSE)

return(res.dt.long)

}

n.boots <- 2000
batch.seeds <- matrix(1:n.boots, nrow=n.cores)
for (j in 1:ncol(batch.seeds)) {
  
  # bootstrap in parallel
  res.parallel <- mclapply(batch.seeds[,j], function(seed) {
    if (seed==1) {
      # observed data
      Data <- Data.OBS
    } else {
      # resample with replacement for bootstrap
      ids <- unique(Data.OBS$id)
      boot.i <- sort(sample(x=ids,size=length(ids),replace=TRUE))
      Data <- do.call(rbind,lapply(1:length(boot.i), function(ib) {
        D.ib <- Data.OBS[id==boot.i[ib]]
        D.ib[,id := ib]
        return(D.ib)
      }))
    }
    setkey(Data)
    OneDataEst(Data)
  }, mc.cores=n.cores)
  
  # save results intermittently
  for (seed.i in 1:length(res.parallel)) {
    seed <- batch.seeds[seed.i,j]
    res.dt.long <- res.parallel[[seed.i]]
    myfilename <- "data-tvem/bootstrap_results"
    save(res.dt.long,file=paste0(myfilename,"-seed_",seed,".Rdata"))
    rm(seed,res.dt.long)
  }
  rm(res.parallel)
}

q()
