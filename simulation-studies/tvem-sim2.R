rm(list=ls())
library("data.table")
library("tvem")
source("../time_varying-helper_funs.R")
source("../time_varying-lagged_effect-wrapper.R")

# simulation settings
simsettings <- expand.grid("uly"=c(TRUE),
                           "tes"=c(TRUE),
                           "n_trt"=c(10,20)
)
simsettings

# initialize for parallel MC jobs
args <- 1
# arguments from command line 
if(Sys.getenv("RSTUDIO") != "1" && !grepl("apple",Sys.getenv("R_PLATFORM"))) {
  args <- commandArgs(trailingOnly=TRUE) # for CMD BATCH '--args 1'
}
(seed <- as.integer(args[1]))

(simsettings[seed,])

# generate data ###############################################################
N <- 5000L # sample size
T_ <- as.integer(simsettings[seed,"n_trt"]) # treatment time points
REP <- 500L # number of sims

(uly <- simsettings[seed,"uly"]*2) # Ut --> Lt,Yt coefficients
# Xt --> Ys effects
delta.dt <- data.table(expand.grid("d"=1:2,t=1:T_,"deltaX"=NA_real_,
                                   "deltaXY"=NA_real_))
setkey(delta.dt)
if (simsettings[seed,"tes"]==FALSE) {
  # stable time-invariant effect
  delta.lag1 <- 0.5
  delta.lag2 <- 0.25
} else {
  # time-varying effects
  delta.lag1 <- sapply(1:T_, function(x) exp((T_/2)-x)/(1+exp((T_/2)-x)) )
  delta.lag2 <- sapply(1:(T_-1), function(x) 2*(x/T_-0.6)^2 )
}
delta.dt[d==1,deltaXY := delta.lag1]
delta.dt[d==2 & (d+t)<=(T_+1L),deltaX := delta.lag2]
delta.dt[d==1,deltaX := 0.0]
delta.dt[d==2 & (d+t)<=(T_+1L),deltaXY := 0.0]
setkey(delta.dt)

OneData <- function(n=10) {
  # unmeasured confounder
  U <- rnorm(n)
  # measured time-invariant confounder
  C <- rnorm(n)
  # initialize
  data_long.t <- data.table(expand.grid(
    id=1:n,time=1:(T_+1L),L=NA_real_,Y=NA_real_,X=NA_integer_))
  data_long.t[, C := rep(C,T_+1L)]
  setkey(data_long.t)
  for (tt in 1:(T_+1L)) {
    Ut <- (U^2) + rnorm(n)
    Lt <- uly*Ut + rnorm(n, sd=tt/T_)
    Yt <- uly*Ut + rnorm(n, sd=tt/T_)
    Yt <- Yt + log(abs(C)) + (C*Ut)
    if (tt>1) {
      Xt_minus1 <- data_long.t[time==tt-1, X]
      Yt_minus1 <- data_long.t[time==tt-1, Y]
      # lag 1 treatment effects
      Lt <- Lt + 4*(Xt_minus1-0.5)
      Yt <- Yt + delta.dt[d==1L & t==tt-1, deltaXY]*Yt_minus1*Xt_minus1
    }
    if (tt>2) {
      Xt_minus2 <- data_long.t[time==tt-2, X]
      # lag 2 treatment effects
      Yt <- Yt + delta.dt[d==2L & t==tt-2, deltaX]*Xt_minus2
    }
    Yt <- Yt-mean(Yt) # mean-center Yt at each time t
    Xt.ast <- min(0.35,1/var(Lt))*Lt + min(0.35,1/var(Yt))*Yt + 0.35*C
    Xt <- rbinom(n,1,exp(Xt.ast)/(1+exp(Xt.ast)))
    rm(Xt.ast)
    data_long.t[time==tt, c("L","Y","X"):=list(Lt,Yt,Xt)]
    
    rm(Lt,Yt,Xt)
  }
  setcolorder(data_long.t, c("id","time","C","L","Y","X"))
  setkey(data_long.t)
  return(data_long.t)
}

meths <- c("true_effect","tvem_only","time_int_lin","all_separate","gest_tvem")
simres <- NULL
ptm <- proc.time()[3]
for (i in 1:REP) {
  Data <- OneData(n=N)
  setkey(Data)
  
  res <- NULL
  for (j in 1:length(meths)) {
    # apply each method in turn to the same dataset
    meth.j <- meths[j]
    # same format for saving estimates using each method
    res.j <- data.table(delta.dt)
    res.j[, c("deltaX","deltaXY") := NA_real_]
    
    if (meth.j=="true_effect") {
      res.j <- data.table(delta.dt)
    } else if (meth.j=="tvem_only") {
      # shift by lag 2 then run tvem
      Data_lag2 <- data.table(Data)
      Data_lag2[, c("Xtplus1","Ltplus1","Ytplus1","Ytplus2") := list(
        c(.SD[-1,X],NA),
        c(.SD[-1,L],NA),
        c(.SD[-1,Y],NA),
        c(.SD[-(1:2),Y],NA,NA)), by=id]
      setkey(Data_lag2)
      Data_lag2[, XtYtplus1 := Xtplus1*Ytplus1]
      Data_lag2[, XtYt := X*Y]
      res.lag2 <- tvem(data=Data_lag2,
                       formula=Ytplus2~Xtplus1+X+Ltplus1+L+Ytplus1+Y+
                         XtYtplus1+XtYt,
                       invar_effect=~C,
                       num_knots=T_/2,penalize=TRUE,grid=1:T_,
                       id=id,time=time,use_naive_se=TRUE)
      
      # lag 1
      res.j[d==1, deltaX := 
              res.lag2$grid_fitted_coefficients$Xtplus1[,"estimate"]]
      res.j[d==1, deltaXY := 
              res.lag2$grid_fitted_coefficients$XtYtplus1[,"estimate"]]
      # lag 2
      res.j[d==2, deltaX := 
              res.lag2$grid_fitted_coefficients$X[,"estimate"]]
      res.j[d==2, deltaXY := 
              res.lag2$grid_fitted_coefficients$XtYt[,"estimate"]]
      
    } else {
      # set up data for g-estimation
      Data.ys <- data.table(Data)
      Data.ys[,time:=time-1] # g-estimation function indexes first wave as t=0
      setnames(Data.ys,"time","t")
      
      if (meth.j=="time_int_lin") {
        # linear interaction term with time
        res.all_lags <- OneEstLagged(
          Data=Data.ys, T_,
          Lnames=c("C","L"), Treat_name="X", Outcome_name="Y", 
          Znames="Y", use.DR=TRUE, time.coef.fun="interaction")
      
        # populate effects for each lag-treatment time combination
        for (dt_i in 1:nrow(res.j)) {
          # lag
          dt_i.coefs <- res.all_lags$est[[res.j[dt_i,d]]]
          main.eff <- grep("Y:",names(dt_i.coefs),
                           value=TRUE,fixed=TRUE,invert=TRUE)
          ints.eff <- grep("Y:",names(dt_i.coefs),
                           value=TRUE,fixed=TRUE,invert=FALSE)
          
          # treatment time: note that first time from g-est function is t=0
          res.j[dt_i,deltaX := crossprod(c(1,res.j[dt_i,t-1L]), 
                                         dt_i.coefs[main.eff])[1,1]]
          res.j[dt_i,deltaXY := crossprod(c(1,res.j[dt_i,t-1L]),
                                          dt_i.coefs[ints.eff])[1,1]]
          rm(dt_i.coefs)
        }
      } else if (meth.j=="all_separate") {
        # each outcome is handled as a separate end-of-study outcome
        res.j <- lapply(1:T_, function(ys) {
          res.ys <- OneEst(Data.ys,ys,Lnames=c("C","L"),
                           Treat_name="X",Outcome_name="Y",
                           Znames="Y", use.DR=TRUE)
          delta.est <- tail(lapply(res.ys$est, function(dt_i.coefs) {
            main.eff <- grep("Y_.",names(dt_i.coefs),
                             value=TRUE,fixed=TRUE,invert=TRUE)
            ints.eff <- grep("Y_.",names(dt_i.coefs),
                             value=TRUE,fixed=TRUE,invert=FALSE)
            c(dt_i.coefs[main.eff],tail(dt_i.coefs[ints.eff],n=1))
          }),n=2) # only lag 1 and 2 effects
          if (length(delta.est)==1L) {
            delta.est[[2]] <- delta.est[[1]]
            delta.est[[1]] <- c(NA,NA)
          }
          return(cbind("d"=2:1,"t"=c(ys-1L,ys),"s"=ys+1L,
                       do.call(rbind,delta.est)))
        })
        res.j <- data.table(do.call(rbind,res.j))
        res.j[, s := NULL]
        res.j <- res.j[t>0]
        setnames(res.j,names(delta.dt))
        
      } else if (meth.j=="gest_tvem") {
        # g-estimation + TVEM
        res.all_lags <- OneEstLagged(
          Data=Data.ys,T_,
          Lnames=c("C","L"), Treat_name="X", Outcome_name="Y", 
          Znames="Y", use.DR=TRUE, time.coef.fun="tvem")
        
        res.j <- data.table(do.call(rbind,res.all_lags$est[1:2]))
        setnames(res.j,names(delta.dt))
      }
    }
    setkey(res.j)
    # long format with main or interaction effect as a column
    res.j <- melt(
      res.j,id.vars=colnames(res.j)[1:2],
      measure.vars=c("deltaX","deltaXY"),value.name=meth.j,
      variable.name="mainORint",variable.factor=FALSE)
    setkey(res.j)
    res[[j]] <- res.j
  }
  
  # merge all the estimates into a single data.table
  res.dt <- Reduce(function(...)
    merge(..., by=c("d","t","mainORint"),all=TRUE), 
    res[!unlist(lapply(res,is.null))])
  setkey(res.dt)
  res.dt <- res.dt[!is.na(true_effect)]
  setkey(res.dt)
  
  # long format with method type as a column
  res.dt.long <- melt(
    res.dt,id.vars=colnames(res.dt)[1:4],
    measure.vars=meths[-1],value.name="est",
    variable.name="meth",variable.factor=FALSE)
  
  simres[[i]] <- res.dt.long
  
  cat(i,"|",round((proc.time()[3]-ptm)/60),"mins \n")
}

simres.dt <- rbindlist(simres)
simres.dt <- cbind(data.table(simsettings[seed,]),simres.dt)
setkey(simres.dt)
## mean absolute bias for each lag, across all treatment times
simres.dt[, bias := est-true_effect]
simres.dt[, mean(abs(bias)), by=eval(key(simres.dt)[c(1:4,6,8)])]
simres.dt[, min(abs(bias)), by=eval(key(simres.dt)[c(1:4,6,8)])]
simres.dt[, max(abs(bias)), by=eval(key(simres.dt)[c(1:4,6,8)])]

save(simres.dt, file=paste0("tvem-sim2-",seed,".Rdata"))

q()

# load results ################################################################
rm(list=ls())
library("data.table")
library("xtable")
files_to_load <- list.files()[grep(".Rdata",list.files())]
files_to_load <- grep("tvem-sim2",files_to_load,fixed=TRUE,value=TRUE)
simres <- NULL
seed <- 1L
for (myfile in files_to_load) {
  load(file=myfile)
  simres[[seed]] <- simres.dt
  rm(simres.dt)
  seed <- seed+1L
}
simres.dt <- rbindlist(simres)
setkey(simres.dt)

# number of simulations
simres.dt[, .N, by=eval(key(simres.dt)[1:8])][,unique(N)]

# results for each lag averaged over all time points
res.table <- NULL
res.table[["bias.abs"]] <- simres.dt[,mean(abs(bias),na.rm=TRUE),
                                     by=eval(key(simres.dt)[c(1:8)])][
                                       ,mean(V1), by=eval(key(simres.dt)[c(1:4,6,8)])]
setnames(res.table[["bias.abs"]],"V1","bias.abs")
res.table[["bias"]] <- simres.dt[,mean(bias,na.rm=TRUE),
                                 by=eval(key(simres.dt)[c(1:8)])][
                                   ,mean(V1), by=eval(key(simres.dt)[c(1:4,6,8)])]
setnames(res.table[["bias"]],"V1","bias")
res.table[["ese"]] <- simres.dt[,sd(est,na.rm=TRUE),
                                by=eval(key(simres.dt)[c(1:8)])][
                                  ,mean(V1), by=eval(key(simres.dt)[c(1:4,6,8)])]
setnames(res.table[["ese"]],"V1","ese")
res.table[["mse"]] <- simres.dt[,sqrt(mean(bias*bias,na.rm=TRUE)),
                                by=eval(key(simres.dt)[c(1:8)])][
                                  ,mean(V1), by=eval(key(simres.dt)[c(1:4,6,8)])]
setnames(res.table[["mse"]],"V1","rmse")
(res.table <- Reduce(function(...) merge(..., all = TRUE), res.table))

# rename methods
res.table[meth=="tvem_only", meth := "1.tvem_only"]
res.table[meth=="time_int_lin", meth := "2.time_int_lin"]
res.table[meth=="all_separate", meth := "3.all_separate"]
res.table[meth=="gest_tvem", meth := "4.gest_tvem"]
setkey(res.table)

print(xtable(res.table,include.rownames=FALSE))
