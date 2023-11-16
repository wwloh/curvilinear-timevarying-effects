OneEstLagged <- function(
  Data, # long format
  T_,   # time point for last outcome
  Lnames, # covariate names
  Treat_name,
  Outcome_name,
  Znames=NULL, # effect modifiers in SNMM
  use.DR=TRUE,
  time.coef.fun="interaction"
  ) {
  ## names cannot be "P" or "Ry" or contain "_."
  
  # prepare dataset for fitting ###############################################
  D <- data.table(Data)
  setkey(D)
  
  # create rows with NAs for observations with missing waves
  if (D[,.N,by=id][,any(N<T_+1)]) {
    for (i in D[,.N,by=id][N<T_+1,id]) {
      # create new rows
      D <- rbind(D,data.table("id"=i,"t"=(0:T_)[!((0:T_) %in% D[id==i,t])]),
                 fill=TRUE)
      setkey(D)
    }
  }
  
  # set repeated covariates to value at baseline and to NA at other times
  Lbaseline <- apply(D[,lapply(.SD, function(x) 
    length(unique(x[!is.na(x)]))==1L), by="id"][
      ,!(colnames(D)%in%c("id","t")),with=FALSE],2,all)
  Lbaseline <- Lnames[Lnames %in% names(Lbaseline[Lbaseline])]
  D[t>0,(Lbaseline) := NA]
  rm(Lbaseline)
  
  # reshape to wide format
  timesep <- "_." # used to indicate time point for a variable
  D <- data.table::dcast(data = D, 
                         formula = id ~ t, 
                         value.var = names(D)[!(names(D) %in% c("id","t"))],
                         sep = timesep,
                         fill = NA)
  setkey(D)
  
  ## ignore variables (besides outcome) observed at or after time T
  var_T <- grep(paste0(timesep,T_),colnames(D),value=TRUE)
  D[, var_T[var_T!=paste0(Outcome_name,timesep,T_)]] <- NULL
  rm(var_T)
  
  ## remove columns with all NAs
  D <- D[, colSums(is.na(D))<nrow(D), with=FALSE]
  setkey(D)
  
  ## baseline-only covariates
  Lbaseline <- Lnames[
    colSums(sapply(sapply(Lnames,paste0,timesep),grepl,colnames(D)))==1]
  Lbaseline <- colnames(D)[sapply(
    sapply(Lbaseline,paste0,timesep),grep,colnames(D))]
  
  res.coef <- res.formulae <- res.fit <- tt.predictors.list <- NULL
  for (tt in (T_-1L):0) {
    # treatment name
    A_t <- paste0(Treat_name,timesep,tt)
    # PS name
    P_t <- paste0("P",timesep,tt)
    
    # variable history
    tt.predictors <- c(
      sapply(Lnames,paste,c(0,(tt-2L):tt),sep=timesep), # covariate history
      paste0(Outcome_name,timesep,c(0,(tt-2L):tt)),     # outcome history
      paste0(Treat_name,timesep,pmax(0,tt-(2:1)))       # treatment history
    )
    names(tt.predictors) <- NULL
    tt.predictors <- tt.predictors[(tt.predictors %in% colnames(D)) & 
                                     (tt.predictors != A_t)]
    tt.predictors <- sort(tt.predictors)
    
    tt.ps <- as.formula(paste0(A_t,"~",paste(tt.predictors,collapse="+")))
    
    # fit PS model
    ps.fit <- NULL
    if (use.DR==TRUE) {
      ps.fit <- glm(formula=tt.ps,data=D,family=binomial("logit"))
      ps.hat <- ps.fit$fitted.values
      
      if (nrow(D)>length(ps.hat)) {
        ps.hat.all <- rep(NA,nrow(D))
        ps.hat.all[as.integer(names(ps.hat))] <- ps.hat
        ps.hat <- ps.hat.all
        rm(ps.hat.all)
      }
      D[, paste0("P",timesep,tt) := ps.hat]
    }
    
    tt.predictors.list[[tt+1]] <- tt.predictors
    res.formulae[[tt+1]] <- list("PS"=tt.ps)
    res.fit[[tt+1]] <- list("PS"=ps.fit)
  }
    
  # consider each lag in turn
  res_lag <- NULL
  for (lag_d in 1:T_) {
    out.mm.list <- NULL
    for (ss in lag_d:T_) {
      # outcome name
      if (lag_d==1L) {
        Rt <- paste0(Outcome_name,timesep,ss)
      } else {
        Rt <- paste0("Ry",timesep,ss)
      }
      
      # treatment time
      tt <- ss-lag_d
      # treatment name
      A_t <- paste0(Treat_name,timesep,tt)
      # PS name
      P_t <- paste0("P",timesep,tt)
    
      ## predictors from SNMM
      if(is.null(Znames) || !("1" %in% Znames)) {
        Znames <- c("1",Znames)
      }
      
      tt.ZP <- unlist(lapply(Znames, function(z) {
        if (use.DR==TRUE) {
          z_ap <- c(A_t,P_t)
        } else {
          z_ap <- A_t
        }
        if (z=="1") {
          # constant unconditional effect
          return(z_ap)
        } else {
          # interaction with all previous occurrences of the covariate
          z_t <- grep(z,colnames(D),value=TRUE)
          if (length(z_t)>1) {
            ## keep only contemporaneous occurrence
            z_t <- z_t[as.integer(lapply(
              strsplit(z_t,split=paste0(z,timesep)),"[",2)) == tt]
          }
          z_t <- sapply(z_t, function(zz) paste0(zz,":",z_ap))
          colnames(z_t) <- NULL
          return(z_t)
        }
      }))
      
      # outcome model
      tt.predictors <- tt.predictors.list[[tt+1]]
      ## limit to covariates at same wave as treatment
      tt.predictors.samewave <- which(as.integer(unlist(lapply(
        strsplit(tt.predictors,split=timesep,fixed=TRUE),"[[",2)))==tt)
      tt.predictors.samewave <- tt.predictors[tt.predictors.samewave]
      
      tt.out.predictors <- unique(c(Lbaseline,tt.predictors.samewave,tt.ZP))
      tt.out <- as.formula(paste0(Rt,"~ -1 + ", ## no intercept
                                  paste(tt.out.predictors,collapse="+")))
      ## model matrix with interaction terms and NAs
      #### https://stackoverflow.com/questions/5616210/model-matrix-with-na-action-null
      out.mm <- model.matrix(tt.out, model.frame(
        formula = tt.out, data=D, na.action=na.pass))
      
      ## append outcome
      out.mm <- cbind(unlist(D[,..Rt]),out.mm)
      out.mm <- data.table(out.mm)
      setnames(out.mm,old=1,new=Outcome_name)
      
      # same time-varying confounders
      tt.out.names <- which(!(colnames(out.mm) %in% 
                                c("id",Outcome_name,tt.ZP,Lbaseline)))
      tt.out.names.new <- unlist(lapply(
        strsplit(colnames(out.mm)[tt.out.names],split=timesep),"[",1))
      setnames(out.mm,old=tt.out.names,
               new=paste0(tt.out.names.new,"_lag_",lag_d))
      
      ## relabel treatments in SNMM to common names for common lagged effects
      snmm.names <- lapply(strsplit(tt.ZP,":"),function(zp) {
        paste(lapply(sapply(zp,strsplit,timesep,fixed=TRUE),"[",1),collapse=":")
      })
      snmm.names <- paste0(snmm.names,"_lag_",lag_d)
      setnames(out.mm,old=which(colnames(out.mm) %in% tt.ZP),new=snmm.names)
      
      ## add ID and time point for treatment and lag
      out.mm[, id := D$id]
      out.mm[, trt_time := tt]
      setcolorder(out.mm,c(ncol(out.mm)-(1:0),1:(ncol(out.mm)-2)))
      setkey(out.mm)
      
      out.mm.list[[ss]] <- out.mm
      
      rm(tt.predictors,tt.predictors.samewave,tt.out.predictors,tt.out,tt.ZP,
         out.mm,snmm.names)
    }
    
    # construct response vector and design matrix
    D.stacked <- rbindlist(out.mm.list,use.names=TRUE,fill=FALSE)
    setkey(D.stacked)
    
    ## inspect a single individual at random
    if (FALSE) {
      some.id <- sample(D$id,1)
      D[id==some.id]
      D.stacked[id==some.id]
      Data[id==some.id]
      rm(some.id)
    }
    
    # fit outcome model
    out.preds <- colnames(D.stacked)[
      !(colnames(D.stacked) %in% c("id",Outcome_name,"trt_time"))]
    try.lm <- FALSE
    if (time.coef.fun=="tvem") {
      # quirks of tvem:
      #### no periods or colons in column names
      #### times must be positive (non-zero)
      out.preds.timevary <- out.preds[!(out.preds %in% Lbaseline)]
      out.preds.tvem <- sapply(out.preds.timevary,gsub,pattern=":",
                               replacement="_int_")
      Lbaseline.tvem <- sapply(Lbaseline,gsub,pattern=timesep,replacement="")
      
      # formula for time-varying covariates
      tvc.formula <- as.formula(paste0(Outcome_name, "~", paste(
        out.preds.tvem, collapse="+")))
      # formula for time-fixed covariates
      tfc.formula <- as.formula(paste0("~", paste(
        Lbaseline.tvem, collapse="+")))
      
      D.stacked.tvem <- data.table(D.stacked)
      D.stacked.tvem[, trt_time := trt_time+1]
      setnames(D.stacked.tvem,old=Lbaseline,new=Lbaseline.tvem)
      setnames(D.stacked.tvem,old=out.preds.timevary,new=out.preds.tvem)
      setkey(D.stacked.tvem)
      tt.values <- unique(D.stacked.tvem$trt_time)
      
      ptm <- tryCatch(system.time(
        out.fit <- tvem(data=D.stacked.tvem,
                        formula=tvc.formula,
                        invar_effect=tfc.formula,
                        num_knots=pmin(floor(length(tt.values)/2),10),
                        penalize=TRUE,grid=tt.values,
                        id=id,time=trt_time,use_naive_se=TRUE)
      )[3], error=function(cond) return(NA))
      
      if (!is.na(ptm)) {
        out.fit.coef <- out.fit$grid_fitted_coefficients
        psi.hat <- out.fit.coef[
          grepl(pattern=paste0("_lag_",lag_d),names(out.fit.coef)) &
            grepl(pattern=Treat_name,names(out.fit.coef))]
        psi.hat <- do.call(cbind,lapply(psi.hat, function(psi) 
          psi[,"estimate"]))
        
        # blip down outcome
        blip_down.t <- D.stacked.tvem[,rowSums(psi.hat*(.SD)), by=id, 
                                      .SDcols=colnames(psi.hat)][,V1]
        res.coef[[lag_d]] <- cbind("lag"=lag_d,
                                   "trt_time"=out.fit$time_grid,
                                   psi.hat)
        res.formulae[[lag_d]][["Outcome"]] <- list(tvc.formula,tfc.formula)
        res.fit[[lag_d]][["Outcome"]] <- out.fit
        rm(D.stacked.tvem)
      } else {
        # TVEM doesn't work: use linear interaction with time
        try.lm <- TRUE
      }
    }
    if (time.coef.fun=="interaction" || try.lm==TRUE) {
      out.formula <- as.formula(paste0(Outcome_name, "~", paste(
        c(paste0(out.preds,"*trt_time"),out.preds), collapse="+")))
      out.fit <- lm(out.formula, data=D.stacked)
      
      out.fit.coef <- coef(out.fit)
      psi.hat <- out.fit.coef[
        grepl(pattern=paste0("_lag_",lag_d),names(out.fit.coef)) &
          grepl(pattern=Treat_name,names(out.fit.coef))]
      res.coef[[lag_d]] <- psi.hat
      res.formulae[[lag_d]][["Outcome"]] <- formula(out.fit)
      res.fit[[lag_d]][["Outcome"]] <- out.fit
      
      # blip down outcome
      blip_down.t <- (model.matrix(out.fit)[,names(psi.hat),drop=FALSE] %*% 
                        psi.hat)[,1]
    }
    
    if (nrow(D.stacked)>length(blip_down.t)) {
      blip_down.t.all <- rep(NA,nrow(D.stacked))
      blip_down.t.all[as.integer(names(blip_down.t))] <- blip_down.t
      blip_down.t <- blip_down.t.all
      rm(blip_down.t.all)
    }
    Rt_minus1 <- D.stacked[,..Outcome_name] - blip_down.t
    D.stacked[, paste0("Ry.lag_",lag_d) := Rt_minus1]
    rm(Rt_minus1,blip_down.t)
    for (ss in lag_d:T_) {
      Rt_minus1 <- D.stacked[
        (trt_time+lag_d)==ss, c("id",paste0("Ry.lag_",lag_d)),with=FALSE]
      setnames(Rt_minus1,2,paste0("Ry",timesep,ss))
      setkey(Rt_minus1)
      if(colnames(Rt_minus1)[2] %in% colnames(D)) {
        D[,colnames(Rt_minus1)[2] := NULL]
        setkey(D)
      }
      D <- merge(D,Rt_minus1,all=TRUE)
      setkey(D)
      rm(Rt_minus1)
    }
  }
  
  res_lag <- list("est"=res.coef,
                  "formulae"=res.formulae,
                  "fitted"=res.fit)
  
  return( res_lag )
}

