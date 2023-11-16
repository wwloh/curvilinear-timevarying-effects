rm(list=ls())
library("data.table")

filename <- "Most Recent Data - Released 5.10.2022 (Recommended)/2c_COVID19_daily_survey_ALL_cleaned_deid_2022-02-06_17_39.csv"
if(!file.exists(filename)) {
  library("osfr")
  # download data file from OSF
  osf_download(osf_ls_files(osf_retrieve_node("gpxwa")),conflicts="skip")
}

## prep data ##################################################################
rawData <- data.table::fread(filename)
setkey(rawData)
names(rawData)

# aggregate by week
rawData[, week_elapsed := ceiling((days_elapsed+1)/7)]

# demographic variables
rawData.demo <- data.table::fread("Most Recent Data - Released 5.10.2022 (Recommended)/1c_COVID19_demographics_cleaned_deid_2022-02-06_17_39.csv")
setkey(rawData.demo)
# define dummy variables for categorical responses
rawData.demo[, table(gender,useNA="ifany")]
rawData.demo[, gender_female := (gender==1)*1L]
rawData.demo[, table(marital,useNA="ifany")]
rawData.demo[, marital_single := (marital==1)*1L]
rawData.demo[, table(dependents,useNA="ifany")]
rawData.demo[, dependents := (dependents>0)*1L]
rawData.demo[, table(housing,useNA="ifany")]
rawData.demo[, housing1orfewer := (housing<=1)*1L]
rawData.demo[, table(income,useNA="ifany")]
rawData.demo[, income3orbelow := (income<=3)*1L]

rawData <- merge(rawData,rawData.demo,by="sub_id",all.x=TRUE)
setkey(rawData)
rm(rawData.demo)

# treatment
Xname <- "exercise"
# outcome
Yname <- "PHQ9"

# time-invariant covariates
Cnames <- c(
  "age1","gender_female","race1___3","marital_single","dependents",
  "housing1orfewer","income3orbelow"
)
# time-varying covariates
Lnames <- c(
  "sleepdiary_wakes"
)

Data <- rawData
setnames(Data,"sub_id","id")
setkey(Data)
select_cols <- c("id","week_elapsed",Cnames,Lnames,Xname,Yname)
Data <- Data[, ..select_cols]
setkey(Data)
rm(select_cols)

# aggregate by week
Data <- Data[!is.na(week_elapsed)] # remove missing times
Data <- Data[, lapply(.SD,mean,na.rm=TRUE), by=c("id","week_elapsed"), 
             .SDcols=c(Cnames,Lnames,Xname,Yname)]
Data[, median(exercise), by=id][,median(V1)]
Data[, exercise := (exercise >= 0.5)*1L]
setkey(Data)

# remove observations with missing data
Data <- Data[!apply(Data[, lapply(.SD,is.na),
                         .SDcols=c(Cnames,Lnames,Xname,Yname)],1,any)]
setkey(Data)

# sample size in each week
N.byweek <- Data[,.N,by=week_elapsed]
setkey(N.byweek)
N.byweek

# restrict to selected waves
T_ <- 13
Data <- Data[week_elapsed <= (T_+1L)]
setkey(Data)
setnames(Data,"week_elapsed","t")
setkey(Data)

# check any covariates with singular values
(SINGULAR <- which(Data[,lapply(.SD,function(x) length(unique(x)))]==1L))
if (length(SINGULAR)>0L) {
  Data[, colnames(Data)[SINGULAR] := NULL]  
}
setkey(Data)

# number of unique individuals
length(unique(Data$id))

summary(Data)

# inspect data for a randomly selected observation
Data[id==sample(id,1)]
