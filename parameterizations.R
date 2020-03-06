setwd("C:/Users/Nathan/Downloads/Profits")
data = fread('paneldata.csv')
library(foreach)
library(iterators)
library(snow)
library(doSNOW)
library(readxl)
tfp = data.table(read_excel('C:/Users/Nathan/Downloads/BCLDatabase_online_v2.3.xlsx',sheet='TFP'))
setnames(tfp,1,'Year')
tfp[,Year:=as.integer(Year)]
data[tfp,on = 'Year',tfp:=i.USA]
data[,tfpgrowth:=tfp/shift(tfp)-1]

setnames(data,
         c("PtntppPatentApplicationPopula","CprtppCopyrightRegisteredPopu","TrademarksRegisteredPopulation","privateRD_GDP","publicRD_GDP","ProductivityNonFarmOutputper","GrowthrateGrowthRateofRealG"),
         c("patents","copyrights","trademarks","privateRD","publicRD","productivity","GDP"))
setnames(data,
         c("HouseMedianMedianideologyHouse","SenateMedianmedianideologySena","IdeologyPresmedianideologyPresid","DefenseGDPDefenseExpenditures","mergers_GDP","PolicyUncertainty"),
         c("house","senate","pres","defense","mergers","uncertainty"))
setnames(data,
         c("nonDefenseGDP","UnionmembersofunionNonFarm","WorkersinvolvedinStopagesas","CorptaxCorporateTaxRateofhi","MaxtaxMaximalIndividualMargin","MinWagebiteadjustedforcovera","BOPGDP","CivServPop","PagesofRegulations","AntitrustCivilCases","RestraintofTradecases"),
         c("nondefense","union","stoppages","corptax","maxtax","minwage","BOP","civserv","regs","antitrust","restraintoftrade"))

setnames(data,"Adjustedfouryearmovingaverag","inflation")

data[,logprofits := log(pmax(ProfitShare,0.01))]
data[,logLAGprofitshigh := log(pmax(Profitshare1,0.01))*((Year>1888.5&Year<1901.5)|(Year>1984.5&Year<2017.5))]
data[,logLAGprofitslow := log(pmax(Profitshare1,0.01))*(Year>1931.5&Year<1984.5)]
data[,logLAGprofits := log(pmax(Profitshare1,0.01))]
data[,productivity1:=shift(tfpgrowth,1)]
data[,productivity2:=shift(tfpgrowth,2)]
data[,productivity3:=shift(tfpgrowth,3)]
data[,productivity4:=shift(tfpgrowth,4)]
data[,productivity5:=shift(tfpgrowth,5)]
parameters = lm(logprofits~logLAGprofits+productivity1+productivity2+productivity3+productivity4+productivity5,data=data[Year>1889.5&Year<2017.5])
parameters1 = lm(logprofits~logLAGprofitshigh+productivity1+productivity2+productivity3+productivity4+productivity5,data=data[((Year>1888.5&Year<1901.5)|(Year>1984.5&Year<2017.5))])
parameters2 = lm(logprofits~logLAGprofitslow+productivity1+productivity2+productivity3+productivity4+productivity5,data=data[(Year>1931.5&Year<1984.5)])
test = nls(logprofits~a1*logLAGprofitshigh+a3*logLAGprofitslow+a2*productivity1+a2^2*productivity2+a2^3*productivity3+a2^4*productivity4+a2^5*productivity5,data=data[Year>1889.5&Year<2017.5], start = list(a1 = 0, a2 = 0,a3=0), trace = T,algorithm = 'port',lower=c(-1000,0,-1000), upper=c(1000,1000,1000))
summary(parameters)
summary(parameters1)
summary(parameters2)
summary(test)