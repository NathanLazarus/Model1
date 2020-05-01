setwd('C:/Users/Nathan/Downloads/PerturbationMethods/Model1')
#frenchdata = fread('C:/Users/Nathan/Downloads/Model1_Deprecated/Model1MultiplicativeUtility/FREDdata.csv')[,2]
frenchdata = fread('C:/Users/Nathan/Downloads/Model1_Deprecated/Model1MultiplicativeUtility/longtermproductivity.csv')
setnames(frenchdata,1,'TFP')
frenchdata[,logTFP:=log(TFP)]
frenchdata[,year:=.I+1889]

#frenchdata[,year:=.I+1953]

frenchdata_trend_pre1979 = lm(logTFP~year,frenchdata[year>1945&year<1980])
frenchdata[year>1945&year<1980,detrended:=logTFP - frenchdata_trend_pre1979$coefficients[2]*year - frenchdata_trend_pre1979$coefficients[1]]
frenchdata_trend_post1979 = lm(logTFP~year,frenchdata[year>1979])
frenchdata[year>1979,detrended:=logTFP - frenchdata_trend_post1979$coefficients[2]*year - frenchdata_trend_post1979$coefficients[1]]
frenchdata[,lagged_detrended := shift(detrended)]

frenchdata_trend_whole = lm(logTFP~year,frenchdata[year>1945])
summary(frenchdata_trend_whole)
frenchdata[year>1945,detrended_whole:=logTFP - frenchdata_trend_whole$coefficients[2]*year - frenchdata_trend_whole$coefficients[1]]
frenchdata[,lagged_detrended_whole := shift(detrended_whole)]
reg_whole = lm(detrended_whole~lagged_detrended_whole + 0,data = frenchdata[year>1945&year<1980])
summary(reg_whole)

# frenchdata = frenchdata[year>=1946]
# frenchdata[,shocks:=detrended-lagged_detrended*.9]
# g = exp(frenchdata_trend_whole$coefficients[2])
# frenchdata[year==1946,zeta:=1]
# for(i in which(frenchdata$year==1947):nrow(frenchdata)){
#   frenchdata[i,zeta:=frenchdata$zeta[i-1]^0.9*exp(frenchdata$shocks[i])]
# }

firstyear = 1980
frenchdata = frenchdata[year>=firstyear]
g2 = exp((frenchdata[year==max(year)]$logTFP - frenchdata[year==firstyear]$logTFP)/(max(frenchdata$year)-firstyear))
frenchdata[,fakeTFP := g2^(year-firstyear)*frenchdata[year==firstyear]$TFP]
frenchdata[,zeta:=TFP/fakeTFP]
frenchdata[,logzeta:=log(TFP)-log(fakeTFP)]
frenchdata[,shocks:=zeta^0.9/shift(zeta)]
frenchdata[,logshocks:=logzeta - shift(logzeta)*0.9]
frenchdata[,explogshocks:=exp(logshocks)]
frenchdata[year==firstyear,check:=1]
for(i in which(frenchdata$year==firstyear+1):nrow(frenchdata)){
  frenchdata[i,check:=frenchdata$check[i-1]^0.9*exp(frenchdata$logshocks[i])]
}
# frenchdata[,modeledTFPgrowth:=(g*zeta-shift(zeta,1))/sqrt(zeta*shift(zeta,1))]
# frenchdata[,realTFPgrowth:=(TFP-shift(TFP,1))/sqrt(TFP*shift(TFP,1))]
# sprintf("%.15f",g)
# sprintf("%.15f",g2)
# sprintf("%.15f",reg_whole$coefficients[1])
fwrite(frenchdata[,.(year,logshocks)],'TFPshocks.csv',col.names = F)
# 
# reg = lm(detrended_whole~lagged_detrended_whole + 0,data = frenchdata[year>1945&year<1980])
# summary(reg)
# reg2 = lm(detrended~lagged_detrended + 0,data = frenchdata[year>1980])
# summary(reg2)
# reg3 = lm(detrended~lagged_detrended + 0,data = frenchdata[year>1945&year!=1980])
# summary(reg3)
# 
# fwrite(frenchdata[year>0],file = 'FrenchDataRegressions.csv')