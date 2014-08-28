#####DETERMINE RESPONSE TYPES####

########ATTACH DATA###########

#attach table to use for ANOVAs
ANOVA= read.csv(file=file.choose())


#attach table for correlations
LinReg= read.csv(file=file.choose())


#attach necessary packages
library(nlme) #used for ANOVAs

#autocorrelation test
JulFl = read.csv(file=file.choose())
LeLn = read.csv(file=file.choose())
InLn = read.csv(file=file.choose())
Infl = read.csv(file=file.choose())

#####test for correlation between traits####

JulFl1 = tapply(JulFl[,7], JulFl[,2], data=JulFl, mean)
LeLn1 = tapply(LeLn[,7], LeLn[,2], data=LeLn, mean)
InLn1 = tapply(InLn[,7], InLn[,2], data=InLn, mean)
Infl1 = tapply(Infl[,7], Infl[,2], data=Infl, mean)

cor(JulFl1, Infl1)
cor(JulFl1, LeLn1)
cor(JulFl1, InLn1)
cor(LeLn1, InLn1)
cor(LeLn1, Infl1)
cor(InLn1, Infl1)

summary(lm(JulFl1~Infl1))
summary(lm(JulFl1~LeLn1))
summary(lm(JulFl1~InLn1))
summary(lm(LeLn1~InLn1))
summary(lm(LeLn1~Infl1))
summary(lm(InLn1~Infl1))




###########SHORT TERM############

#subset data for years less than 2001 this will include years 1994-2000
Short_TermANOVA <- subset(ANOVA, strYear < 2001)
Short_TermANOVA_C <- subset(ANOVA, strTrea=="C" & strYear < 2001) #want only control plots for this analysis
Short_TermLinReg_TDD<- subset(LinReg, strYear < 2001)
Short_TermLinReg_Year<- subset(LinReg, strTrea=="C" & strYear < 2001) #want only control plots for this analysis


#preview tables to make sure correct years and treatments are included
Short_TermANOVA$strYear
Short_TermANOVA_C$strTrea 
Short_TermLinReg_TDD$strYear
Short_TermLinReg_Year$strTrea

#determine number of levels for loop (should be the same for all)
str(Short_TermLinReg_TDD$SiteSppSexx)
str(Short_TermANOVA_C$SiteSppSexx)
str(Short_TermANOVA$SiteSppSexx)
str(Short_TermLinReg_Year$SiteSppSexx)

#create an empty data frame to send stats values to

N=28
Stats.ST <- data.frame(SiteSppSexx=numeric(N), linregTDD=numeric(N), ANOVA_C=numeric(N), linregYEAR=numeric(N), ANOVA_TREA=numeric(N), ANOVA_YEAR=numeric(N), ANOVA_INT=numeric(N), meanC=numeric(N), meanW=numeric(N))


AllSppLinReg_TDD <- unique(Short_TermLinReg_TDD$SiteSppSexx)
AllSppANOVA_C <- unique(Short_TermANOVA_C$SiteSppSexx)
AllSppANOVA <- unique(Short_TermANOVA$SiteSppSexx)
AllSppLinReg_Year <- unique(Short_TermLinReg_Year$SiteSppSexx)

for (i in 1:N) {
  SiteSppSexxLR.i <- AllSppLinReg_TDD[i]
  SiteSppSexxLR.ii <- AllSppLinReg_Year[i]
  SiteSppSexxANOVA.c <- AllSppANOVA_C[i]
  SiteSppSexxANOVA.i <- AllSppANOVA[i]
  
  ST_linreg.i <- Short_TermLinReg_TDD[Short_TermLinReg_TDD$SiteSppSexx == SiteSppSexxLR.i,]
  ST_linreg.ii <- Short_TermLinReg_Year[Short_TermLinReg_Year$SiteSppSexx == SiteSppSexxLR.i,]
  ST_anova.c <- Short_TermANOVA_C[Short_TermANOVA_C$SiteSppSexx == SiteSppSexxANOVA.c,]
  ST_anova.i <- Short_TermANOVA[Short_TermANOVA$SiteSppSexx == SiteSppSexxANOVA.i,]
  
  
  #correlation with TDD, first letter
  ST.linreg_TDD <- lm(numResponse ~ numTDDpsum, data=ST_linreg.i)
  
  
  #Relationship with year, second letter
  ST.anova_c<-lme(numResponse~factor(strYear), random=~1|factor(strIdLo), method="ML",data=ST_anova.c)  
  ST.linreg_Year <- lm(numResponse ~ strYear, data=ST_linreg.ii)
  
  
  #Response to experimental warming, third letter
  attach(ST_anova.i)
  ST.anova<-lme(numResponse~factor(strTrea)*factor(strYear), random=~1|factor(strIdLo), method="ML",data=ST_anova.i)
  
  #calculations
  means = tapply(numResponse, list(factor(strTrea)),data=ST_anova.i, mean)#means per group
  
  P.linregTDD=summary(ST.linreg_TDD)$coefficients[,4][2]
  P.anova_C = anova(ST.anova_c)$"p-value"[2]
  P.linregYear=summary(ST.linreg_Year)$coefficients[,4][2]
  P.anova_Trea=anova(ST.anova)$"p-value"[2]
  P.anova_Year=anova(ST.anova)$"p-value"[3]
  P.anova_Int=anova(ST.anova)$"p-value"[4]
  meanC=means[1]
  meanW=means[2]
  
  Stats.ST[i,] <- c(SiteSppSexxLR.i, P.linregTDD, P.anova_C, P.linregYear, P.anova_Trea, P.anova_Year, P.anova_Int, meanC, meanW)
  
  
  #print(SiteSppSexxLR.i)
  #print(summary(ST.linreg_TDD))
  #print(anova(ST.anova_c))
  #print(summary(ST.linreg_Year))
  #print(anova(ST.anova))
  #print(means)
  
}

Stats.STdata <- cbind(AllSppLinReg_TDD, Stats.ST)
####CHANGE FILE PATH WHEN CHANGING TRAITS####
write.csv(Stats.STdata, "c:/Users/Kelseyann/My Documents/2- Hollister Lab/Independent Study (2012)/QLeLn_ST.csv")




###########LONG TERM############

#subset data for years less than 2001 this will include years 1994-2000
Long_TermANOVA <- subset(ANOVA, strYear > 2001)
Long_TermANOVA_C <- subset(ANOVA, strTrea=="C" & strYear > 2001) #want only control plots for this analysis
Long_TermLinReg_TDD<- subset(LinReg, strYear > 2001)
Long_TermLinReg_Year<- subset(LinReg, strTrea=="C" & strYear > 2001) #want only control plots for this analysis


#preview tables to make sure correct years and treatments are included
Long_TermANOVA$strYear
Long_TermANOVA_C$strTrea 
Long_TermLinReg_TDD$strYear
Long_TermLinReg_Year$strTrea

#determine number of levels for loop (should be the same for all)
str(Long_TermLinReg_TDD$SiteSppSexx)
str(Long_TermANOVA_C$SiteSppSexx)
str(Long_TermANOVA$SiteSppSexx)
str(Long_TermLinReg_Year$SiteSppSexx)


#create empty data frame to send results to
N=28
Stats.LT <- data.frame(SiteSppSexx=numeric(N), linregTDD=numeric(N), ANOVA_C=numeric(N), linregYEAR=numeric(N), ANOVA_TREA=numeric(N), ANOVA_YEAR=numeric(N), ANOVA_INT=numeric(N), meanC=numeric(N), meanW=numeric(N))


AllSppLinReg_TDD <- unique(Long_TermLinReg_TDD$SiteSppSexx)
AllSppANOVA_C <- unique(Long_TermANOVA_C$SiteSppSexx)
AllSppANOVA <- unique(Long_TermANOVA$SiteSppSexx)
AllSppLinReg_Year <- unique(Long_TermLinReg_Year$SiteSppSexx)

for (i in 1:N) {
  SiteSppSexxLR.i <- AllSppLinReg_TDD[i]
  SiteSppSexxLR.ii <- AllSppLinReg_Year[i]
  SiteSppSexxANOVA.c <- AllSppANOVA_C[i]
  SiteSppSexxANOVA.i <- AllSppANOVA[i]
  
  LT_linreg.i <- Long_TermLinReg_TDD[Long_TermLinReg_TDD$SiteSppSexx == SiteSppSexxLR.i,]
  LT_linreg.ii <- Long_TermLinReg_Year[Long_TermLinReg_Year$SiteSppSexx == SiteSppSexxLR.i,]
  LT_anova.c <- Long_TermANOVA_C[Long_TermANOVA_C$SiteSppSexx == SiteSppSexxANOVA.c,]
  LT_anova.i <- Long_TermANOVA[Long_TermANOVA$SiteSppSexx == SiteSppSexxANOVA.i,]
  
  
  #correlation with TDD, first letter
  LT.linreg_TDD <- lm(numResponse ~ numTDDpsum, data=LT_linreg.i)
  
  
  #Relationship with year, second letter
  LT.anova_c<-lme(numResponse~factor(strYear), random=~1|factor(strIdLo), method="ML",data=LT_anova.c)  
  LT.linreg_Year <- lm(numResponse ~ strYear, data=LT_linreg.ii)
  
  
  #Response to experimental warming, third letter
  attach(LT_anova.i)
  LT.anova<-lme(numResponse~factor(strTrea)*factor(strYear), random=~1|factor(strIdLo), method="ML",data=LT_anova.i)
  
  #calculations
  means = tapply(numResponse, list(factor(strTrea)),data=LT_anova.i, mean)#means per group
  
  
  P.linregTDD=summary(LT.linreg_TDD)$coefficients[,4][2]
  P.anova_C = anova(LT.anova_c)$"p-value"[2]
  P.linregYear=summary(LT.linreg_Year)$coefficients[,4][2]
  P.anova_Trea=anova(LT.anova)$"p-value"[2]
  P.anova_Year=anova(LT.anova)$"p-value"[3]
  P.anova_Int=anova(LT.anova)$"p-value"[4]
  meanC=means[1]
  meanW=means[2]
  
  Stats.LT[i,] <- c(SiteSppSexxLR.i, P.linregTDD, P.anova_C, P.linregYear, P.anova_Trea, P.anova_Year, P.anova_Int, meanC, meanW)
  
  
  #print(SiteSppSexxLR.i)
  #print(summary(LT.linreg_TDD))
  #print(anova(LT.anova_c))
  #print(summary(LT.linreg_Year))
  #print(anova(LT.anova))
  #print(means)
  
}

Stats.LTdata <- cbind(AllSppLinReg_TDD, Stats.LT)
####CHANGE FILE PATH WHEN CHANGING TRAITS####
write.csv(Stats.LTdata, "c:/Users/Kelseyann/My Documents/2- Hollister Lab/Independent Study (2012)/QLeLn_LT.csv")








###########BOXPLOTS############

#attach table to use for plots
Plot_Infl= read.csv(file=file.choose())
Plot_InLn= read.csv(file=file.choose())
Plot_LeLn= read.csv(file=file.choose())
Plot_JulFl= read.csv(file=file.choose())

str(Plot_Infl$SiteSppSexx)
str(Plot_InLn$SiteSppSexx)
str(Plot_LeLn$SiteSppSexx)
str(Plot_JulFl$SiteSppSexx)

#Infl
pdf("c:/Users/Kelseyann/My Documents/ALASKA/Independent Study (2012)/Consistency All Species/R Outputs/QInfl_graph.pdf")

for (i in 1:22) {
  AllSpp <- unique(Plot_Infl$SiteSppSexx)
  SiteSppSexx.i <- AllSpp[i]
  Species.i <- Plot_Infl[Plot_Infl$SiteSppSexx == SiteSppSexx.i,]
  
  boxplot(numResponse~strYear, main=SiteSppSexx.i, data=Species.i)
  
}
dev.off()

#InLn
pdf("c:/Users/Kelseyann/My Documents/ALASKA/Independent Study (2012)/Consistency All Species/R Outputs/QInLn_graph.pdf")

for (i in 1:21) {
  AllSpp <- unique(Plot_InLn$SiteSppSexx)
  SiteSppSexx.i <- AllSpp[i]
  Species.i <- Plot_InLn[Plot_InLn$SiteSppSexx == SiteSppSexx.i,]
  
  boxplot(numResponse~strYear, main=SiteSppSexx.i, data=Species.i)
  
}
dev.off()

#LeLn
pdf("c:/Users/Kelseyann/My Documents/ALASKA/Independent Study (2012)/Consistency All Species/R Outputs/QLeLn_graph.pdf")

for (i in 1:28) {
  AllSpp <- unique(Plot_LeLn$SiteSppSexx)
  SiteSppSexx.i <- AllSpp[i]
  Species.i <- Plot_LeLn[Plot_LeLn$SiteSppSexx == SiteSppSexx.i,]
  
  boxplot(numResponse~strYear, main=SiteSppSexx.i, data=Species.i)
  
}
dev.off()

#JulFl
pdf("c:/Users/Kelseyann/My Documents/ALASKA/Independent Study (2012)/Consistency All Species/R Outputs/JulFl_graph.pdf")

for (i in 1:19) {
  AllSpp <- unique(Plot_JulFl$SiteSppSexx)
  SiteSppSexx.i <- AllSpp[i]
  Species.i <- Plot_JulFl[Plot_JulFl$SiteSppSexx == SiteSppSexx.i,]
  
  boxplot(numResponse~strYear, main=SiteSppSexx.i, data=Species.i)
  
}
dev.off()

