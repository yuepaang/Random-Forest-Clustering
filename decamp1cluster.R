library(sas7bdat)
dat <- read.sas7bdat("decamp1cluster.sas7bdat")
status_var <- c("sex", "ethnicity", "race", "education", "marital", "income", "deploy")
current_military_status<- c("active", "retired", "veteran", "fammember")
cough_var <- c(dat$cough, dat$cough46, dat$coughmorning, dat$coughnight,dat$cough3mos,
               dat$yearscough)
symptom_var <- c(dat$shortnessbreath,dat$walkslower,dat$stopwalking100,dat$toobreathless,dat$dyspnea,
                 dat$yearsshortnessbreath,dat$numexacerbations,dat$numexacerbationsadmission,
                 dat$chroniccough)
phlegm_var <- c(dat$phlegm,dat$phlegm2,dat$phlegmmorning,dat$phlegmnight,dat$phlegm3mos,
                dat$yearsphlegm,dat$chronicphlegm)
medical_history_var <- c(dat$lcfather, dat$lcmother,dat$lcbrother,dat$lcsister,dat$lcchildren,
                         dat$priorcancer,dat$priorcancerYN)

na_count <- function(x){
  sum(is.na(x))
}
# counting the missing value each columns
col_na <- apply(dat,2,na_count)
dat1 <- dat[,names(col_na[which(col_na < 100)])]
# counting the missing value each rows
row_na <- apply(dat1,1,na_count)
dat2 <- dat1[which(row_na < 30),]
dat3 <- dat2[complete.cases(dat2),]

# dataset to analyze
dat4 <- dat3[,-c(6,7,8,9,22,23,24,25)]  # smoking_var <- c("CFsmoke", "PKyears")
dat5 <- dat4[,-c(1,2,3,4)]

model1 <- princomp(dat5)
summary(model1)
loadings(model1)
model1$scores

library(psych)
# remove "active", "cysfib", "ipulmonaryfib", "obbronch", "centralairobstruct", "TBpneumonia", "pft", "study"
dat6 <- dat5[,-c(9,46,48,49,50,51,53,68)]
fit <- principal(dat6, nfactors=10) # 59%
fit$Structure
fit$scores

fit1 <- principal(dat6,13,rotate = "Promax") # 68%
fit2 <-fa(dat6,10,n.obs = 71, fm = "minrank", rotate = "Promax") #72%

dat_reduction <- fit1$scores

source("RFDist.R")
source("PAM.R")
source("AdjRand.R")

no.forests <- 100
no.trees <- 4000
distRF <- RFdist(dat_reduction, mtry1=3, no.trees, no.forests, addcl1=T,addcl2=T,imp=T, oob.prox1=T)

cmd1 <- cmdscale(as.dist(distRF$cl1),2)
cmd2 <- cmdscale(as.dist(distRF$cl2),2)
iso1 <- isoMDS(as.dist(distRF$cl1),k=2)$points
iso2 <- isoMDS(as.dist(distRF$cl2),k=2)$points

plot(cmd1, xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
plot(cmd2, xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2") # better
plot(iso1, xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
plot(iso2, xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")

RFclusterLabel <- pamNew(cmd2, 3)
RFclusterLabel2 <- pamNew(as.dist(distRF$cl1), 3)
adjustedRandIndex(RFclusterLabel, RFclusterLabel2)
# After clustering analysis
dat6$label <- RFclusterLabel
dat6$FEV1FVC[17]<- round(2.69/3.48,digits = 2)
dat6$FEV1FVC[32]<-round(1.10/3.72,digits = 2)
dat6$FEV1FVC[59]<-round(1.13/1.63,digits = 2)
dat6$FEV1FVC[64]<- round(2.17/3.54,digits = 2)

n1 = sum(dat6$label==1)
n2 = sum(dat6$label==2)
n3 = sum(dat6$label==3)

# continuous variables
col_idx <- c(1,2,37,38,39,40,48,49,50,51,52,53,54,55,56,62)
df1 <- data.frame(variables = names(dat6)[col_idx], label1 = rep(0,length(col_idx)), 
                  label2 = rep(0,length(col_idx)), label3 = rep(0,length(col_idx)))
for (i in 1:length(df1$variables)){
  
  df1$label1[i] <- paste(round(by(dat6[,names(dat6)[col_idx][i]],dat6$label,mean)[1], digits = 1), "(",
                         round(by(dat6[,names(dat6)[col_idx][i]],dat6$label,sd)[1], digits = 1), ")",sep = "")
  df1$label2[i] <- paste(round(by(dat6[,names(dat6)[col_idx][i]],dat6$label,mean)[2], digits = 1), "(",
                         round(by(dat6[,names(dat6)[col_idx][i]],dat6$label,sd)[2], digits = 1), ")", sep = "")
  df1$label3[i] <- paste(round(by(dat6[,names(dat6)[col_idx][i]],dat6$label,mean)[3], digits = 1), "(",
                         round(by(dat6[,names(dat6)[col_idx][i]],dat6$label,sd)[3], digits = 1),")", sep = "")
}

table(dat6$cough[which(dat6$label==1)])[2]/n1 # 0.5 highest
table(dat6$cough[which(dat6$label==2)])[2]/n2 # higher
table(dat6$cough[which(dat6$label==3)])[2]/n3


table(dat6$dyspnea[which(dat6$label==1)])
table(dat6$dyspnea[which(dat6$label==2)]) # mostly yes
table(dat6$dyspnea[which(dat6$label==3)]) # mostly no

table(dat6$chroniccough[which(dat6$label==1)]) # mostly unknown
table(dat6$chroniccough[which(dat6$label==2)]) # mostly No 0.71
table(dat6$chroniccough[which(dat6$label==3)])

table(dat6$chronicphlegm[which(dat6$label==1)])
table(dat6$chronicphlegm[which(dat6$label==2)]) # mostly No 0.75
table(dat6$chronicphlegm[which(dat6$label==3)])

table(dat6$chronicbronchitis[which(dat6$label==1)]) # all unknow
table(dat6$chronicbronchitis[which(dat6$label==2)]) #mostly No 0.83
table(dat6$chronicbronchitis[which(dat6$label==3)]) 

# Interstitial lung disease 
dat6$ILDpt[which(dat6$label==1)]
dat6$ILDpt[which(dat6$label==2)]
dat6$ILDpt[which(dat6$label==3)] # has YES 



sum(dat6$calcCOPD[which(dat6$label==1)]==1)/n1
sum(dat6$calcCOPD[which(dat6$label==2)]==1)/n2
sum(dat6$calcCOPD[which(dat6$label==3)]==1)/n3

sum(dat6$calcCOPDlln[which(dat6$label==1)]==1)/n1
sum(dat6$calcCOPDlln[which(dat6$label==2)]==1)/n2
sum(dat6$calcCOPDlln[which(dat6$label==3)]==1)/n3





# threshold rule
library(rpart)
rp1 = rpart(factor(label) ~ ., dat6[,c(names(dat6)[col_idx],"label")])
plot(rp1, uniform=T, branch=0, margin=0.1, main = "Classification Tree")
text(rp1, all=T, use.n=T, cex=0.7)
summary(rp1)

rp2 <- rpart(factor(label)~., dat6[dat6$AGE<=76.5,c(names(dat6)[col_idx],"label")])
plot(rp2, uniform=T, branch=0, margin=0.1, main = "Classification Tree")
text(rp2, all=T, use.n=T, cex=0.7)
summary(rp2)

# label 1
dat6$label[dat6$AGE<76.5 & dat6$chronicbronchitis>2] # 70%
# label 2
dat6$label[dat6$AGE>=76.5|(dat6$AGE<=76.5 & dat6$COPDlln >= 0.6589)] # 70.83%
# label 3
dat6$label[dat6$AGE<76.5 & dat6$COPDlln < 0.6589] #83.8%
