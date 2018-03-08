dat <- read.csv("Merged body comp and phenotype data 20170425.csv")
dat <- dat[1:211,]
source("RFDist.R")
source("PAM.R")
source("AdjRand.R")
library(MASS)

no.forests <- 100
no.trees <- 4000
distRF <- RFdist(dat[,-c(1,2,3)], mtry1=10, no.trees, no.forests, addcl1=T,addcl2=T,imp=T, oob.prox1=T)

cmd1 <- cmdscale(as.dist(distRF$cl1),2)
cmd2 <- cmdscale(as.dist(distRF$cl2),2)
iso1 <- isoMDS(as.dist(distRF$cl1),k=2)$points
iso2 <- isoMDS(as.dist(distRF$cl2),k=2)$points

dat$label1 <- pamNew(dist(distRF$cl1),2)
dat$label2 <- pamNew(dist(distRF$cl1),3)
dat$label3 <- pamNew(dist(distRF$cl1),4)
dat$label4 <- pamNew(dist(distRF$cl1),5)
dat$label5 <- pamNew(dist(distRF$cl1),6)
dat$label6 <- pamNew(dist(distRF$cl1),7)
dat$label7 <- pamNew(dist(distRF$cl1),8)

par(mfrow=c(1,1))
plot(distRF$imp1[,4]  ,xlab="variables",ylab="Gini Index Measure", main="Addcl1, Variable Imp.")
distRF$imp1[,4][distRF$imp1[,4]>2.5]


library(sas7bdat)
clinical_data1 <- read.sas7bdat("decamp1cluster.sas7bdat")

case_name <- sapply(dat$Case, as.character)
split1 <- function(x){
  strsplit(x, "-")
}
converted_case <- matrix(unlist(split1(case_name)),ncol=2,byrow = T)
dat$study <- converted_case[, 1]
dat$cn <- converted_case[, 2]

clinical_data1 <- clinical_data1[clinical_data1$cn %in% dat[which(dat$study=="4703"),]$cn,] 
clinical_data1$label1 <- dat$label1
clinical_data1$label2 <- dat$label2
clinical_data1$label3 <- dat$label3


k.max <- 8
data <- dat[,-c(2,3,102,103,104,105,106,107,108)]
sil <- rep(0, k.max)
# Compute the average silhouette width for 
# k = 2 to k = 7
total_label <- cbind(rep(0,211),dat$label1, dat$label2, dat$label3, dat$label4,
                     dat$label5, dat$label6, dat$label7)
for(i in 2:k.max){
  ss <- silhouette(total_label[,i], dist(distRF$cl1))
  sil[i] <- mean(ss[, 3])
}
# Plot the  average silhouette width
par(mfrow=(c(1,1)))
plot(1:k.max, sil, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of clusters k", ylab = "Average silhouette width",
     main = "Optimal number of clusters", ylim = c(0,0.3))
abline(v = which.max(sil), lty = 2)


library(fpc)
plotcluster(data[,-c(100,101)], dat$label1) # discriminant projection plot
plotcluster(data[,-c(100,101)], dat$label2) # discriminant projection plot
# 2 cluster
plot(cmd1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
points(cmd1, pch = ifelse(dat$label1==1, 1, 3),
       col = ifelse(dat$label1==1,"red","blue"))

plot(iso1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
points(iso1, pch = ifelse(dat$label1==1, 1, 3),
       col = ifelse(dat$label1==1,"red","blue"))

# 3 cluster
plot(cmd1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
points(cmd1, pch = ifelse(dat$label2==1, 3, ifelse(dat$label2==2,1,2)),
       col = ifelse(dat$label2==1, 3, ifelse(dat$label2==2,1,2)))

plot(iso1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
points(iso1, pch = ifelse(dat$label2==1, 3, ifelse(dat$label2==2,1,2)),
       col = ifelse(dat$label2==1, 3, ifelse(dat$label2==2,1,2)))

clinical_data <- data.frame(AGE=clinical_data1$AGE)
clinical_data$label1 <- clinical_data1$label1
clinical_data$label2 <- clinical_data1$label2
clinical_data$BMI <- clinical_data1$basebmi
clinical_data$FEV1 <- clinical_data1$FEV1pred
clinical_data$PKyears <- clinical_data1$packyears
clinical_data$race <- clinical_data1$race
clinical_data$ethnicity <- clinical_data1$ethnicity
clinical_data$sex <- clinical_data1$sex
clinical_data$smkostat <- clinical_data1$ptsmokstat
clinical_data$calcCOPD <- clinical_data1$calcCOPD

mean1 <- function(x){
  mean(x, na.rm = T)
}
sd1 <- function(x){
  sd(x, na.rm = T)
}


## 2 cluster
table(clinical_data$label1)
mean1(clinical_data$AGE[which(clinical_data$label1==1)])
sd1(clinical_data$AGE[which(clinical_data$label1==1)])
mean1(clinical_data$AGE[which(clinical_data$label1==2)])
sd1(clinical_data$AGE[which(clinical_data$label1==2)])
sum(is.na(clinical_data$AGE[which(clinical_data$label1==1)]))
sum(is.na(clinical_data$AGE[which(clinical_data$label1==2)]))


mean1(clinical_data$BMI[which(clinical_data$label1==1)])
sd1(clinical_data$BMI[which(clinical_data$label1==1)])
mean1(clinical_data$BMI[which(clinical_data$label1==2)])
sd1(clinical_data$BMI[which(clinical_data$label1==2)])
sum(is.na(clinical_data$BMI[which(clinical_data$label1==1)]))
sum(is.na(clinical_data$BMI[which(clinical_data$label1==2)]))


mean1(clinical_data$FEV1[which(clinical_data$label1==1)])
sd1(clinical_data$FEV1[which(clinical_data$label1==1)])
mean1(clinical_data$FEV1[which(clinical_data$label1==2)])
sd1(clinical_data$FEV1[which(clinical_data$label1==2)])
sum(is.na(clinical_data$FEV1[which(clinical_data$label1==1)]))
sum(is.na(clinical_data$FEV1[which(clinical_data$label1==2)]))


mean1(clinical_data$PKyears[which(clinical_data$label1==1)])
sd1(clinical_data$PKyears[which(clinical_data$label1==1)])
mean1(clinical_data$PKyears[which(clinical_data$label1==2)])
sd1(clinical_data$PKyears[which(clinical_data$label1==2)])
sum(is.na(clinical_data$PKyears[which(clinical_data$label1==1)]))
sum(is.na(clinical_data$PKyears[which(clinical_data$label1==2)]))

table(clinical_data$sex,clinical_data$label1)
table(clinical_data$race,clinical_data$label1)
table(clinical_data$ethnicity,clinical_data$label1)
table(clinical_data$smkostat,clinical_data$label1)
sum(is.na(clinical_data$smkostat[which(clinical_data$label1==1)]))
sum(is.na(clinical_data$smkostat[which(clinical_data$label1==2)]))
table(clinical_data$calcCOPD,clinical_data$label1)
sum(is.na(clinical_data$calcCOPD[which(clinical_data$label1==1)]))
sum(is.na(clinical_data$calcCOPD[which(clinical_data$label1==2)]))

####3 cluster
table(clinical_data$label2)
mean1(clinical_data$AGE[which(clinical_data$label2==1)])
sd1(clinical_data$AGE[which(clinical_data$label2==1)])
mean1(clinical_data$AGE[which(clinical_data$label2==2)])
sd1(clinical_data$AGE[which(clinical_data$label2==2)])
mean1(clinical_data$AGE[which(clinical_data$label2==3)])
sd1(clinical_data$AGE[which(clinical_data$label2==3)])
sum(is.na(clinical_data$AGE[which(clinical_data$label2==1)]))
sum(is.na(clinical_data$AGE[which(clinical_data$label2==2)]))
sum(is.na(clinical_data$AGE[which(clinical_data$label2==3)]))


mean1(clinical_data$BMI[which(clinical_data$label2==1)])
sd1(clinical_data$BMI[which(clinical_data$label2==1)])
mean1(clinical_data$BMI[which(clinical_data$label2==2)])
sd1(clinical_data$BMI[which(clinical_data$label2==2)])
mean1(clinical_data$BMI[which(clinical_data$label2==3)])
sd1(clinical_data$BMI[which(clinical_data$label2==3)])
sum(is.na(clinical_data$BMI[which(clinical_data$label2==1)]))
sum(is.na(clinical_data$BMI[which(clinical_data$label2==2)]))
sum(is.na(clinical_data$BMI[which(clinical_data$label2==3)]))


mean1(clinical_data$FEV1[which(clinical_data$label2==1)])
sd1(clinical_data$FEV1[which(clinical_data$label2==1)])
mean1(clinical_data$FEV1[which(clinical_data$label2==2)])
sd1(clinical_data$FEV1[which(clinical_data$label2==2)])
mean1(clinical_data$FEV1[which(clinical_data$label2==3)])
sd1(clinical_data$FEV1[which(clinical_data$label2==3)])
sum(is.na(clinical_data$FEV1[which(clinical_data$label2==1)]))
sum(is.na(clinical_data$FEV1[which(clinical_data$label2==2)]))
sum(is.na(clinical_data$FEV1[which(clinical_data$label2==3)]))


mean1(clinical_data$PKyears[which(clinical_data$label2==1)])
sd1(clinical_data$PKyears[which(clinical_data$label2==1)])
mean1(clinical_data$PKyears[which(clinical_data$label2==2)])
sd1(clinical_data$PKyears[which(clinical_data$label2==2)])
mean1(clinical_data$PKyears[which(clinical_data$label2==3)])
sd1(clinical_data$PKyears[which(clinical_data$label2==3)])
sum(is.na(clinical_data$PKyears[which(clinical_data$label2==1)]))
sum(is.na(clinical_data$PKyears[which(clinical_data$label2==2)]))
sum(is.na(clinical_data$PKyears[which(clinical_data$label2==3)]))

table(clinical_data$sex,clinical_data$label2)
table(clinical_data$race,clinical_data$label2)
table(clinical_data$ethnicity,clinical_data$label2)
table(clinical_data$smkostat,clinical_data$label2)
sum(is.na(clinical_data$smkostat[which(clinical_data$label2==1)]))
sum(is.na(clinical_data$smkostat[which(clinical_data$label2==2)]))
sum(is.na(clinical_data$smkostat[which(clinical_data$label2==3)]))

table(clinical_data$calcCOPD,clinical_data$label2)
sum(is.na(clinical_data$calcCOPD[which(clinical_data$label2==1)]))
sum(is.na(clinical_data$calcCOPD[which(clinical_data$label2==2)]))
sum(is.na(clinical_data$calcCOPD[which(clinical_data$label2==3)]))

table(clinical_data$label1, clinical_data$label2)

t.test(clinical_data$AGE[which(clinical_data$label1==1)],
       clinical_data$AGE[which(clinical_data$label1==2)])
t.test(clinical_data$BMI[which(clinical_data$label1==1)],
       clinical_data$BMI[which(clinical_data$label1==2)])
t.test(clinical_data$FEV1[which(clinical_data$label1==1)],
       clinical_data$FEV1[which(clinical_data$label1==2)])
t.test(clinical_data$PKyears[which(clinical_data$label1==1)],
       clinical_data$PKyears[which(clinical_data$label1==2)])

summary(aov(clinical_data$AGE~clinical_data$label2))
summary(aov(clinical_data$BMI~clinical_data$label2))
summary(aov(clinical_data$FEV1~clinical_data$label2))
summary(aov(clinical_data$PKyears~clinical_data$label2))

fisher.test(table(clinical_data$sex, clinical_data$label1))
fisher.test(table(clinical_data$race, clinical_data$label1))
fisher.test(table(clinical_data$ethnicity, clinical_data$label1))
fisher.test(table(clinical_data$smkostat, clinical_data$label1))
fisher.test(table(clinical_data$calcCOPD, clinical_data$label1))

fisher.test(table(clinical_data$sex, clinical_data$label2))
x<-matrix(c(80,19,10,114,33,18,101,15,13),ncol=3)
fisher.test(x,workspace = 2e7)

fisher.test(table(clinical_data$race, clinical_data$label2))
fisher.test(table(clinical_data$ethnicity, clinical_data$label2))
fisher.test(table(clinical_data$smkostat, clinical_data$label2))
fisher.test(table(clinical_data$calcCOPD, clinical_data$label2))
