# Analysis
dat <- read.csv("data.csv",stringsAsFactors=FALSE)
dat <- na.omit(dat)
source("RFDist.R")
source("PAM.R")
source("AdjRand.R")
library(MASS)
no.forests <- 100
no.trees <- 4000
train_data <- dat[,-c(1,2)]
distRF <- RFdist(train_data, mtry1=10, no.trees, no.forests, addcl1=T,addcl2=T,imp=T, oob.prox1=T)

cmd1 <- cmdscale(as.dist(distRF$cl1),2)
cmd2 <- cmdscale(as.dist(distRF$cl2),2)
iso1 <- isoMDS(as.dist(distRF$cl1),k=2)$points
iso2 <- isoMDS(as.dist(distRF$cl2),k=2)$points

par(mfrow=c(2,3))
plot(cmd1, xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2", main = c("Addcl1, cMDS"))
plot(iso1, xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2", main = c("Addcl1, isoMDS"))
# better
plot(distRF$imp1[,4], xlab = "variables", ylab = "Gini Index Measure",
     main = "Addcl1, Variable Imp.")
plot(cmd2, xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2", main = c("Addcl2, cMDS"))
plot(iso2, xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2",main = c("Addcl2, isoMDS"))
plot(distRF$imp2[,4], xlab = "variables", ylab = "Gini Index Measure",
     main = "Addcl2, Variable Imp.")


RFclusterLabel <- kmeans(dist(distRF$cl2,"manhattan"),2)$cluster
RFclusterLabel2 <- kmeans(dist(distRF$cl2,"manhattan"),3)$cluster
RFclusterLabel3 <- kmeans(dist(distRF$cl2,"manhattan"),4)$cluster
dat$label1 <- RFclusterLabel
dat$label2 <- RFclusterLabel2
dat$label3 <- RFclusterLabel3

dat$label1 <- pamNew(dist(distRF$cl1),2)
dat$label2 <- pamNew(dist(distRF$cl1),3)
dat$label3 <- pamNew(dist(distRF$cl1),4)
dat$label4 <- pamNew(dist(distRF$cl1),5)
dat$label5 <- pamNew(dist(distRF$cl1),6)
dat$label6 <- pamNew(dist(distRF$cl1),7)
dat$label7 <- pamNew(dist(distRF$cl1),8)

par(mfrow=c(1,1))
plot(distRF$imp1[,4]  ,xlab="variables",ylab="Gini Index Measure", main="Addcl1, Variable Imp.")

#
# dat$label1 <- pam(dist(distRF$cl1),2)$cluster
# dat$label2 <- pam(dist(distRF$cl1),3)$cluster
# dat$label3 <- pam(dist(distRF$cl1),4)$cluster
# dat$label4 <- pam(dist(distRF$cl1),5)$cluster
# dat$label5 <- pam(dist(distRF$cl1),6)$cluster
# dat$label6 <- pam(dist(distRF$cl1),7)$cluster
# dat$label7 <- pam(dist(distRF$cl1),8)$cluster

# After clustering
# extract the case
case_name <- sapply(dat$Case, as.character)
split1 <- function(x){
  strsplit(x, "-")
}
converted_case <- matrix(unlist(split1(case_name)),ncol=2,byrow = T)
dat$study <- converted_case[, 1]
dat$cn <- converted_case[, 2]
library(sas7bdat)
clinical_data1 <- read.sas7bdat("data.sas7bdat")
clinical_data2 <- read.sas7bdat("data.sas7bdat")
dat1 <- dat[which(dat$study=="A"),]
dat2 <- dat[which(dat$study=="B"),]

clinical_data1 <- clinical_data1[clinical_data1$cn %in% dat[which(dat$study=="A"),]$cn,] 
clinical_data1$label1 <- dat1$label1
clinical_data1$label2 <- dat1$label2
clinical_data1$label3 <- dat1$label3



clinical_data2 <- clinical_data2[clinical_data2$cn %in% dat[which(dat$study=="B"),]$cn,] 
clinical_data2$label1 <- dat2$label1
clinical_data2$label2 <- dat2$label2
clinical_data2$label3 <- dat2$label3



# visualization
k.max <- 8
data <- dat[,-c(2,3,102,103,104,105,106,107,108)]
sil <- rep(0, k.max)
# Compute the average silhouette width for 
# k = 2 to k = 7
total_label <- cbind(rep(0,414),dat$label1, dat$label2, dat$label3, dat$label4,
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


library(ggplot2)
library(ggfortify)
autoplot(prcomp(data[,-c(100,101)]), data=dat,colour="label1") + scale_y_continuous(limits = c(-0.2,0.3))
autoplot(prcomp(data[,-c(100,101)]), data=dat,colour="label2") + scale_y_continuous(limits = c(-0.2,0.3))
autoplot(prcomp(data[,-c(100,101)]), data=dat,colour="label3") 

autoplot(prcomp(dist(distRF$cl1)), data=dat,colour="label1") 
autoplot(prcomp(dist(distRF$cl1)), data=dat,colour="label2") 
#####
pcs<- princomp(data[,-c(100,101)],scores = T)
pc <- prcomp(data[,-c(100,101)],scale. = T)
plot(pc)
plot(pcs)
M <- pcs$loadings[,1:2]
t(M)


plot(-1*pcs$scores[,1], -1*pcs$scores[,2],col=dat$label1)

library(fpc)
plotcluster(data[,-c(100,101)], dat$label1) # discriminant projection plot
plotcluster(data[,-c(100,101)], dat$label2) # discriminant projection plot


clusplot(data[,-c(100,101)], dat$label1, color=TRUE, shade=TRUE, 
         labels=2, lines=0)
##################################

plot(cmd2, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
text(cmd2, label = ifelse(dat$label1==1, "1", "2"),
     col=dat$label1)

plot(cmd1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
points(cmd1, pch = ifelse(dat$label1==1, 1, 3),
     col = ifelse(dat$study=="4703","red","blue"))
legend(0.23,0.315,c("data1","data2"),pch=15,text.width = 0.1,
       y.intersp=0.5, x.intersp = 1,
       col = c("red","blue"))

############################

plot(iso2, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
text(iso2, label = ifelse(dat$label1==1, "1", "2"),
     col=c(1:4))

plot(iso1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
text(iso1, label = ifelse(dat$label1==1, "1", "2"),
     col=dat$label1)



# 3 cluster
plot(cmd2, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
text(cmd2, label = ifelse(dat$label2==1, "1",
                          ifelse(dat$label2==2,"2","3")),
     col=dat$label2)

plot(cmd1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
text(cmd1, label = ifelse(dat$label2==1, "1",
                          ifelse(dat$label2==2,"2","3")),
     col=dat$label2)
################
plot(cmd1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
points(cmd1, pch = ifelse(dat$label2==1, 3, ifelse(dat$label2==2,1,2)),
     col = ifelse(dat$study=="4703","red","blue"))
legend(0.23,0.315,c("data1","data2"),pch=15,text.width = 0.1,
       y.intersp=0.5, x.intersp = 1,
       col = c("red","blue"))
 ###############
plot(iso2, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
text(iso2, label = ifelse(dat$label2==1, "1",
                          ifelse(dat$label2==2,"2","3")),
     col=dat$label2)

plot(iso1, type = "n", xlab = "Scaling Dimension 1", ylab = "Scaling Dimension 2")
text(iso1, label = ifelse(dat$label2==1, "1",
                          ifelse(dat$label2==2,"2","3")),
     col=dat$label2)




#### Cluster analysis
clinical_data <- data.frame(AGE = c(clinical_data1$AGE, clinical_data2$AGE))
clinical_data$label1 <- c(clinical_data1$label1, clinical_data2$label1)
clinical_data$label2 <- c(clinical_data1$label2, clinical_data2$label2)
clinical_data$label3 <- c(clinical_data1$label3, clinical_data2$label3)
clinical_data$BMI <- c(clinical_data1$basebmi, clinical_data2$basebmi)
clinical_data$FEV1 <- c(clinical_data1$FEV1pred, clinical_data2$FEV1pred)
clinical_data$PKyears <- c(clinical_data1$packyears, clinical_data2$PKyears)
clinical_data$race <- c(clinical_data1$race, clinical_data2$race)
clinical_data$ethnicity <- c(clinical_data1$ethnicity, clinical_data2$ethnicity)
clinical_data$sex <- c(clinical_data1$sex, clinical_data2$sex)
clinical_data$smkostat <- c(clinical_data1$ptsmokstat, clinical_data2$smokstat)
clinical_data$calcCOPD <- c(clinical_data1$calcCOPD, clinical_data2$calcCOPD)

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

######3-cluster################################
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
fisher.test(table(clinical_data$ethnicity, clinical_data$label2))
fisher.test(table(clinical_data$smkostat, clinical_data$label2))
fisher.test(table(clinical_data$calcCOPD, clinical_data$label2))

#####################################

# continuous data
mean1 <- function(x){
  mean(x, na.rm = T)
}
sd1 <- function(x){
  sd(x, na.rm = T)
}
col_idx <- c(5:10,24,25,57:59,75:82,85:90)


# 2 cluster
df1 <- data.frame(variables = names(clinical_data1)[col_idx], label1 = rep(0,length(col_idx)), 
                  label2 = rep(0,length(col_idx)))
for (i in 1:length(df1$variables)){
  df1$total[i] <- paste(round(mean1(clinical_data1[,names(clinical_data1)[col_idx][i]]), digits = 1), "(",
                        round(sd1(clinical_data1[,names(clinical_data1)[col_idx][i]]), digits = 1), ")",sep = "")
  
  df1$label1[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,mean1)[1], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,sd1)[1], digits = 1), ")",sep = "")
  df1$label2[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,mean1)[2], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,sd1)[2], digits = 1), ")", sep = "")
}

# 3 clusters
df2 <- data.frame(variables = names(clinical_data1)[col_idx], label1 = rep(0,length(col_idx)), 
                  label2 = rep(0,length(col_idx)), label3 = rep(0,length(col_idx)))
for (i in 1:length(df1$variables)){
  df2$total[i] <- paste(round(mean1(clinical_data1[,names(clinical_data1)[col_idx][i]]), digits = 1), "(",
                         round(sd1(clinical_data1[,names(clinical_data1)[col_idx][i]]), digits = 1), ")",sep = "")
  
  df2$label1[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,mean1)[1], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,sd1)[1], digits = 1), ")",sep = "")
  df2$label2[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,mean1)[2], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,sd1)[2], digits = 1), ")", sep = "")
  df2$label3[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,mean1)[3], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],clinical_data1$label1,sd1)[3], digits = 1),")", sep = "")
}




# 4 clusters
df3 <- data.frame(variables = names(clinical_data1)[col_idx], label1 = rep(0,length(col_idx)), 
                  label2 = rep(0,length(col_idx)), label3 = rep(0,length(col_idx)),
                  label4 = rep(0,length(col_idx)))
for (i in 1:length(df2$variables)){
  df3$total[i] <- paste(round(mean1(clinical_data1[,names(clinical_data1)[col_idx][i]]), digits = 1), "(",
                        round(sd1(clinical_data1[,names(clinical_data1)[col_idx][i]]), digits = 1), ")",sep = "")
  
  
  df3$label1[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],
                                  clinical_data1$label2,mean1)[1], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],
                                  clinical_data1$label2,sd1)[1], digits = 1), ")",sep = "")
  df3$label2[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],
                                  clinical_data1$label2,mean1)[2], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],
                                  clinical_data1$label2,sd1)[2], digits = 1), ")", sep = "")
  df3$label3[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],
                                  clinical_data1$label2,mean1)[3], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],
                                  clinical_data1$label2,sd1)[3], digits = 1),")", sep = "")
  df3$label4[i] <- paste(round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],
                                  clinical_data1$label2,mean1)[4], digits = 1), "(",
                         round(by(clinical_data1[,names(clinical_data1)[col_idx][i]],
                                  clinical_data1$label2,sd1)[4], digits = 1),")", sep = "")
}





