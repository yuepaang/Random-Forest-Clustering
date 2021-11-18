data1 <- read.csv("data1.csv", na.strings = c("#N/A"))
data2 <- read.csv("data2.csv",na.strings = c("#N/A",'.'))

#separate into different label groups
library(dplyr)

comptable <- function(data){
  
  variable <- names(data)[c(2:16)]
  
  
  #categorical or continuous
  comp <- function(name){
    
    x <- data[,name]
    m <- length(unique(x))
    if (m > 6) {
      mn <- round(mean(x,na.rm = T),digits = 2)
      sderr <- round(sd(x,na.rm = T),digits = 2)
      result <- paste0(mn," ± ", sderr)
      result <- data.frame(variable = name,stat = result)
    }else{
      title <- data.frame(variable = name,stat = "")
      cat <- as.character(unique(data2[,name]))
      for (i in cat) {
        n <- length(x)
        count <- sum(x == i,na.rm = T)
        percentage <- round(count/n*100, digits = 2)
        result <- paste0(count,"(", percentage,")")
        result <- data.frame(variable = i,stat = result)
        title <- rbind(title,result)
      }
      result <- title
    }
    return(result)
  }
  
  stat_table <- c()
  for (name in variable) {
    stat_table <- rbind(stat_table,comp(name))
  }
  return(stat_table)
}

#separate into groups
mean_1 <- data2 %>% filter(X2_means_cluster == 1)
mean_2 <- data2 %>% filter(X2_means_cluster == 2)
rf_1 <- data2 %>% filter(X2_RF_cluster == 1)
rf_2 <- data2 %>% filter(X2_RF_cluster == 2)
hc_1 <- data2 %>% filter(X2_HC_cluster == 1)
hc_2 <- data2 %>% filter(X2_HC_cluster == 2)

sc_1 <- data2 %>% filter(sequencing_cluster == 1)
sc_2 <- data2 %>% filter(sequencing_cluster == 2)

stat_table_whole <- cbind(comptable(mean_1),comptable(mean_2)[,2],comptable(rf_1)[,2],comptable(rf_2)[,2],
                          comptable(hc_1)[,2],comptable(hc_2)[,2], comptable(sc_1)[,2], comptable(sc_2)[,2])
names(stat_table_whole) <- c("variable","mean_1","mean_2","rf_1","rf_2","hc_1",
                             "hc_2", "sc_1", "sc_2")
write.csv(stat_table_whole,file="new.csv")
#separate into groups
mean_1 <- data2 %>% filter(X3_means_cluster == 1)
mean_2 <- data2 %>% filter(X3_means_cluster == 2)
mean_3 <- data2 %>% filter(X3_means_cluster == 3)
rf_1 <- data2 %>% filter(X3_RF_cluster == 1)
rf_2 <- data2 %>% filter(X3_RF_cluster == 2)
rf_3 <- data2 %>% filter(X3_RF_cluster == 3)

stat_table_whole2 <- cbind(comptable(mean_1),comptable(mean_2)[,2],comptable(mean_3)[,2],
                           comptable(rf_1)[,2],comptable(rf_2)[,2],comptable(rf_3)[,2])
names(stat_table_whole2) <- c("variable","mean_1","mean_2","mean_3","rf_1","rf_2","rf_3")

################################################data1
#separate into groups
mean_1 <- data1 %>% filter(X2_means_cluster == 1)
mean_2 <- data1 %>% filter(X2_means_cluster == 2)
rf_1 <- data1 %>% filter(X2_RF_cluster == 1)
rf_2 <- data1 %>% filter(X2_RF_cluster == 2)

stat_table_whole3 <- cbind(comptable(mean_1),comptable(mean_2)[,2],comptable(rf_1)[,2],comptable(rf_2)[,2])
names(stat_table_whole3) <- c("variable","mean_1","mean_2","rf_1","rf_2")

#separate into groups
mean_1 <- data1 %>% filter(X3_means_cluster == 1)
mean_2 <- data1 %>% filter(X3_means_cluster == 2)
mean_3 <- data1 %>% filter(X3_means_cluster == 3)
rf_1 <- data1 %>% filter(X3_RF_cluster == 1)
rf_2 <- data1 %>% filter(X3_RF_cluster == 2)
rf_3 <- data1 %>% filter(X3_RF_cluster == 3)

stat_table_whole4 <- cbind(comptable(mean_1),comptable(mean_2)[,2],comptable(mean_3)[,2],comptable(rf_1)[,2],comptable(rf_2)[,2],comptable(rf_3)[,2])
names(stat_table_whole4) <- c("variable","mean_1","mean_2","mean_3","rf_1","rf_2","rf_3")

write.csv(stat_table_whole, "1.csv")
write.csv(stat_table_whole2, "2.csv")
write.csv(stat_table_whole3, "3.csv")
write.csv(stat_table_whole4, "4.csv")

library(irr)
kappa2(cbind(data2$X2_HC_cluster, data$X2_means_cluster),'unweighted')
kappa2(cbind(data2$X2_RF_cluster, data$X2_means_cluster),'unweighted')
kappa2(cbind(data2$X2_HC_cluster, data$X2_RF_cluster),'unweighted')
kappam.fleiss(cbind(data$X2_HC_cluster, data$X2_RF_cluster,data$X2_means_cluster), exact = T)

idx <- which(data$sequencing_cluster=="1" | data$sequencing_cluster=="2")
kappa2(cbind(data$X2_HC_cluster[idx],
             as.numeric(as.character(data$sequencing_cluster[idx]))),'unweighted')
kappa2(cbind(data$X2_RF_cluster[idx],
             as.numeric(as.character(data$sequencing_cluster[idx]))),'unweighted')
kappa2(cbind(data$X2_means_cluster[idx],
             as.numeric(as.character(data$sequencing_cluster[idx]))),'unweighted')


# 3 KM vs 2 HC
kappa2(cbind(data$X2_HC_cluster, data$X3_means_cluster),'unweighted')
chisq.test(as.matrix(table(data$X2_HC_cluster,data$X3_means_cluster)))
kappa2(cbind(data$X2_HC_cluster, data$X3_RF_cluster),'unweighted')
chisq.test(as.matrix(table(data$X2_HC_cluster,data$X3_RF_cluster)))

chisq.test(as.matrix(table(data$X2_means_cluster[idx],
                           as.numeric(as.character(data$sequencing_cluster[idx])))))
chisq.test(as.matrix(table(data$X2_RF_cluster[idx],
                           as.numeric(as.character(data$sequencing_cluster[idx])))))
chisq.test(as.matrix(table(data$X2_HC_cluster[idx],
                           as.numeric(as.character(data$sequencing_cluster[idx])))))


table(data$X2_HC_cluster,data$X2_means_cluster)
chisq.test(as.matrix(table(data$X2_HC_cluster,data$X2_means_cluster)))
chisq.test(as.matrix(table(data$X2_HC_cluster,data$X2_RF_cluster)))
chisq.test(as.matrix(table(data$X2_RF_cluster,data$X2_means_cluster)))
# 3 vs 3

kappa2(cbind(data$X3_RF_cluster, data$X3_means_cluster),'unweighted')
chisq.test(as.matrix(table(data$X3_RF_cluster,data$X3_means_cluster)))

