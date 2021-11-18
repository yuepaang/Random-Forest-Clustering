# imaging data analysis
airway_data_1 <- read.csv(".csv")
airway_data_2 <- read.csv(".csv")
body_composition_data_1 <- read.csv(".csv")
body_composition_data_2 <- read.csv(".csv")
densitometry_data_1 <- read.csv("csv")
densitometry_data_2 <- read.csv("csv")
phenotype_data_1 <- read.csv("csv")
phenotype_data_2 <- read.csv("csv")

# combine the data according to "cn"
idx1 <- intersect(intersect(intersect(airway_data_1$cn, body_composition_data_1$cn), 
                    densitometry_data_1$cn), phenotype_data_1$cn)
idx2 <- intersect(intersect(intersect(airway_data_2$cn, body_composition_data_2$cn), 
                            densitometry_data_2$cn), phenotype_data_2$cn)

data1 <- cbind(airway_data_1[which(airway_data_1$cn %in% idx1),],
                 body_composition_data_1[which(body_composition_data_1$cn %in% idx1),],
                 densitometry_data_1[which(densitometry_data_1$cn %in% idx1),],
                 phenotype_data_1[which(phenotype_data_1$cn %in% idx1),])

data2 <- cbind(airway_data_2[which(airway_data_2$cn %in% idx2),],
                 body_composition_data_2[which(body_composition_data_2$cn %in% idx2),],
                 densitometry_data_2[which(densitometry_data_2$cn %in% idx2),],
                 phenotype_data_2[which(phenotype_data_2$cn %in% idx2),])
# extract the features to analyze
imaging_data <- rbind(data1, data2)[,-c(1:6,79:86, 88:92)]

# remove the missing value 
na_count <- function(x){
  sum(is.na(x))
}

col_na <- apply(imaging_data, 2,na_count)

row_na <- apply(imaging_data[,names(col_na[which(col_na < 50)])], 1, na_count)
dat <- imaging_data[,names(col_na[which(col_na < 50)])][which(row_na < 30),]
dat <- dat[complete.cases(dat),]

# clustering
source("RFDist.R")
source("PAM.R")
source("AdjRand.R")


no.forests <- 100
no.trees <- 4000
distRF <- RFdist(dat[,-c(27,28)], mtry1=20, no.trees, no.forests, addcl1=T,addcl2=T,imp=T, oob.prox1=T)
