# two-way hierarchy clustering
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering visualization
library(dendextend) # for comparing two dendrograms
library(ggplot2)

dat <- read.csv("data.csv",stringsAsFactors=FALSE)
train <- dat[,3:100]

# make sure all are numeric data
for (i in 1:98){
  if (!is.numeric(train[,i])){
    train[,i] <- as.numeric(train[,i])
  }
}
# two NA line 4703-116, 4703-172
train <- na.omit(train)
# Dissimilarity matrix
d <- dist(train, method = "euclidean")
hc <- hclust(d, method = "ward.D")
plot(hc, cex = 0.6, hang=-1)

# Cut tree into 2 groups
sub_grp <- cutree(hc, k = 2)

# Number of members in each cluster
table(sub_grp)
tmp <- cmdscale(d, k=2)
x <- tmp[,1]
y <- tmp[,2]
p <- ggplot(data.frame(x,y),aes(x,y))
p + geom_point(size=3, alpha=0.8,aes(color=factor(sub_grp)))
# get the agglomerative coefficient
# which measures the amount of clustering structure found (values closer to 1 suggest strong clustering structure).
hc2 <- agnes(train, method = "complete")
hc2$ac


x  <- t(as.matrix(scale(train)))
dd.row <- as.dendrogram(hclust(dist(x)))
row.ord <- order.dendrogram(dd.row)

dd.col <- as.dendrogram(hclust(dist(t(x))))
col.ord <- order.dendrogram(dd.col)


myTheme <- modifyList(custom.theme(region = brewer.pal(11, "RdBu")),
                      list(
                        strip.background=list(col='gray'),
                        panel.background=list(col='gray')
                      ))
library(lattice)
library(latticeExtra)
levelplot(t(x[row.ord, col.ord]),
          par.settings=myTheme,
          aspect = "fill",
          xlab="cases",
          ylab="features",
          scales = list(x = list(rot = 90)),
          colorkey = list(space = "left"),
          labels = list(cex=0.02),
          at=pretty(c(-20,25), n=8),
          legend =
            list(right =
                   list(fun = dendrogramGrob,
                        args =
                          list(x = dd.row, ord = row.ord,
                               side = "right",
                               size = 10
                               )),
                 top =
                   list(fun = dendrogramGrob,
                        args =
                          list(x = dd.col, 
                               side = "top",
                               size = 10
                               ))))

rownames(train)[col.ord[1:144]]
rownames(train)[col.ord[145:414]]
