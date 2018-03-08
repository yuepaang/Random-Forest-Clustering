# PAM based on silhouette width
pamNew <- function (x, k, diss1 = inherits(x, "dist"), metric1 = "euclidean")
{
  
  #############################################################################################################
  # A new pam clustering function which corrects the clustering membership based on the sillhouette strength. #
  # The clustering membership of an observation with a negative sillhouette strength is reassigned to its     #
  # neighboring cluster.                                                                                      #
  # The inputs of the function are similar to the original 'pam' function.                                    #
  # The function returns a vector of clustering labels.                                                       #
  # Copyright 2003 Tao Shi and Steve Horvath (last modified 10/31/03)                                         #
  #############################################################################################################
  
  if (diss1)
  {
    if (!is.null(attr(x, "Labels"))) { original.row.names <- attr(x, "Labels")}
    names(x) <- as.character(c(1:attr(x, "Size")))
  } 
  else
  {
    if(!is.null(dimnames(x)[[1]])) { original.row.names <- dimnames(x)[[1]]}
    row.names(x) <- as.character(c(1:dim(x)[[1]]))
  }
  pam1 <- pam(x,k,diss=diss1, metric=metric1)
  label2 <- pam1$clustering
  silinfo1 <- pam1$silinfo$widths
  index1 <- as.numeric(as.character(row.names(silinfo1)))
  silinfo2 <- silinfo1[order(index1),]
  labelnew <- ifelse(silinfo2[,3]<0, silinfo2[,2], silinfo2[,1])
  names(labelnew) <- original.row.names
  labelnew    
}