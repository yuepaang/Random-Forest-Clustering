# Adjusted Rand index between two partitions
Rand <- function(tab,adjust=T) {
  
  ##########################################################################
  # The function computes the (adjusted) Rand index between two partitions #
  # Copyright Steve Horvath and Luohua Jiang, UCLA, 2003                   #
  ##########################################################################
  
  # helper function
  choosenew <- function(n,k) {
    n <- c(n); out1 <- rep(0,length(n));
    for (i in c(1:length(n)) ){
      if ( n[i]<k ) {out1[i] <- 0}
      else {out1[i] <- choose(n[i],k) }
    }
    out1
  }
  
  a <- 0; b <- 0; c <- 0; d <- 0; nn <- 0
  n <- nrow(tab)
  for (i in 1:n) {
    for(j in 1:n) {
      a <- a+choosenew(tab[i,j],2)
      nj <- sum(tab[,j])
      c <- c+choosenew(nj,2)
    }
    ni <- sum(tab[i,])
    b <- b+choosenew(ni,2)
    nn <- nn+ni
  }
  if(adjust==T) {
    d <- choosenew(nn,2)
    adrand <- (a-(b*c/n)/d)/(0.5*(b+c/n)-(b*c/n)/d)
    adrand
  } else {
    b <- b-a
    c <- c/n-a
    d <- choosenew(nn,2)-a-b-c
    rand <- (a+d)/(a+b+c+d)
    rand
  }
}


