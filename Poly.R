library(dplyr)
library(polynom)

##########################################################
## Function to find the characteristic polynom of a matrix
##########################################################

make.character.poly <- function(A){
  n <- length(A[, 1])
  C <- numeric(n+1)
  E <- matrix(ncol=n, nrow=n, 0)
  liambda <- eigen(A)$values
  C[1] <- 1
  
  for(i in 1:(n)){
    C[2:(i+1)] <- C[2:(i+1)] - liambda[i]*C[1:i]
  }
  
  C <- rev(C)
  
  char.poly <- C %>% polynomial() 
  char.poly %>% plot(lwd=2, col="blue", xlab="liambda", ylab="Values of polynom")
  abline(h=0, col="red", lwd=2)
  
  return(char.poly)
}

############################################################

## example nr1
K <- matrix(nrow=3, ncol=3)
K[, 1] <- c(1, 4, 3)
K[, 2] <- c(2, 2, 2)
K[, 3] <- c(3, 2, 3)

make.character.poly(K)

## example nr2

A1 <- matrix(ncol=2, nrow=2, c(2, -1, 1, 0))
make.character.poly(A1)


###############################################
### Inverse iteration method implementation ###
###############################################
Tm <- matrix(ncol=4, nrow=4)
Tm[, 1] <- c(2.34, 2, 0, 0)
Tm[, 2] <- c(2, 2.34, 2, 0)
Tm[, 3] <- c(0, 2, 2.34, 2)
Tm[, 4] <- c(0, 0, 2, 2.34)

C <- matrix(0, ncol=4, nrow=4)
diag(C) <- 0.3
E <- matrix(0, ncol=4, nrow=4)
diag(E) <- 1
A <- Tm + 3*C # We will try to find this matrix' eigen value and eigen vector

eps <- 0.0001 #The measure of error.

make.character.poly(A)
eigen(A)
liambda.prad <- (4+5.5)/2 ## We find this value analyzing the previous plot
counter <- 0
i       <- 1
x.prad  <- c(1, 0, 0, 1)

while(counter!=1){
  
  y.prad       <- solve(A - E * as.numeric(liambda.prad), x.prad) 
  x.new        <- y.prad/((y.prad)^2 %>% sum() %>% sqrt())
  liambda.new  <- (t((A %*% x.new)) %*% x.new)/ (x.new %*% x.new)
  
  dist.x       <- ((-1^i)*x.new - x.prad)^2 %>% sum() %>% sqrt()
  dist.liambda <- (liambda.new - liambda.prad)^2 %>% sum() %>% sqrt()
  
  if(dist.x <= eps & dist.liambda <= eps) counter=1
  
  cat("--------------- \n")
  cat(i, dist.x, "\n")         # We want the distances to converge
  cat(i, dist.liambda, "\n")   # We want the distances to converge
  cat("--------------- \n")
  x.prad       <- x.new
  liambda.prad <- liambda.new
  i <- i+1
}

x.prad           # eigen vector
liambda.prad     # eigen value
