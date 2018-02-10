knn <- function(klist,x.train,y.train,x.test) {
  # k-nearest neighbors classification
  # 
  # klist is a list of values of k to be tested
  # x.train, y.train: the training set
  # x.test: the test set
  # Output: a matrix of predictions for the test set (one column for each k in klist)	
  # Number of training and test examples
  n.train <- nrow(x.train)
  n.test <- nrow(x.test)
  
  # Matrix to store predictions
  p.test <- matrix(NA, n.test, length(klist))
  
  # Vector to store the distances of a point to the training points
  dsq <- numeric(n.train)
  
  # Loop on the test instances
  for (tst in 1:n.test)
  {
    # Compute distances to training instances
    for (trn in 1:n.train)
    {
      dsq[trn] <- sum((x.train[trn,] - x.test[tst,])^2)
    }
    
    # Sort distances from smallest to largest
    ord <- order(dsq)
    
    # Make prediction by averaging the k nearest neighbors
    for (ik in 1:length(klist)) {
      p.test[tst,ik] <- mean(y.train[ord[1:klist[ik]]])
    }
  }
  
  # Return the matrix of predictions
  invisible(p.test)
}

knn.cv <- function(klist,x.train,y.train,nfolds) {
  # Cross-validation for kNN
  #
  # Perform nfolds-cross validation of kNN, for the values of k in klist
  
  # Number of instances
  n.train <- nrow(x.train)
  
  # Matrix to store predictions
  p.cv <- matrix(NA, n.train, length(klist))
  
  # Prepare the folds
  s <- split(sample(n.train),rep(1:nfolds,length=n.train))
  
  # Cross-validation
  for (i in seq(nfolds)) {
    p.cv[s[[i]],] <- knn(klist,x.train[-s[[i]],], y.train[-s[[i]]], x.train[s[[i]],])
  }
  
  # Return matrix of CV predictions
  invisible(p.cv)
}




# Download prostate data
con = url ("https://web.stanford.edu/~hastie/ElemStatLearn/datasets/prostate.data")
prost=read.csv(con,row.names=1,sep="\t")
# Alternatively, load the file and read from local file as follows
# prost=read.csv('prostate.data.txt',row.names=1,sep="\t")

# Create training and test sets
x.train <- prost[prost$train,1:8]
x.test <- prost[!prost$train,1:8]
y.train <- prost[prost$train,9]
y.test <- prost[!prost$train,9]

# Make predictions by kNN
klist <- seq(20) # we test all values of k
nfolds <- 5 # we make 5-fold cross-validation
y.pred.train <- knn(klist,x.train,y.train,x.train)
y.pred.test <- knn(klist,x.train,y.train,x.test)
y.pred.cv <- knn.cv(klist,x.train,y.train,nfolds)

# Compute mean-square error (MSE)
mse.train <- apply((y.pred.train - y.train)^2 , 2, mean)
mse.test <- apply((y.pred.test - y.test)^2 , 2, mean)
mse.cv <- apply((y.pred.cv - y.train)^2 , 2, mean)

# Plot MSE as a function of k
plot(mse.train , ylim=c(0,2) , type='l' , xlab='k' , ylab='MSE', col=1 , lwd=2)
lines(mse.test , col=2 , lwd=2)
lines(mse.cv, col=3 , lwd=2)
legend("bottomright",legend=c('Train','Test','CV'),text.col=seq(3) , lty=1 , col=seq(3))
