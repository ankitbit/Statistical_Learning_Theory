# Prostate data for illustrating ridge regression 
# 
# Data downloaded from
# http://statweb.stanford.edu/~tibs/ElemStatLearn/
# 10-04-2016
#
prostate <- read.table("prostate_data.txt", header=TRUE, row.names = 1)
plot(prostate)

train.sample <- which(prostate$train==TRUE)

#use.only <- 1:dim(prostate)[1]
use.only <- train.sample

Y <- scale( prostate$lpsa[use.only], center=TRUE, scale=FALSE)
X <- scale( as.matrix(prostate[use.only,1:8]), center=TRUE, scale=TRUE)
n <- dim(X)[1]
p <- dim(X)[2]

XtX <- t(X)%*%X 
d2 <- eigen(XtX,symmetric = TRUE, only.values = TRUE)$values


lambda.max <- 1e5
n.lambdas <- 25
lambda.v <- exp(seq(0,log(lambda.max+1),length=n.lambdas))-1
#lambda.v <- seq(0,lambda.max,length=n.lambdas)

#############
# estimated coefficients path
#############
beta.path <- matrix(0,nrow=n.lambdas, ncol=p)
diag.H.lambda <- matrix(0,nrow=n.lambdas, ncol=n)
for (l in 1:n.lambdas){ 
  lambda <- lambda.v[l]
  H.lambda.aux <- t(solve(XtX + lambda*diag(1,p))) %*% t(X) 
  beta.path[l,] <-  H.lambda.aux %*% Y
  H.lambda <- X %*% H.lambda.aux 
  diag.H.lambda[l,] <- diag(H.lambda)
} 
plot(c(-1,log(1+lambda.v[n.lambdas])), range(beta.path),type="n",
     xlab="log(1+lambda)",ylab="coefficients")
abline(h=0,lty=2)
for(j in 1:p){
  lines(log(1+lambda.v),beta.path[,j],col=4)
  points(log(1+lambda.v),beta.path[,j],pch=19,cex=.7,col=4)
}
text(0*(1:p), beta.path[1,],names(prostate)[1:p],pos=2)

#############
# effective degrees of freedom
#############
df.v <- numeric(n.lambdas)
for (l in 1:n.lambdas){
  lambda <- lambda.v[l]
  df.v[l] <- sum(d2/(d2+lambda)) 
}
plot(log(1+lambda.v),df.v)
points(0*df.v,df.v,col=2,pch=19)
text(0*df.v,df.v,round(df.v,2),col=2,pch=19,pos=4)

# linear interpolation to obtain the values lambda.vv
# such that the corresponding df are 0,1,...,8 (approx.)
lambda.vv <- approx(x=df.v,y=lambda.v,xout=0:8)$y
lambda.vv[1] <- lambda.v[n.lambdas]
df.vv <- numeric(length(lambda.vv))
for (l in 1:length(lambda.vv)){
  lambda <- lambda.vv[l]
  df.vv[l] <- sum(d2/(d2+lambda)) 
}
print(df.vv)

# another way to compute df's
trace.H.lambda <- apply(diag.H.lambda,1,sum)
print(trace.H.lambda - df.v)

# estimated coefficients path against effective degrees of freedom
plot(c(0,p+1), range(beta.path),type="n",xlab="df(lambda)",ylab="coefficients")
abline(h=0,lty=2)
for(j in 1:p){
  lines(df.v,beta.path[,j],col=4)
  points(df.v,beta.path[,j],pch=19,cex=.7,col=4)
}
text(p+0*(1:p), beta.path[1,],names(prostate)[1:p],pos=4)


#############
# choosing lambda by leave-one-out cross validation
#############
PMSE.CV <- numeric(n.lambdas)
for (l in 1:n.lambdas){
  lambda <- lambda.v[l]
  PMSE.CV[l] <- 0
  for (i in 1:n){
#   m.Y.i <- mean(Y[-i])
    m.Y.i <- 0
    X.i <- X[-i,]; Y.i <- Y[-i]-m.Y.i
    Xi <- X[i,]; Yi <- Y[i]
    beta.i <- solve(t(X.i)%*%X.i + lambda*diag(1,p)) %*% t(X.i) %*% Y.i
    hat.Yi <- Xi %*% beta.i + m.Y.i
    PMSE.CV[l] <- PMSE.CV[l] + (hat.Yi-Yi)^2
  }
  PMSE.CV[l] <- PMSE.CV[l]/n
}
lambda.CV <- lambda.v[which.min(PMSE.CV)]
df.CV <- df.v[which.min(PMSE.CV)]

plot(log(1+lambda.v), PMSE.CV)
abline(v=log(1+lambda.CV),col=2,lty=2)

plot(df.v, PMSE.CV)
abline(v=df.CV,col=2,lty=2)


#############
#### computing PMSE.CV using the diagonal of H.lambda matrices 
#############
PMSE.CV.H.lambda <- numeric(n.lambdas)
for (l in 1:n.lambdas){
  lambda <- lambda.v[l]
  hat.Y <- X %*% beta.path[l,]
  PMSE.CV.H.lambda[l] <- sum( ((Y-hat.Y)/(1-diag.H.lambda[l,]))^2 )/n
}
lambda.CV.H.lambda <- lambda.v[which.min(PMSE.CV.H.lambda)]
df.CV.H.lambda <- df.v[which.min(PMSE.CV.H.lambda)]

plot(df.v, PMSE.CV.H.lambda)
points(df.v, PMSE.CV,col=3,pch=19,cex=.5)
abline(v=df.CV.H.lambda,col=2,lty=2)



#############
#### computing PMSE.GCV in ridge regression
#############
PMSE.GCV <- numeric(n.lambdas)
for (l in 1:n.lambdas){
  lambda <- lambda.v[l]
  hat.Y <- X %*% beta.path[l,]
  nu <- sum(diag.H.lambda[l,])
  PMSE.GCV[l] <- sum( ((Y-hat.Y)/(1-nu/n))^2 )/n
}
lambda.GCV <- lambda.v[which.min(PMSE.GCV)]
df.GCV <- df.v[which.min(PMSE.GCV)]

plot(df.v, PMSE.GCV)
points(df.v, PMSE.CV,col=6,pch=19,cex=.75)
abline(v=df.GCV,col=1,lty=2,lwd=3)
abline(v=df.CV.H.lambda,col=6,lty=6)
legend("top",c("PMSE.GCV","PMSE.CV","lambda.GCV","lambda.CV"),
       pch=c(1,19,NA,NA),lty=c(0,0,2,6),lwd=c(0,0,3,1),col=c(1,6,1,6))
