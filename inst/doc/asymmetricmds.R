## ----asymscal, echo=FALSE, fig.width = 10, fig.height = 8----------------
library(asymmetry)
# m <- matrix(c(3,1,2,2,4,2.2,1.5,3.3,2.2,5.5,3,4), 6, 2, byrow=TRUE)
# plot(m)
# sk <- m%*%t(matrix(c(m[,2],-1*m[,1]),6,2))
# e<-matrix(1,6,1)
# sk - e%*%t(colSums(sk)) - as.matrix(rowSums(sk))%*%t(e)

x<-matrix(c(-2,-1,0,1,2, -2,-1,0,1,1.75, 2.5,1,0,1.2,2.3, -2,1,0,-1,2),nrow=10,ncol=2)
par(mfrow = c(1,3),pty="s")
plot(x, xlim =c(-2.5,2.5),ylim =c(-2.5,2.5),xlab='Dimension 1',ylab='Dimension 2', main = "Group Space")
points(x=-2,y=2.5,col = "red", pch = 19)
points(x=2,y=2.3,col = "blue", pch = 19)

plot(x%*%diag(c(.25,1)), xlim =c(-2.5,2.5),ylim =c(-2.5,2.5),xlab='Dimension 1',ylab='Dimension 2',main="Object 1 (.25,1)")
points(x=-2*.25,y=2.5,col = "red", pch = 19)
points(x=2*.25,y=2.3,col = "blue", pch = 19)

plot(x%*%diag(c(1,.25)), xlim =c(-2.5,2.5),ylim =c(-2.5,2.5),xlab='Dimension 1',ylab='Dimension 2', main = "Object 2 (1,.25)")
points(x=-2,y=2.5*.25,col = "red", pch = 19)
points(x=2,y=2.3*.25,col = "blue", pch = 19)



## ----slide-explain, echo=FALSE, fig.width = 10, fig.height = 8-----------
par(pty="s",mfrow=c(1,2))
rows <- 6
m <- matrix(c(2,0,-1.5,3,4,5,4,3,-1,-2,4,-1.5),nrow=6,2)
plot(m, xlab = "dimension 1", ylab = "dimension 2", xlim = c(-2,7), ylim = c(-2,7))
slope <- -.45
intercept <- 0
negslope <- -1/slope
abline(intercept,slope, col ="pink",lty=2)
arrows(0,0,2,intercept+slope*2,col ="red")
#abline(4-negslope*2,negslope)
neginter <- m[6,2]-negslope*m[6,1]
#abline(m[6,2]-negslope*m[6,1],negslope)
s1 <- (intercept-neginter)/(negslope-slope)
s2 <- intercept + s1*slope
segments(s1,s2,m[6,1],m[6,2])
#segments(4,3,2,4)
for (i in 1:rows)
{
neginter <- m[i,2]-negslope*m[i,1]
s1 <- (intercept-neginter)/(negslope-slope)
s2 <- intercept + s1*slope
segments(s1,s2,m[i,1],m[i,2],col = "blue")
	
}
labels <- c("A","B","C","D","E","F")
for (i in 1:rows)
{

text(-.15+m[i,1],.3+m[i,2],labels[i])
}
#
# second application
#
plot(m, xlab = "dimension 1", ylab = "dimension 2", xlim = c(-2,7), ylim = c(-2,7))
slope <- .75
intercept <- 0
negslope <- -1/slope
abline(intercept,slope, col ="pink",lty=2)
arrows(0,0,2,intercept+slope*2,col ="red")
#abline(4-negslope*2,negslope)
neginter <- m[6,2]-negslope*m[6,1]
#abline(m[6,2]-negslope*m[6,1],negslope)
s1 <- (intercept-neginter)/(negslope-slope)
s2 <- intercept + s1*slope
segments(s1,s2,m[6,1],m[6,2])
#segments(4,3,2,4)
for (i in 1:rows)
{
neginter <- m[i,2]-negslope*m[i,1]
s1 <- (intercept-neginter)/(negslope-slope)
s2 <- intercept + s1*slope
segments(s1,s2,m[i,1],m[i,2],col = "blue")
	
}
labels <- c("A","B","C","D","E","F")
for (i in 1:rows)
{

text(.05+m[i,1],.3+m[i,2],labels[i])
}

## ----englishtowns, fig.width = 8, fig.height = 8-------------------------
data(Englishtowns)
v<-slidevector(Englishtowns, ndim = 2, itmax = 250, eps = .001)
plot(v,col="blue",ylim=c(-150,300),xlim=c(-150,300))

