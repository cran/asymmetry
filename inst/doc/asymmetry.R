## ----echo = FALSE, message = FALSE, warning=FALSE-----------------------------
library(knitr)
knitr::opts_chunk$set(
   cache = FALSE,
   dpi = 75,
   fig.width = 6, fig.height = 6,
   fig.keep = 'high',
  # comment = "#>",
  tidy = FALSE)


## ----dataexample--------------------------------------------------------------
library("asymmetry")
data("studentmigration")
idx <- c(3,4,25,27,31) #select five countries
studentmigration[idx,idx]

## ----example------------------------------------------------------------------
q1 <- skewsymmetry(studentmigration[idx,idx])
q1$A

## ----symm---------------------------------------------------------------------
q1$S

## ----example2-----------------------------------------------------------------
summary(q1)

## ----clustering---------------------------------------------------------------
clus <- hclust(as.dist(1/q1$S))
plot(clus,xlab=NA,sub=NA)

## ----linear-------------------------------------------------------------------
   q1$linear

## ----example3, echo=TRUE------------------------------------------------------
library(RColorBrewer)
# creates a color palette from red to blue
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
col_breaks = c(seq(-4000,-.001,length=100),  # negative values are red
  seq(-.001,0.01,length=100),                # zeroes are white
  seq(0.01,4000,length=100))                 # positive values are blue

hmap(q1, col = my_palette)

## ----example4-----------------------------------------------------------------
data(studentmigration)
idx <- c(18,22,27,2,13,31) #select 6 countries
q1 <- skewsymmetry(studentmigration[idx,idx])
q1$A


## ----studentmigration1 , fig.width = 8, fig.height = 8------------------------
# creates a color palette from red to blue
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 299)
col_breaks = c(seq(-4000,-.001,length=100),  # negative values are red
  seq(-.001,0.01,length=100),                # zeroes are white
  seq(0.01,4000,length=100))                 # positive values are blue
data(studentmigration)
hmap(studentmigration[idx,idx], dominance = FALSE, col = my_palette, key = FALSE, xlab = "Destination country", ylab = "Home country", colsep = c(1:6), rowsep = c(1:6))


## ----englishtowns, fig.width = 8, fig.height = 8------------------------------
data(Englishtowns)
v<-slidevector(Englishtowns, ndim = 2, itmax = 2500, eps = .0000001, verbose = FALSE)
plot(v,col="blue",ylim=c(-300,300),xlim=c(-300,300))


## ----englishtownsdecomp-------------------------------------------------------

q2 <- skewsymmetry(v$resid)
summary(q2)


## ----migrationunique----------------------------------------------------------

data("studentmigration")
mm<-studentmigration
mm[mm==0]<-.5          # replace zeroes by a small number
mm <- -log(mm/sum(mm)) # convert similarities to dissimilarities
v<-mdsunique(mm, ndim = 2, itmax = 2100, verbose=FALSE, eps = .0000000001)
plot(v, yplus = .3, ylim = c(-4.5, 4), xlim = c(-4.5, 4))



## ----slide-explain, echo=FALSE, fig.width = 10, fig.height = 8----------------
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

