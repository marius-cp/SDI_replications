#####################################################################
#
# Various helper functions for date handling etc
#
#####################################################################

# Like "read.table", but with changed defaults
read.table2 <- function(s) read.table(s,sep=",",header=TRUE,stringsAsFactors=FALSE)  

# Functions for dealing with character dates
tdiff <- function(d1, d2, freq){
  y1 <- substr(d1,1,4)
  y2 <- substr(d2,1,4)
  if (freq == 4) aux <- 6 else aux <- 7
  m1 <- substr(d1,6,aux) 
  m2 <- substr(d2,6,aux)
  return((as.numeric(y2)-as.numeric(y1))*freq+as.numeric(m2)-as.numeric(m1))
}
plust <- function(d1, t, freq){
  if (freq == 4){
   aux1 <- 6
  } else if (freq == 12) {
   aux1 <- 7
  }
  y1 <- as.numeric(substr(d1,1,4))
  p1 <- as.numeric(substr(d1,6,aux1))
  aux <- y1 + (p1-1)/freq + t/freq
  y2 <- floor(aux)
  p2 <- (aux-y2)*freq+1
  if (freq == 4){
	return(paste0(y2,"Q",p2))
  } else if (freq == 12){
	return(paste0(y2,"M", gsub(" ", "0", format(p2, width = 2))))
  }
}
plusq <- function(d1, t) plust(d1, t, 4)
plusm <- function(d1, t) plust(d1, t, 12)
qdiff <- function(d1, d2) tdiff(d1, d2, 4)
mdiff <- function(d1, d2) tdiff(d1, d2, 12)

# Function to make numerical axis in time series plots
# x is a vector of character strings, in the form "1998Q4" for the fourth quarter of 1998
stof <- function(x){
  as.numeric(substr(x, 1, 4)) + (as.numeric(substr(x, 6, 6))-1)/4
}