#' @title Bhattacharyya Coefficient of two distributions
#'
#' @description Calculates the Bhattacharyya Coefficient (BC) between two distributions. BC is the sum of the square root of the multiplication of the elements in bin i for both distributions. It ranges from 0 (no overlap) to 1 (full overlap). v.0.1. Used within \code{eval_brt}. Included to ensure functionality.
#'
#' @source From Camrin Brawn (WHOI): https://zenodo.org/records/7971532.
#'
#' @param x a numeric vector of length >= 2.
#' @param y a numeric vector of length >= 2.
#' @param bw can be either a fixed value of bins or a function to calculate the bandwidth of each bin (see ?bw.nrd). Default is bw.nrd0.
#' @param ...  any optional arguments to be passed to the given bw function.
##########################c
#----
#guillert(at)tcd.ie - 28/11/2014
##########################
#Requirements:
#-R 3
##########################

bhatt_coeff_brt <-function(x,y, bw=bw.nrd0, ...) {
  #SANITIZING
  #x
  if(!is.numeric(x)) {
    stop("'x' must be numeric.")
  }
  if(length(x) < 2) {
    stop("'x' need at least two data points.")
  }

  #y
  if(!is.numeric(y)) {
    stop("'y' must be numeric.")
  }
  if(length(y) < 2) {
    stop("'y' need at least two data points.")
  }

  #bw
  if(length(bw) != 1) {
    stop("'bw' must be either a single numeric value or a single function.")
  }
  if(!is.function(bw)) {
    if(!is.numeric(bw)) {
      stop("'bw' must be either a single numeric value or a single function.")
    }
  }
  #Avoiding non-entire numbers
  if(is.numeric(bw)) {
    bw<-round(bw)
  }

  #BHATTACHARYYA COEFFICIENT
  #sum(sqrt(x relative counts in bin_i * y relative counts in bin_i))

  #Setting the right number of bins (i)
  if(is.function(bw)) {
    #Bin width
    band.width<-bw(c(x,y), ...)
    #Bin breaks
    bin.breaks<-seq(from=min(c(x,y)), to=max(c(x,y)+band.width), by=band.width) #adding an extra bandwith to the max to be sure to include all the data
    #Number of bins
    bin.n<-length(bin.breaks)-1
  } else {
    #Bin breaks
    bin.breaks<-hist(c(x,y), breaks=bw, plot=F)$breaks
    #Bin width
    band.width<-diff(bin.breaks)[1]
    #Number of bins
    bin.n<-bw
  }

  #Counting the number of elements per bin
  histx<-hist(x, breaks=bin.breaks, plot=FALSE)[[2]]
  histy<-hist(y, breaks=bin.breaks, plot=FALSE)[[2]]
  #Relative counts
  rel.histx<-histx/sum(histx)
  rel.histy<-histy/sum(histy)

  #Calculating the Bhattacharyya Coefficient (sum of the square root of the multiple of the relative counts of both distributions)
  bhatt.coeff<-sum(sqrt(rel.histx*rel.histy))
  return(bhatt.coeff)
  #End
}
