
## This function creates a mean-difference plot as described in Schiffman et al. 
## Points in red correspond to high quality features
## Points in black correspond to low quality features
## filled - a filled dataframe (i.e. resulting from the 'fillChromPeaks' function in xcms) of logged feature abundances. The dataframe should have features in the rows and samples in the columns (including blank samples), with the names of features being the rownames of the matrix and the sample names being the column names of the matrix. 
## blanks - chacater string of names corresponding to the blank samples in the fill dataframe.
## biosamples - chacater string of names corresponding to the biological samples in the fill dataframe.
## high - character string of feature names corresponding to the high quality features. 
## low - character string of feature names corresponding to the low quality features.

MDplot <- function(filled,blanks,biosamples,high,low){
  obsMean <- as.vector(apply(filled[,biosamples], 1, mean))
  names(obsMean) <- rownames(filled)
  blankMean <- apply(filled[, blanks], 1, mean)
  names(blankMean) <- rownames(filled)
  smoothScatter((blankMean+obsMean)/2,(obsMean)-(blankMean),
                xlab="Mean",ylab="Difference",main="Mean-Difference Plot",
                cex.lab=1.4,cex.main=1.5)
  abline(h=0,lwd=1,col="blue")
  legend('bottomright',legend=c('Good Quality','Poor Quality'),col=c('red','black'),lwd=2,cex=1.2)

  blanks.low <- blankMean[names(blankMean)%in%low]
  obs.low <- (obsMean)[names(obsMean)%in%low]
  points((blanks.low+obs.low)/2,obs.low-blanks.low,pch=19,cex=.6)
  blanks.high <- blankMean[names(blankMean)%in%high]
  obs.high <- obsMean[names(obsMean)%in%high]
  points((blanks.high+obs.high)/2,obs.high-blanks.high,col='red',pch=19,cex=.7)
  
}

