
## 'boxplot.na' is a function that creates boxplots of percent missing values and provides boxplot statistics
## unfilled - an unfilled dataframe (i.e. before using 'fillChromPeaks' in xcms) of logged feature abundances, with missing values represented as NA. Features should be in the rows and samples in the columns, with rownames equal to feature ID's/names and column names equal to sample names.
## biosamples - chacater string of names corresponding to the biological samples in the unfilled dataframe.
## high - character string of feature names corresponding to the high quality features. 
## low - character string of feature names corresponding to the low quality features.

boxplot.na <- function(unfilled,biosamples,high,low){
  missing <- apply(unfilled[,biosamples],1,function(x) mean(is.na(x)))
  names(missing) <- rownames(unfilled)
  bp.low <- boxplot(missing[names(missing)%in%low]*100,main="Percent missing, Low quality")
  bp.high <- boxplot(missing[names(missing)%in%high]*100,main="Percent missing, High quality")
  boxplot(missing[names(missing)%in%c(high,low)]*100 ~ as.factor(names(missing[names(missing)%in%c(high,low)])%in%high),
          notch=T,names=c('Low Quality','High Quality'),main='Box Plot of Percent Missing Values',
          cex.lab=1.4,cex.main=1.5,ylab='Percent missing',cex.axis=1.4,varwidth=T)
  return(list(LowQuality=bp.low$stats,HighQuality=bp.high$stats))
}

## 'density.na' is a function that creates density plots of percent missing values for high and low quality features

density.na <- function(unfilled,biosamples,high,low){
  missing <- apply(unfilled[,biosamples],1,function(x) mean(is.na(x)))
  names(missing) <- rownames(unfilled)
  plot(density(missing[names(missing)%in%low]*100),main='Density Plot of Percent Missing Values',
          cex.lab=1.4,cex.main=1.5,ylab='Percent missing',cex.axis=1.4,xlab='')
  lines(density(missing[names(missing)%in%high]*100),col='red')
  legend('topleft',col=c('red','black'),lwd=2,legend=c('High Quality','Low Quality'),cex=1.2)
}


## 'filter.na' is a function that returns features that pass a provided percent missing threshold
## threshold - fraction indicating the proportion of missing values threshold. Features with a proportion of missing values less than the threshold are retained
filter.na <- function(threshold,unfilled,biosamples){
  missing <- apply(unfilled[,biosamples],1,function(x) mean(is.na(x)))
  names(missing) <- rownames(unfilled)
  return(names(missing[missing<threshold]))
}

## 'diff.miss.fish' is a function that calculates fisher exact p-values for each feature to test for an association between missing values and the binary biological factor of interest (e.g. case/control)
## unfilled - an unfilled dataframe (i.e. before using 'fillChromPeaks' in xcms) of logged feature abundances, with missing values represented as NA. Features should be in the rows and samples in the columns, with rownames equal to feature ID's/names and column names equal to sample names.
## biosamples - chacater string of names corresponding to the biological samples in the unfilled dataframe.
## bio - factor indicating the binary biological factor of interest

diff.miss.fish <- function(unfilled,biosamples,bio){
  percent.zero <- apply(unfilled[,biosamples],1,function(x) mean(is.na(x)))
  fish.pvals<- c()
  for (i in 1:nrow(unfilled)){
    if(percent.zero[i]==0 | percent.zero[i]==1){
      fish.pvals[i] <- 1
    } else{
      fish.pvals[i] <- fisher.test(as.numeric(is.na(unfilled[i,biosamples])),bio)$p.value
    }
  }
  names(fish.pvals) <- rownames(unfilled)
  return(fish.pvals)
}

## 'diff.miss.chi' is a function that calculates Chi-sqaure p-values for each feature to test for an association between missing values and the multi-level biological factor of interest (e.g. low/medium/high exposure)
## unfilled - an unfilled dataframe (i.e. before using 'fillChromPeaks' in xcms) of logged feature abundances, with missing values represented as NA. Features should be in the rows and samples in the columns, with rownames equal to feature ID's/names and column names equal to sample names.
## biosamples - chacater string of names corresponding to the biological samples in the unfilled dataframe.
## bio - factor indicating the multi-level biological factor of interest


diff.miss.chi <- function(unfilled,biosamples,bio){
  percent.zero <- apply(unfilled[,biosamples],1,function(x) mean(is.na(x)))
  chi.pvals<- c()
  for (i in 1:nrow(unfilled)){
    if(percent.zero[i]==0 | percent.zero[i]==1){
      chi.pvals[i] <- 1
    } else{
      chi.pvals[i] <- chisq.test(table(is.na(unfilled[i,biosamples]),bio))$p.value
    }
  }
  names(chi.pvals) <- rownames(unfilled)
  return(chi.pvals)
}


## 'diff.miss.wilc' is a function that calculates Wilcoxon rank sum p-values for each feature to test for an association between missing values and the continuous biological variable of interest (e.g. BMI)
## unfilled - an unfilled dataframe (i.e. before using 'fillChromPeaks' in xcms) of logged feature abundances, with missing values represented as NA. Features should be in the rows and samples in the columns, with rownames equal to feature ID's/names and column names equal to sample names.
## biosamples - chacater string of names corresponding to the biological samples in the unfilled dataframe.
## bio - numeric, indicating the continuous biological variable of interest


diff.miss.wilc <- function(unfilled,biosamples,bio){
  percent.zero <- apply(unfilled[,biosamples],1,function(x) mean(is.na(x)))
  wilc.pvals<- c()
  for (i in 1:nrow(unfilled)){
    if(percent.zero[i]==0 | percent.zero[i]==1){
      wilc.pvals[i] <- 1
    } else{
      wilc.pvals[i] <- wilcox.test(bio[is.na(unfilled[i,biosamples])],bio[!is.na(unfilled[i,biosamples])])$p.value
    }
  }
  names(wilc.pvals) <- rownames(unfilled)
  return(wilc.pvals)
}


## 'filter.pvals' is a function that returns features that pass the p-value threshold, i.e., features where there is apparent association between missingness and the biology of interest
## threshold - percentile of the p-value distribution. Features that have p-values below this percentile are retained.
## pvals - pvalues from running either 'diff.miss.fish', 'diff.miss.chi' or 'diff.miss.wilc'

filter.pvals <- function(threshold,pvals){
  names(pvals[pvals<quantile(pvals,threshold)])
}

