
## 'boxplot.na' is a function that creates boxplots of percent missing values and provides boxplot statistics
## unfilled - a dataframe before feature filling in xcms, with missing values represented as NA
## biosamples - chacater string of names corresponding to the biological samples in the unfilled dataframe.
## high - character string of feature names corresponding to the high quality features. 
## low - character string of feature names corresponding to the low quality features.

boxplot.na <- function(unfilled,biosamples,high,low){
  missing <- apply(unfilled[,biosamples],1,function(x) mean(is.na(x)))
  names(missing) <- rownames(unfilled)
  bp.low <- boxplot(missing[names(missing)%in%low]*100)
  bp.high <- boxplot(missing[names(missing)%in%high]*100)
  boxplot(missing[names(missing)%in%c(high,low)]*100 ~ as.factor(names(missing[names(missing)%in%c(high,low)])%in%high),
          notch=T,names=c('Low Quality','High Quality'),main='Box Plot of Percent Missing Values',
          cex.lab=1.4,cex.main=1.5,ylab='Percent missing',cex.axis=1.4,varwidth=T)
  return(list(LowQuality=bp.low$stats,HighQuality=bp.high$stats))
}


## 'filter.na' is a function that returns features that pass a provided percent missing threshold
## threshold - fraction indicating the proportion of missing values threshold. Features with a proportion of missing values less than the threshold are retained
filter.na <- function(threshold,unfilled,biosamples){
  missing <- apply(unfilled[,biosamples],1,function(x) mean(is.na(x)))
  names(missing) <- rownames(unfilled)
  return(names(missing[missing<threshold]))
}

## 'diff.miss.fish' is a function that calculates fisher exact p-values for each feature to test for an association between missing values and the biology of interest
## unfilled - a dataframe before feature filling in xcms, with missing values represented as NA
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


## 'filter.fish' is a function that returns features that pass the threshold for fisher exact p-values, i.e., features where there is apparent association between missingness and the biology of interest
## threshold - percentile of the p-value distribution. Features that have p-values below this percentile are retained.
## fish.pvals - pvalues from running 'diff.miss.fish'

filter.fish <- function(threshold,fish.pvals){
  names(fish.pvals[fish.pvals<quantile(fish.pvals,threshold)])
}

