source("reCAT_master/R/get_hmm.R")
library(doParallel)

get_start <- function(bayes_score, mean_score, ordIndex, cls_num, rdata, nthread = 3)
{
	cl <- makeCluster(nthread)
	registerDoParallel(cl)
	le = length(ordIndex)
	start = 0
	p = -Inf
	
	for(i in 1:le)
	{
	  if (i == le)
	  {
	    myord = c(le:1)
	  }
		myord = c(i:1, length(ordIndex):i)
		myord <- head(myord, -1)
		re = get_hmm_order(bayes_score, mean_score, ordIndex, cls_num, myord, rdata)
		re <- as.data.frame(re)
		if (max(re) > p)
		{
		  start = i
		  p = max(re)
		}
	}

	#p = apply(log_lk, 2, max)
	#start = which(p == max(p))

	return(start)
}
