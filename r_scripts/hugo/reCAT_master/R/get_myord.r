# Made by: Hugo Swenson, 2019-01-25

get_myord <- function(cycle_start, ordIndex)
{
  # Obtains the real time-series
	myord = c(cycle_start:length(ordIndex), 1:cycle_start)
	myord <- head(myord, -1)
	return(myord)
}