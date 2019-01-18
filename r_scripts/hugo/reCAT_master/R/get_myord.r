get_myord <- function(cycle_start, ordIndex)
{
	# Checks the direction of the cycle (the order and arranges it accordingly)
  
	# Looks up if the start of the cycle is larger than, ro equal ot half the length of the ordIndex
	# If so, it sets the direction in reverse direction, if not, it is set to forward direction
  
	myord = c(cycle_start:length(ordIndex), 1:cycle_start)
	myord <- head(myord, -1)
	return(myord)
}