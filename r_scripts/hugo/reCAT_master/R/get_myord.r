get_myord <- function(cycle_start, ordIndex)
{
  # Checks the direction of the cycle (the order and arranges it accordingly)
  
  # Looks up if the start of the cycle is larger than, ro equal ot half the length of the ordIndex
  # If so, it sets the direction in reverse direction, if not, it is set to forward direction
  
	#if(cycle_start >= length(ordIndex)/2)
	#{
    # Reverse_direction
	# myord = c(cycle_start:1, length(ordIndex):cycle_start)
	# myord <- head(myord, -1)
	#} else {
    # Forward direction
    myord = c(1:cycle_start, cycle_start:length(ordIndex))
    myord <- head(myord, -1)
	}
  return(myord)
}