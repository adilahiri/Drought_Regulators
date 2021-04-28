#Function to normalize the data
minmax_normalize <- function(x) {
  
  return ((x - min(x)) / (max(x) - min(x)))
}