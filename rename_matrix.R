#Function to name rows with appropriate gene/node names
rename_matrix <- function(shape_mat,Gene_Names,diff_index) {
  
  rownames(shape_mat)[diff_index]<-Gene_Names
  return(shape_mat)
}