#Function to binarize the data by mean or median
binarize_mean_median <- function (Mat_Data,method){
  
  if(method=="mean")
    means<-apply(Mat_Data,2,mean)
  else
    means<-apply(Mat_Data,2,median)
  for(iter in 1:ncol(Mat_Data)){
    
    pos0<- which(Mat_Data[,iter]< means[iter])
    pos1<- which(Mat_Data[,iter]>= means[iter])
    Mat_Data[pos0,iter] <-0
    Mat_Data[pos1,iter] <-1
  }
  return(Mat_Data)
}