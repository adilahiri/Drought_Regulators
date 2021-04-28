#Function to calculate the shape parameters and eventually the node parameters
Shap_Param_Calc<-function(alpha_init,beta_init,Node_Data,Parents) {
  # Case when no parents are there 
  if (is.null(Parents)){
    alpha1 <- alpha_init+ length(which(Node_Data==1))
    beta1 <- beta_init + length(which(Node_Data==0))
    return(cbind(alpha1,beta1))
  }
  
  # Case with one parent
  else if(dim(Parents)[2] == 1 ){
    Truth_Table_One <- expand.grid(c(0,1))#Create a single variable truth table to enlist all the possiblities
    alpha_beta_matrix<-matrix(nrow=2,ncol=2)
    for  (iter in 1:nrow(Truth_Table_One)){
      alpha_val=alpha_init+length(which(Parents==Truth_Table_One[iter,]& Node_Data ==1))
      beta_val =beta_init + length(which(Parents==Truth_Table_One[iter,] & Node_Data ==0))
      alpha_beta_matrix[iter,]=c(alpha_val,beta_val)
    }
    return(alpha_beta_matrix)
  }
  
  # Case with two parents
  else if(dim(Parents)[2] == 2 ){
    Truth_Table_two <- expand.grid(c(0,1),c(0,1))
    alpha_beta_matrix<-matrix(nrow=4,ncol=2)
    for  (iter in 1:nrow(Truth_Table_two)){
      alpha_val=alpha_init+length(which(Parents[,1]==Truth_Table_two[iter,1]& 
                                          Parents[,2]==Truth_Table_two[iter,2]& Node_Data ==1))
      
      beta_val =beta_init + length(which(Parents[,1]==Truth_Table_two[iter,1]& 
                                           Parents[,2]==Truth_Table_two[iter,2] & Node_Data ==0))
      alpha_beta_matrix[iter,]=c(alpha_val,beta_val)
    }
    return(alpha_beta_matrix)
    
  }
  
  
  # Case with three parents
  else if(dim(Parents)[2] == 3 ){
    Truth_Table_three <- expand.grid(c(0,1),c(0,1),c(0,1))
    alpha_beta_matrix<-matrix(nrow=8,ncol=2)
    for  (iter in 1:nrow(Truth_Table_three)){
      alpha_val=alpha_init+length(which(Parents[,1]==Truth_Table_three[iter,1]& 
                                          Parents[,2]==Truth_Table_three[iter,2]& 
                                          Parents[,3]==Truth_Table_three[iter,3]&
                                          Node_Data ==1))
      
      beta_val =beta_init + length(which(Parents[,1]==Truth_Table_three[iter,1]& 
                                           Parents[,2]==Truth_Table_three[iter,2]&
                                           Parents[,3]==Truth_Table_three[iter,3]& 
                                           Node_Data ==0))
      alpha_beta_matrix[iter,]=c(alpha_val,beta_val)
    }
    return(alpha_beta_matrix)
    
  }
  
  
  # Case with four parents
  else if(dim(Parents)[2] == 4 ){
    Truth_Table_four <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1))
    alpha_beta_matrix<-matrix(nrow=16,ncol=2)
    for  (iter in 1:nrow(Truth_Table_four)){
      alpha_val=alpha_init+length(which(Parents[,1]==Truth_Table_four[iter,1]& 
                                          Parents[,2]==Truth_Table_four[iter,2]& 
                                          Parents[,3]==Truth_Table_four[iter,3]&
                                          Parents[,4]==Truth_Table_four[iter,4]&
                                          Node_Data ==1))
      
      beta_val =beta_init + length(which(Parents[,1]==Truth_Table_four[iter,1]& 
                                           Parents[,2]==Truth_Table_four[iter,2]&
                                           Parents[,3]==Truth_Table_four[iter,3]& 
                                           Parents[,4]==Truth_Table_four[iter,4]&
                                           Node_Data ==0))
      alpha_beta_matrix[iter,]=c(alpha_val,beta_val)
    }
    return(alpha_beta_matrix)
    
  }
  
  
  
  # case with five parents 
  else if(dim(Parents)[2] == 5 ){
    Truth_Table_five <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
    alpha_beta_matrix<-matrix(nrow=32,ncol=2)
    for  (iter in 1:nrow(Truth_Table_five)){
      alpha_val=alpha_init+length(which(Parents[,1]==Truth_Table_five[iter,1]& 
                                          Parents[,2]==Truth_Table_five[iter,2]& 
                                          Parents[,3]==Truth_Table_five[iter,3]&
                                          Parents[,4]==Truth_Table_five[iter,4]&
                                          Parents[,5]==Truth_Table_five[iter,5]&
                                          Node_Data ==1))
      
      beta_val =beta_init + length(which(Parents[,1]==Truth_Table_five[iter,1]& 
                                           Parents[,2]==Truth_Table_five[iter,2]&
                                           Parents[,3]==Truth_Table_five[iter,3]& 
                                           Parents[,4]==Truth_Table_five[iter,4]&
                                           Parents[,5]==Truth_Table_five[iter,5]&
                                           Node_Data ==0))
      alpha_beta_matrix[iter,]=c(alpha_val,beta_val)
    }
    return(alpha_beta_matrix)
    
  }
  
  
  
  
  
  # Case with six parents
  else if(dim(Parents)[2] == 6 ){
    Truth_Table_six <- expand.grid(c(0,1),c(0,1),c(0,1),c(0,1),c(0,1),c(0,1))
    alpha_beta_matrix<-matrix(nrow=64,ncol=2)
    for  (iter in 1:nrow(Truth_Table_six)){
      alpha_val=alpha_init+length(which(Parents[,1]==Truth_Table_six[iter,1]& 
                                          Parents[,2]==Truth_Table_six[iter,2]& 
                                          Parents[,3]==Truth_Table_six[iter,3]&
                                          Parents[,4]==Truth_Table_six[iter,4]&
                                          Parents[,5]==Truth_Table_six[iter,5]&
                                          Parents[,6]==Truth_Table_six[iter,6]&
                                          Node_Data ==1))
      
      beta_val =beta_init + length(which(Parents[,1]==Truth_Table_six[iter,1]& 
                                           Parents[,2]==Truth_Table_six[iter,2]&
                                           Parents[,3]==Truth_Table_six[iter,3]& 
                                           Parents[,4]==Truth_Table_six[iter,4]&
                                           Parents[,5]==Truth_Table_six[iter,5]&
                                           Parents[,6]==Truth_Table_six[iter,6]&
                                           Node_Data ==0))
      alpha_beta_matrix[iter,]=c(alpha_val,beta_val)
    }
    return(alpha_beta_matrix)
    
  }
  
  
  # Case with nine parents
  else if(dim(Parents)[2] == 9 ){
    Truth_Table_nine <- expand.grid(c(0,1),c(0,1),c(0,1),
                                   c(0,1),c(0,1),c(0,1),
                                   c(0,1),c(0,1),c(0,1))
    alpha_beta_matrix<-matrix(nrow=512,ncol=2)
    for  (iter in 1:nrow(Truth_Table_nine)){
      alpha_val=alpha_init+length(which(Parents[,1]==Truth_Table_nine[iter,1]& 
                                          Parents[,2]==Truth_Table_nine[iter,2]& 
                                          Parents[,3]==Truth_Table_nine[iter,3]&
                                          Parents[,4]==Truth_Table_nine[iter,4]&
                                          Parents[,5]==Truth_Table_nine[iter,5]&
                                          Parents[,6]==Truth_Table_nine[iter,6]&
                                          Parents[,7]==Truth_Table_nine[iter,7]&
                                          Parents[,8]==Truth_Table_nine[iter,8]&
                                          Parents[,9]==Truth_Table_nine[iter,9]&
                                          Node_Data ==1))
      
      beta_val =beta_init + length(which(Parents[,1]==Truth_Table_nine[iter,1]& 
                                           Parents[,2]==Truth_Table_nine[iter,2]&
                                           Parents[,3]==Truth_Table_nine[iter,3]& 
                                           Parents[,4]==Truth_Table_nine[iter,4]&
                                           Parents[,5]==Truth_Table_nine[iter,5]&
                                           Parents[,6]==Truth_Table_nine[iter,6]&
                                           Parents[,7]==Truth_Table_nine[iter,7]&
                                           Parents[,8]==Truth_Table_nine[iter,8]&
                                           Parents[,9]==Truth_Table_nine[iter,9]&
                                           Node_Data ==0))
      alpha_beta_matrix[iter,]=c(alpha_val,beta_val)
    }
    return(alpha_beta_matrix)
    
  }
  
  
  
}