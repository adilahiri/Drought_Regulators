source('minmax_normalize.R')
source('binarize_mean_median.R')
source('shap_Param_Calc.R')
source('rename_matrix.R')

library(dplyr)
library(ggplot2)
library(plotly)

# Load the Data and prepare the data matrix 
data_set<- as.matrix(read.csv("GSE_Subset.csv",sep=",",header=FALSE))
Gene_Names<- data_set[1,] # Get the Gene Names
data_set<-data_set[-c(1,2),] # Delete the Gene Names and Accession Ids 

data_matrix<-apply(data_set,2,as.numeric) # Convert the dataset to numeric from character
colnames(data_matrix)<-Gene_Names

toDelete_wet <- seq(1, length(data_matrix), 2)
data_matrix <-  data_matrix[-toDelete_wet, ]

# Data matrix will be normalized and then discretized
data_norm<- apply(data_matrix,2,minmax_normalize)
data_binary <- binarize_mean_median(data_matrix,method="mean")

# Visualize the binary data and get the activation and inhibition counts
Activ_Inhib_Mat <- matrix(nrow=32,ncol=3)
colnames(Activ_Inhib_Mat)<-c("Inhibit","Active","total")

for (iter in 1:32){
  Activ_Inhib_Mat [iter,1] <- table(data_binary[,iter])[1]
  Activ_Inhib_Mat [iter,2] <- table(data_binary[,iter])[2]
  Activ_Inhib_Mat [iter,3] <- table(data_binary[,iter])[1]+
    table(data_binary[,iter])[2]
  
}

Activ_Inhib_Mat<-as.data.frame(Activ_Inhib_Mat)
Activ_Inhib_Mat$GeneNames<-Gene_Names

p1<-plot_ly(Activ_Inhib_Mat,x=~GeneNames,y=~Active,type='bar',
            name='Activated')%>% 
  add_trace(y=~Inhibit, name='Inhibited')%>%
  layout(xaxis=list(title='Nodes'),yaxis=list(title='count'),barmode='stack')

show(p1)

# Now we will Calculate the conditional probabilities
# associated with the graph model 

Nodes <- 32
obs<- nrow(data_binary)

alpha_init <- 1
beta_init <- 1
shape_param <- list()
Parents<-matrix()
shape_mat<- matrix(ncol=2)

for (iter_node in 1:Nodes){
  if (iter_node>=1 & iter_node <= 13)
    parents <- NULL
  
  else if (iter_node == 14)## MKK4
    parents<-cbind(data_binary[,13]) ## Depends on MAP3K15
  
  else if (iter_node == 15) ## MPK6
    parents<-cbind(data_binary[,14]) ## Depends on MKK4
  else if (iter_node == 16) ## WRKY59
    parents<-cbind(data_binary[,15]) ## Depends on MPK6)
  
  else if(iter_node ==17) ## RD20 
    parents<-cbind(data_binary[,2]) ## Depends on ANAC072
  
  else if (iter_node == 18) ## ICE1
    parents<-cbind(data_binary[,4],data_binary[,5]) ## Depends on HOS1 and SIZ1
  
  else if (iter_node == 19) ## MYB15
    parents<-cbind(data_binary[,18]) ## Depends on ICE1
  
  else if (iter_node == 20) ## DREB1A
    parents<-cbind(data_binary[,18],data_binary[,19]) ## Depends on ICE1, MYB15
  else if (iter_node == 21) ## DREB1B
    parents<-cbind(data_binary[,18],data_binary[,19]) ## Depends on ICE1, MYB15
  else if (iter_node == 22) ## DREB1C
    parents<-cbind(data_binary[,18],data_binary[,19]) ## Depends on ICE1, MYB15
  
  else if (iter_node == 23) ## WRKY60
    parents<-cbind(data_binary[,6],data_binary[,7]) ## Depends on WRKY18, WRKY40
  
  else if (iter_node == 24) ## MYB2
    parents<-cbind(data_binary[,23]) ## Depends on WRKY60
  
  else if (iter_node == 25) ## ANAC019
    parents<-cbind(data_binary[,1],data_binary[,24]) ## Depends on MYC2 , MYB2
  else if (iter_node == 26) ## ANAC055
    parents<-cbind(data_binary[,1],data_binary[,24]) ## Depends on MYC2 , MYB2
  else if (iter_node == 27) ## RD22
    parents<-cbind(data_binary[,1],data_binary[,24]) ## Depends on MYC2 , MYB2
  else if (iter_node == 28) ## ATAF1
    parents<-cbind(data_binary[,1],data_binary[,24]) ## Depends on MYC2 , MYB2
  
  else if (iter_node==29) ## DREB1D or CBF4
    parents <- cbind(data_binary[,8],data_binary[,9],data_binary[,10],data_binary[,11]) # AREB1,AREB2,ABF1,ABF3
 
   else if (iter_node==30) ## ERD1 
    parents <- cbind(data_binary[,2],data_binary[,3],data_binary[,25],
                     data_binary[,26],data_binary[,28]) #ANAC072,ZFHD1,ANAC019,ANAC055,ATAF1
   
   else if (iter_node ==31) ##DREB2A
     parents <- cbind(data_binary[,8],data_binary[,9],
                      data_binary[,10],data_binary[,11],
                      data_binary[,12],data_binary[,16]) #AREB1,AREB2,ABF1,ABF3, DRIP1,WRK59
   
   else if (iter_node ==32) ##RD29A
     parents <- cbind(data_binary[,6],data_binary[,7],data_binary[,20], 
                      data_binary[,21],data_binary[,22],data_binary[,23],
                      data_binary[,28],data_binary[,29],data_binary[,31]) #WRKY18, WRKY40,DREB1A,DREB1B,DREB1C,
                                                                    #WRKY60,ATAF1,DREB1D,DREB2A

  
  shape_param [[iter_node]] <- Shap_Param_Calc(alpha_init,beta_init,
                                               data_binary[,iter_node],parents)
  initial_index<-dim(shape_mat)[1]
  shape_mat <- rbind(shape_mat,shape_param[[iter_node]])
  diff_index <- initial_index+1:(dim(shape_mat)[1]-initial_index)
  shape_mat <-rename_matrix(shape_mat,Gene_Names[iter_node],diff_index)

}
shape_mat <-shape_mat[-1,]

Expected_Values <- signif(cbind(shape_mat[,1]/(shape_mat[,1]+shape_mat[,2])),digits = 3)

shape_mat2<-cbind(shape_mat,Expected_Values,(1-Expected_Values))

colnames(shape_mat2)<-c("alpha","beta","Activated","Inhibited")






## Define the CPT Tables for each nodes
status <- c("yes", "no")
cptMYC2 <- matrix(c(shape_mat2[1,3], shape_mat2[1,4]), ncol=2, dimnames=list(NULL, status))
cptANAC072 <-matrix(c(shape_mat2[2,3], shape_mat2[2,4]), ncol=2, dimnames=list(NULL, status))
cptZFHD1<-matrix(c(shape_mat2[3,3], shape_mat2[3,4]), ncol=2, dimnames=list(NULL, status))
cptHOS1<- matrix(c(shape_mat2[4,3], shape_mat2[4,4]), ncol=2, dimnames=list(NULL, status))
cptSIZ1<- matrix(c(shape_mat2[5,3], shape_mat2[5,4]), ncol=2, dimnames=list(NULL, status))
cptWRKY18<- matrix(c(shape_mat2[6,3], shape_mat2[6,4]), ncol=2, dimnames=list(NULL, status))
cptWRKY40<- matrix(c(shape_mat2[7,3], shape_mat2[7,4]), ncol=2, dimnames=list(NULL, status))
cptAREB1<- matrix(c(shape_mat2[8,3], shape_mat2[8,4]), ncol=2, dimnames=list(NULL, status))
cptAREB2 <- matrix(c(shape_mat2[9,3], shape_mat2[9,4]), ncol=2, dimnames=list(NULL, status))
cptABF1 <- matrix(c(shape_mat2[10,3], shape_mat2[10,4]), ncol=2, dimnames=list(NULL, status))
cptABF3 <- matrix(c(shape_mat2[11,3], shape_mat2[11,4]), ncol=2, dimnames=list(NULL, status))
cptDRIP1 <- matrix(c(shape_mat2[12,3], shape_mat2[12,4]), ncol=2, dimnames=list(NULL, status))
cptMAP3K15 <- matrix(c(shape_mat2[13,3], shape_mat2[13,4]), ncol=2, dimnames=list(NULL, status))

# MKK4
cptMKK4<-c(shape_mat2[14,3],shape_mat2[14,4],shape_mat2[15,3],shape_mat2[15,4])
dim(cptMKK4)<-c(2,2)
dimnames(cptMKK4)=list("MKK4"=status,"MAP3K15"=status)

# MPK6
cptMPK6<-c(shape_mat2[16,3],shape_mat2[16,4],shape_mat2[17,3],shape_mat2[17,4])
dim(cptMPK6)<-c(2,2)
dimnames(cptMPK6)=list("MPK6"=status,"MKK4"=status)

# WRKY59
cptWRKY59<-c(shape_mat2[18,3],shape_mat2[18,4],shape_mat2[19,3],shape_mat2[19,4])
dim(cptWRKY59)<-c(2,2)
dimnames(cptWRKY59)=list("WRKY59"=status,"MPK6"=status)

#RD20
cptRD20<-c(shape_mat2[20,3],shape_mat2[20,4],shape_mat2[21,3],shape_mat2[21,4])
dim(cptRD20)<-c(2,2)
dimnames(cptRD20)=list("RD20"=status,"ANAC072"=status)

#ICE1
cptICE1<-c(shape_mat2[22,3],shape_mat2[22,4],shape_mat2[23,3],shape_mat2[23,4],
           shape_mat2[24,3],shape_mat2[24,4],shape_mat2[25,3],shape_mat2[25,4])
dim(cptICE1) = c(2, 2, 2)
dimnames(cptICE1) = list("ICE1" =status, "HOS1" =status, "SIZ1" =status)



#MYB15
cptMYB15<- c(shape_mat2[26,3],shape_mat2[26,4],shape_mat2[27,3],shape_mat2[27,4])
dim(cptMYB15)=c(2,2)
dimnames(cptMYB15) =list("MYB15"=status, "ICE1"=status)

#DREB1A
cptDREB1A<-c(t(shape_mat2[28:31,3:4]))
dim(cptDREB1A) = c(2, 2, 2)
dimnames(cptDREB1A) = list("DREB1A" =status, "ICE1" =status, "MYB15" =status)

#DREB1B
cptDREB1B<-c(t(shape_mat2[32:35,3:4]))
dim(cptDREB1B) = c(2, 2, 2)
dimnames(cptDREB1B) = list("DREB1B" =status, "ICE1" =status, "MYB15" =status)

#DREB1C
cptDREB1C<-c(t(shape_mat2[36:39,3:4]))
dim(cptDREB1C) = c(2, 2, 2)
dimnames(cptDREB1C) = list("DREB1C" =status, "ICE1" =status, "MYB15" =status)


#WRKY60
cptWRKY60<-c(t(shape_mat2[40:43,3:4]))
dim(cptWRKY60) = c(2, 2, 2)
dimnames(cptWRKY60) = list("WRKY60" =status, "WRKY18" =status, "WRKY40" =status)


#MYB2
cptMYB2 <- c(shape_mat2[44,3],shape_mat2[44,4],shape_mat2[45,3],shape_mat2[45,4])
dim(cptMYB2)=c(2,2)
dimnames(cptMYB2)= list("MYB2"=status, "WRKY60"=status)

#ANAC019
cptANAC019 <- c(t(shape_mat2[46:49,3:4]))
dim(cptANAC019) = c(2, 2, 2)
dimnames(cptANAC019) = list("ANAC019" =status, "MYC2" =status, "MYB2" =status)

#ANAC055 
cptANAC055 <- c(t(shape_mat2[50:53,3:4]))
dim(cptANAC055) = c(2, 2, 2)
dimnames(cptANAC055) = list("ANAC055" =status, "MYC2" =status, "MYB2" =status)

#RD22
cptRD22 <- c(t(shape_mat2[54:57,3:4]))
dim(cptRD22) = c(2, 2, 2)
dimnames(cptRD22) = list("RD22" =status, "MYC2" =status, "MYB2" =status)

#ATAF1
cptATAF1 <- c(t(shape_mat2[58:61,3:4]))
dim(cptATAF1) = c(2, 2, 2)
dimnames(cptATAF1) = list("ATAF1" =status, "MYC2" =status, "MYB2" =status)

#DREB1D
cptDREB1D <-c(t(shape_mat2[62:77,3:4]))
dim(cptDREB1D) = c(2,2,2,2,2)
dimnames(cptDREB1D) = list("DREB1D" =status, "AREB1" =status, "AREB2" =status,
                           "ABF1" =status, "ABF3" =status)
#ERD1
cptERD1 <-c(t(shape_mat2[78:109,3:4]))
dim(cptERD1) = c(2,2,2,2,2,2)
dimnames(cptERD1) = list("ERD1" =status, "ANAC072" =status, "ZFHD1" =status,
                           "ANAC019" =status, "ANAC055" =status,"ATAF1"=status)

# DREB2A
cptDREB2A <-c(t(shape_mat2[110:173,3:4]))
dim(cptDREB2A) = c(2,2,2,2,2,2,2)
dimnames(cptDREB2A) = list("DREB2A" =status, "AREB1" =status, "AREB2" =status,
                         "ABF1" =status, "ABF3" =status,"DRIP1"=status, "WRKY59"=status)

#RD29A
cptRD29A <-c(t(shape_mat2[174:685,3:4]))
dim(cptRD29A) = c(2,2,2,2,2,2,2,2,2,2)
dimnames(cptRD29A) = list("RD29A" =status, "WRKY18" =status, "WRKY40" =status,
                           "DREB1A" =status, "DREB1B" =status,"DREB1C"=status, 
                           "WRKY60"=status, "ATAF1"=status, "DREB1D"=status, 
                            "DREB2A"=status)


### Initialize and setup the network

library(bnlearn)
network_model<-model2network("[MYC2][ANAC072][ZFHD1][HOS1][SIZ1][WRKY18][WRKY40][AREB1][AREB2][ABF1][ABF3][DRIP1][MAP3K15][MKK4|MAP3K15][MPK6|MKK4][WRKY59|MPK6][RD20|ANAC072][ICE1|HOS1:SIZ1][MYB15|ICE1][DREB1A|ICE1:MYB15][DREB1B|ICE1:MYB15][DREB1C|ICE1:MYB15][WRKY60|WRKY18:WRKY40][MYB2|WRKY60][ANAC019|MYC2:MYB2][ANAC055|MYC2:MYB2][RD22|MYC2:MYB2][ATAF1|MYC2:MYB2][DREB1D|AREB1:AREB2:ABF1:ABF3][ERD1|ANAC072:ZFHD1:ANAC019:ANAC055:ATAF1][DREB2A|AREB1:AREB2:ABF1:ABF3:DRIP1:WRKY59][RD29A|WRKY18:WRKY40:DREB1A:DREB1B:DREB1C:WRKY60:ATAF1:DREB1D:DREB2A]")
par(cex=2)
show(graphviz.plot(network_model))


# Now lets create the BN Learn object and initialize the graph with probabilities
dfit <- custom.fit(network_model, dist=list(MYC2=cptMYC2,ANAC072=cptANAC072,
                                            ZFHD1=cptZFHD1,HOS1=cptHOS1,SIZ1=cptSIZ1,
                                            WRKY18=cptWRKY18,WRKY40=cptWRKY40,AREB1=cptAREB1,
                                            AREB2=cptAREB2,ABF1=cptABF1,ABF3=cptABF3,DRIP1=cptDRIP1,
                                            MAP3K15=cptMAP3K15,MKK4=cptMKK4,MPK6=cptMPK6,WRKY59=cptWRKY59,
                                            RD20=cptRD20,ICE1=cptICE1,MYB15=cptMYB15,DREB1A=cptDREB1A,
                                            DREB1B=cptDREB1B,DREB1C=cptDREB1C,WRKY60=cptWRKY60,
                                            MYB2=cptMYB2,ANAC019=cptANAC019,ANAC055=cptANAC055,
                                            RD22=cptRD22,ATAF1=cptATAF1,DREB1D=cptDREB1D,
                                            ERD1=cptERD1,DREB2A=cptDREB2A,RD29A=cptRD29A))



# Likelihood weighting for inference. 
k=600000
Probability_Matrix_Active<-matrix(nrow=28,ncol=4)
Probability_Matrix_In<-matrix(nrow=28,ncol=4)

colnames(Probability_Matrix_Active)<- c("RD29A","RD20","RD22","ERD1")
colnames(Probability_Matrix_In)<- c("RD29A","RD20","RD22","ERD1")

remove<-c("RD29A","RD20","RD22","ERD1")
Interevention_List<-Gene_Names
Interevention_List<-Interevention_List[! Interevention_List %in% remove]
rownames(Probability_Matrix_Active)<-Interevention_List
rownames(Probability_Matrix_In)<-Interevention_List


iter_RD29A=2
iter_RD20=2
iter_RD22=2
iter_ERD1=1

set.seed(4)

Probability_Matrix_Active[1,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MYC2="yes"),method='lw',n=k)
Probability_Matrix_Active[1,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MYC2="yes"),method='lw',n=k)
Probability_Matrix_Active[1,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MYC2="yes"),method='lw',n=k)
Probability_Matrix_Active[1,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MYC2="yes"),method='lw',n=k)

Probability_Matrix_Active[2,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ANAC072="yes"),method='lw',n=k)
Probability_Matrix_Active[2,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ANAC072="yes"),method='lw',n=k)
Probability_Matrix_Active[2,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ANAC072="yes"),method='lw',n=k)
Probability_Matrix_Active[2,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ANAC072="yes"),method='lw',n=k)


Probability_Matrix_Active[3,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ZFHD1="yes"),method='lw',n=k)
Probability_Matrix_Active[3,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ZFHD1="yes"),method='lw',n=k)
Probability_Matrix_Active[3,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ZFHD1="yes"),method='lw',n=k)
Probability_Matrix_Active[3,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ZFHD1="yes"),method='lw',n=k)


Probability_Matrix_Active[4,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(HOS1="yes"),method='lw',n=k)
Probability_Matrix_Active[4,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(HOS1="yes"),method='lw',n=k)
Probability_Matrix_Active[4,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(HOS1="yes"),method='lw',n=k)
Probability_Matrix_Active[4,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(HOS1="yes"),method='lw',n=k)


Probability_Matrix_Active[5,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(SIZ1="yes"),method='lw',n=k)
Probability_Matrix_Active[5,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(SIZ1="yes"),method='lw',n=k)
Probability_Matrix_Active[5,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(SIZ1="yes"),method='lw',n=k)
Probability_Matrix_Active[5,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(SIZ1="yes"),method='lw',n=k)


Probability_Matrix_Active[6,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(WRKY18="yes"),method='lw',n=k)
Probability_Matrix_Active[6,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(WRKY18="yes"),method='lw',n=k)
Probability_Matrix_Active[6,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(WRKY18="yes"),method='lw',n=k)
Probability_Matrix_Active[6,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(WRKY18="yes"),method='lw',n=k)


Probability_Matrix_Active[7,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(WRKY40="yes"),method='lw',n=k)
Probability_Matrix_Active[7,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(WRKY40="yes"),method='lw',n=k)
Probability_Matrix_Active[7,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(WRKY40="yes"),method='lw',n=k)
Probability_Matrix_Active[7,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(WRKY40="yes"),method='lw',n=k)


Probability_Matrix_Active[8,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(AREB1="yes"),method='lw',n=k)
Probability_Matrix_Active[8,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(AREB1="yes"),method='lw',n=k)
Probability_Matrix_Active[8,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(AREB1="yes"),method='lw',n=k)
Probability_Matrix_Active[8,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(AREB1="yes"),method='lw',n=k)


Probability_Matrix_Active[9,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(AREB2="yes"),method='lw',n=k)
Probability_Matrix_Active[9,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(AREB2="yes"),method='lw',n=k)
Probability_Matrix_Active[9,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(AREB2="yes"),method='lw',n=k)
Probability_Matrix_Active[9,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(AREB2="yes"),method='lw',n=k)


Probability_Matrix_Active[10,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ABF1="yes"),method='lw',n=k)
Probability_Matrix_Active[10,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ABF1="yes"),method='lw',n=k)
Probability_Matrix_Active[10,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ABF1="yes"),method='lw',n=k)
Probability_Matrix_Active[10,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ABF1="yes"),method='lw',n=k)


Probability_Matrix_Active[11,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ABF3="yes"),method='lw',n=k)
Probability_Matrix_Active[11,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ABF3="yes"),method='lw',n=k)
Probability_Matrix_Active[11,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ABF3="yes"),method='lw',n=k)
Probability_Matrix_Active[11,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ABF3="yes"),method='lw',n=k)


Probability_Matrix_Active[12,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DRIP1="yes"),method='lw',n=k)
Probability_Matrix_Active[12,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DRIP1="yes"),method='lw',n=k)
Probability_Matrix_Active[12,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DRIP1="yes"),method='lw',n=k)
Probability_Matrix_Active[12,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DRIP1="yes"),method='lw',n=k)


Probability_Matrix_Active[13,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MAP3K15="yes"),method='lw',n=k)
Probability_Matrix_Active[13,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MAP3K15="yes"),method='lw',n=k)
Probability_Matrix_Active[13,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MAP3K15="yes"),method='lw',n=k)
Probability_Matrix_Active[13,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MAP3K15="yes"),method='lw',n=k)


Probability_Matrix_Active[14,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MKK4="yes"),method='lw',n=k)
Probability_Matrix_Active[14,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MKK4="yes"),method='lw',n=k)
Probability_Matrix_Active[14,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MKK4="yes"),method='lw',n=k)
Probability_Matrix_Active[14,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MKK4="yes"),method='lw',n=k)


Probability_Matrix_Active[15,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MPK6="yes"),method='lw',n=k)
Probability_Matrix_Active[15,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MPK6="yes"),method='lw',n=k)
Probability_Matrix_Active[15,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MPK6="yes"),method='lw',n=k)
Probability_Matrix_Active[15,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MPK6="yes"),method='lw',n=k)


Probability_Matrix_Active[16,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(WRKY59="yes"),method='lw',n=k)
Probability_Matrix_Active[16,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(WRKY59="yes"),method='lw',n=k)
Probability_Matrix_Active[16,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(WRKY59="yes"),method='lw',n=k)
Probability_Matrix_Active[16,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(WRKY59="yes"),method='lw',n=k)


Probability_Matrix_Active[17,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ICE1="yes"),method='lw',n=k)
Probability_Matrix_Active[17,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ICE1="yes"),method='lw',n=k)
Probability_Matrix_Active[17,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ICE1="yes"),method='lw',n=k)
Probability_Matrix_Active[17,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ICE1="yes"),method='lw',n=k)


Probability_Matrix_Active[18,1]=prob_MYC2_act<-cpquery(dfit,(RD29A==status[iter_RD29A]),list(MYB15="yes"),method='lw',n=k)
Probability_Matrix_Active[18,2]=prob_MYC2_act<-cpquery(dfit,(RD20==status[iter_RD20]),list(MYB15="yes"),method='lw',n=k)
Probability_Matrix_Active[18,3]=prob_MYC2_act<-cpquery(dfit,(RD22==status[iter_RD22]),list(MYB15="yes"),method='lw',n=k)
Probability_Matrix_Active[18,4]=prob_MYC2_act<-cpquery(dfit,(ERD1==status[iter_ERD1]),list(MYB15="yes"),method='lw',n=k)


Probability_Matrix_Active[19,1]=prob_MYC2_act<-cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB1A="yes"),method='lw',n=k)
Probability_Matrix_Active[19,2]=prob_MYC2_act<-cpquery(dfit,(RD20==status[iter_RD20]),list(DREB1A="yes"),method='lw',n=k)
Probability_Matrix_Active[19,3]=prob_MYC2_act<-cpquery(dfit,(RD22==status[iter_RD22]),list(DREB1A="yes"),method='lw',n=k)
Probability_Matrix_Active[19,4]=prob_MYC2_act<-cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB1A="yes"),method='lw',n=k)


Probability_Matrix_Active[20,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB1B="yes"),method='lw',n=k)
Probability_Matrix_Active[20,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DREB1B="yes"),method='lw',n=k)
Probability_Matrix_Active[20,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DREB1B="yes"),method='lw',n=k)
Probability_Matrix_Active[20,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB1B="yes"),method='lw',n=k)


Probability_Matrix_Active[21,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB1C="yes"),method='lw',n=k)
Probability_Matrix_Active[21,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DREB1C="yes"),method='lw',n=k)
Probability_Matrix_Active[21,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DREB1C="yes"),method='lw',n=k)
Probability_Matrix_Active[21,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB1C="yes"),method='lw',n=k)


Probability_Matrix_Active[22,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(WRKY60="yes"),method='lw',n=k)
Probability_Matrix_Active[22,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(WRKY60="yes"),method='lw',n=k)
Probability_Matrix_Active[22,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(WRKY60="yes"),method='lw',n=k)
Probability_Matrix_Active[22,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(WRKY60="yes"),method='lw',n=k)


Probability_Matrix_Active[23,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MYB2="yes"),method='lw',n=k)
Probability_Matrix_Active[23,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MYB2="yes"),method='lw',n=k)
Probability_Matrix_Active[23,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MYB2="yes"),method='lw',n=k)
Probability_Matrix_Active[23,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MYB2="yes"),method='lw',n=k)

Probability_Matrix_Active[24,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ANAC019="yes"),method='lw',n=k)
Probability_Matrix_Active[24,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ANAC019="yes"),method='lw',n=k)
Probability_Matrix_Active[24,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ANAC019="yes"),method='lw',n=k)
Probability_Matrix_Active[24,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ANAC019="yes"),method='lw',n=k)


Probability_Matrix_Active[25,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ANAC055="yes"),method='lw',n=k)
Probability_Matrix_Active[25,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ANAC055="yes"),method='lw',n=k)
Probability_Matrix_Active[25,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ANAC055="yes"),method='lw',n=k)
Probability_Matrix_Active[25,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ANAC055="yes"),method='lw',n=k)


Probability_Matrix_Active[26,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ATAF1="yes"),method='lw',n=k)
Probability_Matrix_Active[26,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ATAF1="yes"),method='lw',n=k)
Probability_Matrix_Active[26,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ATAF1="yes"),method='lw',n=k)
Probability_Matrix_Active[26,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ATAF1="yes"),method='lw',n=k)


Probability_Matrix_Active[27,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB1D="yes"),method='lw',n=k)
Probability_Matrix_Active[27,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DREB1D="yes"),method='lw',n=k)
Probability_Matrix_Active[27,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DREB1D="yes"),method='lw',n=k)
Probability_Matrix_Active[27,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB1D="yes"),method='lw',n=k)

Probability_Matrix_Active[28,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB2A="yes"),method='lw',n=k)
Probability_Matrix_Active[28,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DREB2A="yes"),method='lw',n=k)
Probability_Matrix_Active[28,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DREB2A="yes"),method='lw',n=k)
Probability_Matrix_Active[28,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB2A="yes"),method='lw',n=k)


## INACTIVE Case
Probability_Matrix_In[1,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MYC2="no"),method='lw',n=k)
Probability_Matrix_In[1,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MYC2="no"),method='lw',n=k)
Probability_Matrix_In[1,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MYC2="no"),method='lw',n=k)
Probability_Matrix_In[1,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MYC2="no"),method='lw',n=k)

Probability_Matrix_In[2,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ANAC072="no"),method='lw',n=k)
Probability_Matrix_In[2,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ANAC072="no"),method='lw',n=k)
Probability_Matrix_In[2,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ANAC072="no"),method='lw',n=k)
Probability_Matrix_In[2,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ANAC072="no"),method='lw',n=k)


Probability_Matrix_In[3,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ZFHD1="no"),method='lw',n=k)
Probability_Matrix_In[3,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ZFHD1="no"),method='lw',n=k)
Probability_Matrix_In[3,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ZFHD1="no"),method='lw',n=k)
Probability_Matrix_In[3,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ZFHD1="no"),method='lw',n=k)


Probability_Matrix_In[4,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(HOS1="no"),method='lw',n=k)
Probability_Matrix_In[4,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(HOS1="no"),method='lw',n=k)
Probability_Matrix_In[4,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(HOS1="no"),method='lw',n=k)
Probability_Matrix_In[4,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(HOS1="no"),method='lw',n=k)


Probability_Matrix_In[5,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(SIZ1="no"),method='lw',n=k)
Probability_Matrix_In[5,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(SIZ1="no"),method='lw',n=k)
Probability_Matrix_In[5,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(SIZ1="no"),method='lw',n=k)
Probability_Matrix_In[5,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(SIZ1="no"),method='lw',n=k)


Probability_Matrix_In[6,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(WRKY18="no"),method='lw',n=k)
Probability_Matrix_In[6,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(WRKY18="no"),method='lw',n=k)
Probability_Matrix_In[6,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(WRKY18="no"),method='lw',n=k)
Probability_Matrix_In[6,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(WRKY18="no"),method='lw',n=k)


Probability_Matrix_In[7,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(WRKY40="no"),method='lw',n=k)
Probability_Matrix_In[7,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(WRKY40="no"),method='lw',n=k)
Probability_Matrix_In[7,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(WRKY40="no"),method='lw',n=k)
Probability_Matrix_In[7,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(WRKY40="no"),method='lw',n=k)


Probability_Matrix_In[8,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(AREB1="no"),method='lw',n=k)
Probability_Matrix_In[8,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(AREB1="no"),method='lw',n=k)
Probability_Matrix_In[8,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(AREB1="no"),method='lw',n=k)
Probability_Matrix_In[8,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(AREB1="no"),method='lw',n=k)


Probability_Matrix_In[9,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(AREB2="no"),method='lw',n=k)
Probability_Matrix_In[9,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(AREB2="no"),method='lw',n=k)
Probability_Matrix_In[9,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(AREB2="no"),method='lw',n=k)
Probability_Matrix_In[9,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(AREB2="no"),method='lw',n=k)


Probability_Matrix_In[10,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ABF1="no"),method='lw',n=k)
Probability_Matrix_In[10,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ABF1="no"),method='lw',n=k)
Probability_Matrix_In[10,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ABF1="no"),method='lw',n=k)
Probability_Matrix_In[10,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ABF1="no"),method='lw',n=k)


Probability_Matrix_In[11,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ABF3="no"),method='lw',n=k)
Probability_Matrix_In[11,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ABF3="no"),method='lw',n=k)
Probability_Matrix_In[11,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ABF3="no"),method='lw',n=k)
Probability_Matrix_In[11,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ABF3="no"),method='lw',n=k)


Probability_Matrix_In[12,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DRIP1="no"),method='lw',n=k)
Probability_Matrix_In[12,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DRIP1="no"),method='lw',n=k)
Probability_Matrix_In[12,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DRIP1="no"),method='lw',n=k)
Probability_Matrix_In[12,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DRIP1="no"),method='lw',n=k)


Probability_Matrix_In[13,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MAP3K15="no"),method='lw',n=k)
Probability_Matrix_In[13,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MAP3K15="no"),method='lw',n=k)
Probability_Matrix_In[13,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MAP3K15="no"),method='lw',n=k)
Probability_Matrix_In[13,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MAP3K15="no"),method='lw',n=k)


Probability_Matrix_In[14,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MKK4="no"),method='lw',n=k)
Probability_Matrix_In[14,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MKK4="no"),method='lw',n=k)
Probability_Matrix_In[14,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MKK4="no"),method='lw',n=k)
Probability_Matrix_In[14,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MKK4="no"),method='lw',n=k)


Probability_Matrix_In[15,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MPK6="no"),method='lw',n=k)
Probability_Matrix_In[15,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MPK6="no"),method='lw',n=k)
Probability_Matrix_In[15,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MPK6="no"),method='lw',n=k)
Probability_Matrix_In[15,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MPK6="no"),method='lw',n=k)


Probability_Matrix_In[16,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(WRKY59="no"),method='lw',n=k)
Probability_Matrix_In[16,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(WRKY59="no"),method='lw',n=k)
Probability_Matrix_In[16,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(WRKY59="no"),method='lw',n=k)
Probability_Matrix_In[16,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(WRKY59="no"),method='lw',n=k)


Probability_Matrix_In[17,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ICE1="no"),method='lw',n=k)
Probability_Matrix_In[17,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ICE1="no"),method='lw',n=k)
Probability_Matrix_In[17,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ICE1="no"),method='lw',n=k)
Probability_Matrix_In[17,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ICE1="no"),method='lw',n=k)


Probability_Matrix_In[18,1]=prob_MYC2_act<-cpquery(dfit,(RD29A==status[iter_RD29A]),list(MYB15="no"),method='lw',n=k)
Probability_Matrix_In[18,2]=prob_MYC2_act<-cpquery(dfit,(RD20==status[iter_RD20]),list(MYB15="no"),method='lw',n=k)
Probability_Matrix_In[18,3]=prob_MYC2_act<-cpquery(dfit,(RD22==status[iter_RD22]),list(MYB15="no"),method='lw',n=k)
Probability_Matrix_In[18,4]=prob_MYC2_act<-cpquery(dfit,(ERD1==status[iter_ERD1]),list(MYB15="no"),method='lw',n=k)


Probability_Matrix_In[19,1]=prob_MYC2_act<-cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB1A="no"),method='lw',n=k)
Probability_Matrix_In[19,2]=prob_MYC2_act<-cpquery(dfit,(RD20==status[iter_RD20]),list(DREB1A="no"),method='lw',n=k)
Probability_Matrix_In[19,3]=prob_MYC2_act<-cpquery(dfit,(RD22==status[iter_RD22]),list(DREB1A="no"),method='lw',n=k)
Probability_Matrix_In[19,4]=prob_MYC2_act<-cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB1A="no"),method='lw',n=k)


Probability_Matrix_In[20,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB1B="no"),method='lw',n=k)
Probability_Matrix_In[20,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DREB1B="no"),method='lw',n=k)
Probability_Matrix_In[20,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DREB1B="no"),method='lw',n=k)
Probability_Matrix_In[20,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB1B="no"),method='lw',n=k)


Probability_Matrix_In[21,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB1C="no"),method='lw',n=k)
Probability_Matrix_In[21,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DREB1C="no"),method='lw',n=k)
Probability_Matrix_In[21,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DREB1C="no"),method='lw',n=k)
Probability_Matrix_In[21,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB1C="no"),method='lw',n=k)


Probability_Matrix_In[22,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(WRKY60="no"),method='lw',n=k)
Probability_Matrix_In[22,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(WRKY60="no"),method='lw',n=k)
Probability_Matrix_In[22,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(WRKY60="no"),method='lw',n=k)
Probability_Matrix_In[22,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(WRKY60="no"),method='lw',n=k)


Probability_Matrix_In[23,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(MYB2="no"),method='lw',n=k)
Probability_Matrix_In[23,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(MYB2="no"),method='lw',n=k)
Probability_Matrix_In[23,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(MYB2="no"),method='lw',n=k)
Probability_Matrix_In[23,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(MYB2="no"),method='lw',n=k)

Probability_Matrix_In[24,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ANAC019="no"),method='lw',n=k)
Probability_Matrix_In[24,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ANAC019="no"),method='lw',n=k)
Probability_Matrix_In[24,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ANAC019="no"),method='lw',n=k)
Probability_Matrix_In[24,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ANAC019="no"),method='lw',n=k)


Probability_Matrix_In[25,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ANAC055="no"),method='lw',n=k)
Probability_Matrix_In[25,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ANAC055="no"),method='lw',n=k)
Probability_Matrix_In[25,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ANAC055="no"),method='lw',n=k)
Probability_Matrix_In[25,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ANAC055="no"),method='lw',n=k)


Probability_Matrix_In[26,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ATAF1="no"),method='lw',n=k)
Probability_Matrix_In[26,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(ATAF1="no"),method='lw',n=k)
Probability_Matrix_In[26,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(ATAF1="no"),method='lw',n=k)
Probability_Matrix_In[26,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ATAF1="no"),method='lw',n=k)


Probability_Matrix_In[27,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB1D="no"),method='lw',n=k)
Probability_Matrix_In[27,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DREB1D="no"),method='lw',n=k)
Probability_Matrix_In[27,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DREB1D="no"),method='lw',n=k)
Probability_Matrix_In[27,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB1D="no"),method='lw',n=k)

Probability_Matrix_In[28,1]=cpquery(dfit,(RD29A==status[iter_RD29A]),list(DREB2A="no"),method='lw',n=k)
Probability_Matrix_In[28,2]=cpquery(dfit,(RD20==status[iter_RD20]),list(DREB2A="no"),method='lw',n=k)
Probability_Matrix_In[28,3]=cpquery(dfit,(RD22==status[iter_RD22]),list(DREB2A="no"),method='lw',n=k)
Probability_Matrix_In[28,4]=cpquery(dfit,(ERD1==status[iter_ERD1]),list(DREB2A="no"),method='lw',n=k)



Probability_Matrix_Active_df<-as.data.frame(Probability_Matrix_Active)
Probability_Matrix_Active_df$GeneNames<-rownames(Probability_Matrix_Active)
p2 <- plot_ly(Probability_Matrix_Active_df, y = ~GeneNames, x = ~RD29A, type = 'bar', name = 'RD29A',orientation = 'h')%>% 
  add_trace(x = ~RD20, name = 'RD20')%>%
  add_trace(x = ~RD22, name = 'RD22') %>%  
  add_trace(x = ~ERD1, name = 'ERD1') %>% 
  layout(xaxis = list(title = 'probability',range = c(0, 1)), barmode = 'group')

show(p2)

Probability_Matrix_In_df<-as.data.frame(Probability_Matrix_In)
Probability_Matrix_In_df$GeneNames<-rownames(Probability_Matrix_In)
p3 <- plot_ly(Probability_Matrix_In_df, y = ~GeneNames, x = ~RD29A, type = 'bar', name = 'RD29A',orientation = 'h')%>% 
  add_trace(x = ~RD20, name = 'RD20')%>%
  add_trace(x = ~RD22, name = 'RD22') %>%  
  add_trace(x = ~ERD1, name = 'ERD1') %>% 
  layout(xaxis = list(title = 'probability',range = c(0, 1)), barmode = 'group')

scale_metric_active <- matrix(nrow=28,ncol=1)
scale_metric_in <- matrix(nrow=28,ncol=1)
show(p3)

for(iter in 1:28){
  
  scale_metric_active[iter]<-prod((Probability_Matrix_Active[iter,1:4]))
  scale_metric_in[iter]<-prod((Probability_Matrix_In[iter,1:4]))
  
}


Probability_Matrix_Active_df$Scale_Metric<-scale_metric_active
Probability_Matrix_In_df$Scale_Metric<-scale_metric_in

Probability_Matrix_Active_df$GeneNames <- factor(Probability_Matrix_Active_df$GeneNames, 
                                                 levels = unique(Probability_Matrix_Active_df$GeneNames)
                                                 [order(Probability_Matrix_Active_df$Scale_Metric, decreasing = TRUE)])


Probability_Matrix_In_df$GeneNames <- factor(Probability_Matrix_In_df$GeneNames, 
                                                 levels = unique(Probability_Matrix_In_df$GeneNames)
                                                 [order(Probability_Matrix_In_df$Scale_Metric, decreasing = TRUE)])



p4 <- plot_ly(Probability_Matrix_Active_df, x = ~GeneNames, y = ~Scale_Metric,
              name='Activation scores under optimal conditions',
              text=round(Probability_Matrix_Active_df$Scale_Metric,4), textposition='auto',type = 'bar') %>% 
  layout(xaxis=list(title='Nodes'),yaxis = list(title = 'scores',range=c(0,0.1)))
show(p4)

p5 <- plot_ly(Probability_Matrix_In_df, x = ~GeneNames, y = ~Scale_Metric,
              name='Inhibition scores under optimal conditions',
              marker=list(color=c('#ff8000')),
              text=round(Probability_Matrix_In_df$Scale_Metric,4), textposition='auto',type = 'bar') %>% 
  layout(xaxis=list(title='Nodes'),
         yaxis = list(title = 'scores',range=c(0,0.1)))
show(p5)

set.seed(4)
a_combo=cpquery(dfit,(RD29A==status[iter_RD29A]),list(ATAF1="no",MYC2="yes"),method='lw',n=k)
b_combo=cpquery(dfit,(RD20==status[iter_RD20]),list(ATAF1="no",MYC2="yes"),method='lw',n=k)
c_combo=cpquery(dfit,(RD22==status[iter_RD22]),list(ATAF1="no",MYC2="yes"),method = 'lw',n=k)
d_combo=cpquery(dfit,(ERD1==status[iter_ERD1]),list(ATAF1="no",MYC2="yes"),method = 'lw',n=k)
combo_result=prod(a_combo,b_combo,c_combo,d_combo)
combo_bar <- c(combo_result,Probability_Matrix_In_df[26,6],Probability_Matrix_Active_df[1,6])
df_combo <-data.frame (Interventions= c("MYC2=1 & ATAF1=0", "ATAF1=0", "MYC2=1"), Scores= combo_bar)

p6 <- plot_ly(df_combo, x = ~Interventions, y = ~Scores,
              name='Single vs multi node intervention under optimal conditions',
              marker=list(color=c('#2ca02c','#ff7f0e','#1f77b4')), 
              text=round(df_combo$Scores,4), textposition='auto',textfont = list(color = '#000000'),type = 'bar') %>% 
  layout(xaxis=list(title='Nodes'),yaxis = list(title = 'scores',range=c(0,0.12)))
show(p6)
