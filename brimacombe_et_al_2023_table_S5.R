########################################################
##
##
##
##  For recreating table S5 of Brimacombe et al. (2024):
##  On the nature of structure in collections  
##  of freely available food webs
##
##
##
########################################################

#clearing environment
rm(list = ls())

#loading necessary packages
library(ggplot2)
library(igraph)
library(dplyr)
library(Hmisc)

#set working directory
setwd(".")

#reading in GCD-11 matrix
ecological <- as.matrix(read.csv("gcd11.csv",row.names = 1))

#the associated metadata of all food webs
metadata <- as.matrix(read.csv("metaData.csv"))
metadata <- as.data.frame(metadata)
metadata$author <- as.factor(metadata$author)
metadata$number_of_nodes <- as.numeric(metadata$number_of_nodes)
metadata$year_published <- as.numeric(metadata$year_published)
metadata$decade_published <- as.factor(metadata$decade_published)
metadata$Ecopath <- as.factor(metadata$Ecopath)
metadata$Specialized_aquatic_type <- as.factor(metadata$Specialized_aquatic_type)

########################################################

#storage for mean and standard deviation in pairwise GCD-11
average_distance_per_realization_network_size<- matrix(,nrow=4,ncol=5)
colnames(average_distance_per_realization_network_size) <- c("network_size_quantile_number","mean_network_size","SD_network_size","mean_distance","number_of_networks")

#obtaining all food webs sourced from a publication that produced only a single network
realization <- metadata[metadata$author=="One_network_per_publication",]

realization$number_of_nodes <- as.numeric(realization$number_of_nodes)

#identifying the quartile that each food web belongs to (based on number of nodes)
Quartile<-cut(realization$number_of_nodes,quantile(realization$number_of_nodes,type=8),include.lowest=TRUE,labels=FALSE)

#adding quartile information as its own column
realization <- cbind(realization,Quartile)

#obtaining only quartile and food web name
network_name_and_quartile <- as.data.frame(cbind(realization$name,realization$Quartile))

#naming the columns
colnames(network_name_and_quartile) <- c("name","network_size_quartile")

network_name_and_quartile$network_size_quartile <- as.numeric(network_name_and_quartile$network_size_quartile)

all_networks_and_quartiles <- merge(metadata,network_name_and_quartile,by="name",all=T)

#looping through the pairwise GCD-11 matrix to find food webs associated 
#with specific author/publication (here "4" is the number of quartiles)
for (k in 1:4) {
  
  #list of rows and columns to delete which are not of the 
  #author/publication grouping of interest right now
  rows_columns_delete <- c()
  
  #looping through the GCD-11 matrix
  for (j in 1:ncol(ecological)) {
    
    #the name of the row we are looking through
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(all_networks_and_quartiles[rowNumber,])
      
      #if the row name in the matrix is the food web we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not an author's food web we want, or is NA, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the GCD-11 matrix)
        if (is.na(rowOfInterest$network_size_quartile)) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        } else if (rowOfInterest$network_size_quartile != k) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
    }
    
  }
  
  #deleting food webs that are not the author/publication grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #quartile number
  average_distance_per_realization_network_size[k,1] <- k
  
  #mean food web size between specific author/publication's food webs
  average_distance_per_realization_network_size[k,2] <- round(mean(realization[realization$Quartile==k,]$number_of_nodes),2)
  
  #standard deviation in food web size between specific author/publication's food webs
  average_distance_per_realization_network_size[k,3] <- round(sd(realization[realization$Quartile==k,]$number_of_nodes),2)
  
  #mean pairwise GCD-11 between specific author/publication's food webs
  average_distance_per_realization_network_size[k,4] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific author/publication grouping of interest
  average_distance_per_realization_network_size[k,5] <- ncol(ecological_intermediate)
  
}

##
##  All information for table S5
##

average_distance_per_realization_network_size
