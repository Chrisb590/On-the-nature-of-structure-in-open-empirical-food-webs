########################################################
##
##
##
##  For recreating table 2 of Brimacombe et al. (2024):
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

##
##  Mean pairwise GCD-11 between (i) food webs sourced from publications
##  that each produced a single network and (ii) food webs sourced from
##  publications that produced multiple networks
##

#storage for mean pairwise GCD-11 between food webs sourced
#from the same publication grouping
each_average_distance <- matrix(,nrow=(length(levels(metadata$author))),ncol=3)
colnames(each_average_distance) <- c("network_number","mean_distance","number_of_networks")

#the different publication groupings
for (i in 1:nrow(each_average_distance)) {
  
  each_average_distance[i,1] <- levels(metadata$author)[i]
  
}

#looping through the pairwise GCD-11 matrix to find 
#food webs associated with specific publication groupings
for (i in 1:nrow(each_average_distance)) {
  
  #list of rows and columns to delete which are not of the 
  #publication grouping of interest right now
  rows_columns_delete <- c()
  
  for (j in 1:ncol(ecological)) {
    
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata[rowNumber,])
      
      #if the row name in the matrix is the food web we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not a publication's food web we want, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$author != each_average_distance[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
    }
    
  }
  
  #deleting food webs that are not sourced from the publication grouping of interest
  #(pairwise GCD-11 matrix is symmetric so we can delete the same rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #mean pairwise GCD-11 between a specific publication's food webs
  each_average_distance[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific publication
  each_average_distance[i,3] <- ncol(ecological_intermediate)
  
}

each_average_distance <- as.data.frame(each_average_distance)
each_average_distance$mean_distance <- as.numeric(each_average_distance$mean_distance)
each_average_distance$number_of_networks <- as.numeric(each_average_distance$number_of_networks)

##
##  binning food webs sourced from the same publication that produced 
##  multiple networks by decade
##

#storage for mean pairwise GCD-11 between food webs sourced from the same publication that produced multiple networks
each_average_distance_for_publication <- matrix(,nrow=(length(levels(metadata$author))),ncol=6)
colnames(each_average_distance_for_publication) <- c("publication_name","publication_year","mean_distance","SD_distance","number_of_networks","Identity")

#obtaining information about all publications that provided multiple networks
for (i in 1:nrow(each_average_distance_for_publication)) {
  
  each_average_distance_for_publication[i,1] <- levels(metadata$author)[i]
  each_average_distance_for_publication[i,2] <- metadata[which(metadata$author==levels(metadata$author)[i])[1],]$year_published
  each_average_distance_for_publication[i,6] <- "Publication"
  
}

each_average_distance_for_publication <- as.data.frame(each_average_distance_for_publication)
each_average_distance_for_publication <- each_average_distance_for_publication[!(each_average_distance_for_publication$publication_name=="One_network_per_publication"),]

#looping through the pairwise GCD-11 matrix to find food webs 
#associated with specific publications that provided multiple
#networks
for (i in 1:nrow(each_average_distance_for_publication)) {
  
  #list of rows and columns to delete which are not of the 
  #publication grouping of interest right now
  rows_columns_delete <- c()
  
  for (j in 1:ncol(ecological)) {
    
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata[rowNumber,])
      
      #if the row name in the matrix is the food web we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not the food web we want, we will need 
        #to delete it (here, we just keep the location of this 
        #food web in the pairwise GCD-11 matrix)
        if (rowOfInterest$author != each_average_distance_for_publication[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
    }
    
  }
  
  #deleting food webs that are not the publication grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #mean pairwise GCD-11 between specific publication's food webs
  each_average_distance_for_publication[i,3] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #standard deviation in the pairwise GCD-11 between specific publication's food webs
  each_average_distance_for_publication[i,4] <- round(sd(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific publication grouping of interest
  each_average_distance_for_publication[i,5] <- ncol(ecological_intermediate)
  
}

each_average_distance_for_publication <- as.data.frame(each_average_distance_for_publication)
each_average_distance_for_publication$mean_distance <- as.numeric(each_average_distance_for_publication$mean_distance)
each_average_distance_for_publication$publication_year <- as.numeric(each_average_distance_for_publication$publication_year)
each_average_distance_for_publication$number_of_networks <- as.numeric(each_average_distance_for_publication$number_of_networks)
each_average_distance_for_publication$SD_distance <- as.numeric(each_average_distance_for_publication$SD_distance)

##
##  All information for table 2
##

#mean pairwise GCD-11 of food webs sourced from publications that only produced a single network
webs_from_pubs_with_single_network <- each_average_distance[each_average_distance$network_number=="One_network_per_publication",]
webs_from_pubs_with_single_network

#weighted mean pairwise GCD-11 of food webs sourced from the same publication that produced multiple networks
webs_from_pubs_with_multiple_networks <- each_average_distance[!(each_average_distance$network_number=="One_network_per_publication"),]
round(weighted.mean(webs_from_pubs_with_multiple_networks$mean_distance,webs_from_pubs_with_multiple_networks$number_of_networks),2)

#weighted mean of mean pairwise GCD-11 from publications before or during 1990s
each_average_distance_for_publication_before_1990 <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<2000,]
round(weighted.mean(each_average_distance_for_publication_before_1990$mean_distance,each_average_distance_for_publication_before_1990$number_of_networks),2)

#weighted mean of mean pairwise GCD-11 from publications after 1990s
each_average_distance_for_publication_after_1990 <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year>=2000,]
round(weighted.mean(each_average_distance_for_publication_after_1990$mean_distance,each_average_distance_for_publication_after_1990$number_of_networks),2)
