########################################################
##
##
##
##  For recreating figure S9 of Brimacombe et al. (2024):
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

#storage for the mean number of nodes in networks from the same publication
#and the standard deviation in the number of nodes from the same publication
storage_for_mean_and_sd_network_size <- matrix(,nrow(each_average_distance),3)
colnames(storage_for_mean_and_sd_network_size) <- c("network_number","mean_number_of_nodes","sd_number_of_nodes")

#obtaining the mean number of nodes in networks from the same publication
#and the standard deviation in the number of nodes from the same publication
for (i in 1:nrow(each_average_distance)) {
  
  storage_for_mean_and_sd_network_size[i,1] <- levels(metadata$author)[i]
  intermediate <- metadata[metadata$author==levels(metadata$author)[i],]
  storage_for_mean_and_sd_network_size[i,2] <- round(mean(intermediate$number_of_nodes),2)
  storage_for_mean_and_sd_network_size[i,3] <- round(sd(intermediate$number_of_nodes),2)
  
}

#merging information about the mean pairwise GCD-11 between networks from the same publication
#and their respective mean number of nodes and standard deviation in the number of nodes
for_regression <- merge(each_average_distance,storage_for_mean_and_sd_network_size,by="network_number")
for_regression <- for_regression[!(for_regression$network_number=="One_network_per_publication"),]
for_regression$mean_distance <- as.numeric(for_regression$mean_distance)
for_regression$mean_number_of_nodes <- as.numeric(for_regression$mean_number_of_nodes)
for_regression$sd_number_of_nodes <- as.numeric(for_regression$sd_number_of_nodes)
for_regression <- droplevels(for_regression)

#
# Figure S9
# 

summary(lm(for_regression$mean_distance~for_regression$mean_number_of_nodes))

ggplot(data=for_regression,aes(x=mean_number_of_nodes,y=mean_distance)) +
  geom_point(alpha = 0.5, color=c("#011288"), size=4) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw() +
  xlab("Mean number of nodes per publication") +
  ylab("Mean pairwise GCD-11") +
  theme(axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=20, colour ="black"),
        axis.title.y = element_text(size=20, colour ="black"),
        text = element_text(size = 20,color="black"),
        axis.text.x = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color="black")) +
  geom_smooth(method='lm', formula= y~x)

summary(lm(for_regression$mean_distance~for_regression$sd_number_of_nodes))

ggplot(data=for_regression,aes(x=sd_number_of_nodes,y=mean_distance)) +
  geom_point(alpha = 0.5, color=c("#011288"), size=4) +
  coord_cartesian(ylim = c(0, 3.5)) +
  theme_bw() +
  xlab("Standard deviation of the number of nodes per publication") +
  ylab("Mean pairwise GCD-11") +
  theme(axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=20, colour ="black"),
        axis.title.y = element_text(size=20, colour ="black"),
        text = element_text(size = 20,color="black"),
        axis.text.x = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color="black")) +
  geom_smooth(method='lm', formula= y~x)
