########################################################
##
##
##
##  For recreating figure S10 of Brimacombe et al. (2024):
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

#storage for pairwise GCD-11
distances<- matrix(,nrow=37401,ncol=3)
colnames(distances) <- c("absolute_network_size","SD_network_size","pairwise_distance")

#obtaining all food webs
metadata_intermediate <- metadata

metadata_intermediate$number_of_nodes <- as.numeric(metadata_intermediate$number_of_nodes)

counter <- 1

#looping through the rows of the GCD-11 matrix
for (k in 2:nrow(ecological)) {
  
  #the name of the row in the GCD-11 matrix
  row_name <- rownames(ecological)[k]
  
  found_row = FALSE
  rowNumber = 1
  
  #looping to find the number of nodes that the row food web has
  while (found_row != TRUE) {
    
    rowOfInterest <- metadata_intermediate[rowNumber,]
    
    #if the row name in the matrix is the food web we are looking for
    if (as.character(rowOfInterest$name) == row_name) {
      
      row_size <- rowOfInterest$number_of_nodes
      found_row = TRUE
      
      
    #otherwise we have not found the food web in the metadata list
    } else {
      
      rowNumber = rowNumber + 1
      
    }
    
  }
  
  #looping through the columns of the GCD-11 matrix
  for (i in 1:(k-1)){
    
    #the name of the column in the GCD-11 matrix
    colmn_name <- colnames(ecological)[i]
    
    found_column = FALSE
    columnNumber = 1
    
    #looping to find the number of nodes that the column food web has
    while (found_column != TRUE) {
      
      columnOfInterest <- metadata_intermediate[columnNumber,]
      
      #if the column name in the matrix is the food web we are looking for
      if (as.character(columnOfInterest$name) == colmn_name) {
        
        column_size <- columnOfInterest$number_of_nodes
        found_column = TRUE
        
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        columnNumber = columnNumber + 1
        
      }
      
    }
    
    distances[counter,1] <- abs(row_size-column_size)
    distances[counter,2] <- sd(c(row_size,column_size))
    distances[counter,3] <- as.numeric(ecological[k,i])
    counter <- counter + 1
    
  }
  
}

distances <- as.data.frame(distances)
colnames(distances) <- c("absolute_network_size","SD_network_size","pairwise_distance")

#
# Figure S10
# 

#absolute food web size and pairwise GCD-11s
summary(lm(distances$absolute_network_size~distances$pairwise_distance))

ggplot(data=distances,aes(x=absolute_network_size,y=pairwise_distance)) +
  geom_point(alpha = 0.5, color=c("#011288"), size=4) +
  coord_cartesian(ylim = c(0, 7.5)) +
  theme_bw() +
  xlab("Absolute difference in network size") +
  ylab("Pairwise GCD-11") +
  theme(axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=20, colour ="black"),
        axis.title.y = element_text(size=20, colour ="black"),
        text = element_text(size = 20,color="black"),
        axis.text.x = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color="black")) +
  geom_smooth(method='lm', formula= y~x)

#standard deviation of food web size and pairwise GCD-11s
summary(lm(distances$SD_network_size~distances$pairwise_distance))

ggplot(data=distances,aes(x=SD_network_size,y=pairwise_distance)) +
  geom_point(alpha = 0.5, color=c("#011288"), size=4) +
  coord_cartesian(ylim = c(0, 7.5)) +
  theme_bw() +
  xlab("Standard deviation in network size") +
  ylab("Pairwise GCD-11") +
  theme(axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=20, colour ="black"),
        axis.title.y = element_text(size=20, colour ="black"),
        text = element_text(size = 20,color="black"),
        axis.text.x = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color="black")) +
  geom_smooth(method='lm', formula= y~x)
