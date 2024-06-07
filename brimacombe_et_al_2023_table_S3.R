########################################################
##
##
##
##  For recreating table S3 of Brimacombe et al. (2024):
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
##  Used to evaluate the median pairwise GCD-11 between
##  food webs from the same ecological system
##

#storage for median pairwise GCD-11 between food webs from the same ecosystem
each_median_distance_same_ecosystem <- matrix(,nrow=(length(levels(as.factor(metadata$Primary_type)))),ncol=3)
colnames(each_median_distance_same_ecosystem) <- c("Ecosystem_type","median_distance","number_of_networks")

for (i in 1:nrow(each_median_distance_same_ecosystem)) {
  
  each_median_distance_same_ecosystem[i,1] <- levels(as.factor(metadata$Primary_type))[i]
  
}

#looping through the pairwise GCD-11 matrix to find food webs associated 
#with specific ecosystems
for (i in 1:nrow(each_median_distance_same_ecosystem)) {
  
  #list of rows and columns to delete which are not of the 
  #ecosystem grouping of interest right now
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
        
        #if the row is not an ecosystem's food web we want, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$Primary_type != each_median_distance_same_ecosystem[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
      #if food web not in metadata, we do not want it
      if (rowNumber > nrow(metadata)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
    }
    
  }
  
  #deleting food webs that are not the ecosystem grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #median pairwise GCD-11 between specific ecosystem's food webs
  each_median_distance_same_ecosystem[i,2] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific ecosystem grouping of interest
  each_median_distance_same_ecosystem[i,3] <- ncol(ecological_intermediate)
  
}

each_median_distance_same_ecosystem <- as.data.frame(each_median_distance_same_ecosystem)
each_median_distance_same_ecosystem$median_distance <- as.numeric(each_median_distance_same_ecosystem$median_distance)
each_median_distance_same_ecosystem$number_of_networks <- as.numeric(each_median_distance_same_ecosystem$number_of_networks)

##
##  Median pairwise GCD-11 between "terrestrial" food webs with 
##  "Digel et al. (2014)" food webs removed
##

#storage for median pairwise GCD-11 between food webs from the same ecosystem
each_median_distance_terrestrial_no_digel <- matrix(,nrow=1,ncol=3)
colnames(each_median_distance_terrestrial_no_digel) <- c("Ecosystem_type","mean_distance","number_of_networks")
each_median_distance_terrestrial_no_digel[1,1] <- "Terrestrial"

#looping through the pairwise GCD-11 matrix to find food webs associated 
#with specific ecosystems
for (i in 1:nrow(each_median_distance_terrestrial_no_digel)) {
  
  #list of rows and columns to delete which are not of the 
  #ecosystem grouping of interest right now
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
        
        #if the row is not an ecosystem's food web we want, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$Primary_type != each_median_distance_terrestrial_no_digel[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        #if the row is Digel_et_al_2014 remove
        if (rowOfInterest$author == "Digel_et_al_2014") {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
      #if food web not in metadata, we do not want it
      if (rowNumber > nrow(metadata)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
    }
    
  }
  
  #deleting food webs that are not the ecosystem grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #median pairwise GCD-11 between specific ecosystem's food webs
  each_median_distance_terrestrial_no_digel[i,2] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific ecosystem grouping of interest
  each_median_distance_terrestrial_no_digel[i,3] <- ncol(ecological_intermediate)
  
}

##
##  Used to evaluate the median pairwise GCD-11 between
##  food webs from different ecological systems
##

##
##  Median pairwise GCD-11 between "aquatic" and 
##  "aquatic and terrestrial" food webs
##

#list of rows and columns to delete which are not of the 
#type we are looking for (i.e., "terrestrial")
rows_columns_delete <- c()

#columns to delete which are the food web type we are 
#trying to compare too (i.e., "aquatic") so as not to
#compare food webs twice (i.e, they are in both the rows
#and columns)
columns_delete <- c()

#rows to delete which are the food web type we are 
#trying to compare too (i.e., "aquatic and terrestrial") 
#so as not to compare food webs twice (i.e, they are in 
#both the rows and columns)
rows_delete <- c()

for (j in 1:ncol(ecological)) {
  
  network_name <- row.names(ecological)[j]
  
  found = FALSE
  rowNumber = 1
  
  #looping to find what number the food web is
  while (found != TRUE) {
    
    rowOfInterest <- as.data.frame(metadata[rowNumber,])
    
    #if the row name in the matrix is the food web we are looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is a "terrestrial" food web, we will need 
      #to delete it
      if (rowOfInterest$Primary_type == "Terrestrial") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      #if the row is an "aquatic" food web, we will need 
      #to delete it from the column   
      } else if (rowOfInterest$Primary_type == "Aquatic") {
        
        columns_delete <- c(columns_delete,j)
        
      #if the row is an "aquatic and terrestrial" food web, we will need 
      #to delete it from the row         
      } else if (rowOfInterest$Primary_type == "Aquatic and terrestrial") {
        
        rows_delete <- c(rows_delete,j)
        
      }
      
      found = TRUE
      
    #otherwise we have not found the food web in the metadata list
    } else {
      
      rowNumber = rowNumber + 1
      
    }
    
    #if food web not in metadata, we do not want it
    if (rowNumber > nrow(metadata)) {
      
      rows_columns_delete <- c(rows_columns_delete,j)
      
      found = TRUE
      
    }
    
  }
  
}

#deleting food webs that are not of interest
each_median_distance_aquatic_v_aquatic_terrestrial <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]

##
##  Median pairwise GCD-11 between "aquatic" and 
##  "terrestrial" food webs
##

#list of rows and columns to delete which are not of the 
#type we are looking for (i.e., "aquatic and terrestrial")
rows_columns_delete <- c()

#columns to delete which are the food web type we are 
#trying to compare too (i.e., "aquatic") so as not to
#compare food webs twice (i.e, they are in both the rows
#and columns)
columns_delete <- c()

#rows to delete which are the food web type we are 
#trying to compare too (i.e., "terrestrial") so as 
#not to compare food webs twice (i.e, they are in 
#both the rows and columns)
rows_delete <- c()

for (j in 1:ncol(ecological)) {
  
  network_name <- row.names(ecological)[j]
  
  found = FALSE
  rowNumber = 1
  
  #looping to find what number the food web is
  while (found != TRUE) {
    
    rowOfInterest <- as.data.frame(metadata[rowNumber,])
    
    #if the row name in the matrix is the food web we are looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is an "aquatic and terrestrial" food web, we will need 
      #to delete it
      if (rowOfInterest$Primary_type == "Aquatic and terrestrial") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      #if the row is an "aquatic" food web, we will need 
      #to delete it from the column
      } else if (rowOfInterest$Primary_type == "Aquatic") {
        
        columns_delete <- c(columns_delete,j)
        
      #if the row is a "terrestrial" food web, we will need 
      #to delete it from the row 
      } else if (rowOfInterest$Primary_type == "Terrestrial") {
        
        rows_delete <- c(rows_delete,j)
        
      }
      
      found = TRUE
      
    #otherwise we have not found the food web in the metadata list
    } else {
      
      rowNumber = rowNumber + 1
      
    }
    
    #if food web not in metadata, we do not want it
    if (rowNumber > nrow(metadata)) {
      
      rows_columns_delete <- c(rows_columns_delete,j)
      
      found = TRUE
      
    }
    
  }
  
}

#deleting food webs that are not of interest
each_median_distance_aquatic_v_terrestrial <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]

##
##  Median pairwise GCD-11 between "aquatic and terrestrial" 
##  and "terrestrial" food webs
##

#list of rows and columns to delete which are not of the 
#type we are looking for (i.e., "aquatic")
rows_columns_delete <- c()

#columns to delete which are the food web type we are 
#trying to compare too (i.e., "aquatic and terrestrial") 
#so as not to compare food webs twice (i.e, they are in 
#both the rows and columns)
columns_delete <- c()

#rows to delete which are the food web type we are 
#trying to compare too (i.e., "terrestrial") so as 
#not to compare food webs twice (i.e, they are in 
#both the rows and columns)
rows_delete <- c()

for (j in 1:ncol(ecological)) {
  
  network_name <- row.names(ecological)[j]
  
  found = FALSE
  rowNumber = 1
  
  #looping to find what number the food web is
  while (found != TRUE) {
    
    rowOfInterest <- as.data.frame(metadata[rowNumber,])
    
    #if the row name in the matrix is the food web we are looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is an "aquatic" food web, we will need 
      #to delete it
      if (rowOfInterest$Primary_type == "Aquatic") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      #if the row is an "aquatic and terrestrial" food web, we will need 
      #to delete it from the column 
      } else if (rowOfInterest$Primary_type == "Aquatic and terrestrial") {
        
        columns_delete <- c(columns_delete,j)
        
      #if the row is a "terrestrial" food web, we will need 
      #to delete it from the row       
      } else if (rowOfInterest$Primary_type == "Terrestrial") {
        
        rows_delete <- c(rows_delete,j)
        
      }
      
      found = TRUE
      
    #otherwise we have not found the food web in the metadata list
    } else {
      
      rowNumber = rowNumber + 1
      
    }
    
    #if food web not in metadata, we do not want it
    if (rowNumber > nrow(metadata)) {
      
      rows_columns_delete <- c(rows_columns_delete,j)
      
      found = TRUE
      
    }
    
  }
  
}

#deleting food webs that are not of interest
each_median_distance_aquatic_terrestrial_v_terrestrial <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]

##
##  All information for table S3
##

#median pairwise GCD-11 between food webs from the same ecosystem
each_median_distance_same_ecosystem

#median pairwise GCD-11 between "terrestrial" food webs with Digel et al. (2014) removed
each_median_distance_terrestrial_no_digel

#median pairwise GCD-11 between "aquatic" and "aquatic and terrestrial" food webs
round(median(each_median_distance_aquatic_v_aquatic_terrestrial),2)

#median pairwise GCD-11 between "aquatic" and "terrestrial" food webs
round(median(each_median_distance_aquatic_v_terrestrial),2)

#median pairwise GCD-11 between "aquatic and terrestrial" and "terrestrial" food webs
round(median(each_median_distance_aquatic_terrestrial_v_terrestrial),2)
