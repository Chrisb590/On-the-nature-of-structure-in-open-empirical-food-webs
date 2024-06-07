########################################################
##
##
##
##  For recreating table S2 of Brimacombe et al. (2024):
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

#number of times you want to randomly sample/number of realizations
sample_times <- 20

########################################################

##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs that are "lake", "marine", "river", and 
##  "stream"
##

##
##  Mean pairwise GCD-11 between "lake" food webs 
##  with publication removed
##

#Angelini_et_al_2013
Angelini_et_al_2013_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Angelini_et_al_2013_random_samples) <- c(colnames(metadata))
Angelini_et_al_2013 <-  metadata[metadata$author=="Angelini_et_al_2013",]

for (i in 1:sample_times) {
  
  Angelini_et_al_2013_random_samples[i,] <- as.matrix(sample_n(Angelini_et_al_2013, 1))
  
}

#Stewart_and_Sprules_2011
Stewart_and_Sprules_2011_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Stewart_and_Sprules_2011_random_samples) <- c(colnames(metadata))
Stewart_and_Sprules_2011 <-  metadata[metadata$author=="Stewart_and_Sprules_2011",]

for (i in 1:sample_times) {
  
  Stewart_and_Sprules_2011_random_samples[i,] <- as.matrix(sample_n(Stewart_and_Sprules_2011, 1))
  
}

#Alcorlo_et_al_2001
Alcorlo_et_al_2001_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Alcorlo_et_al_2001_random_samples) <- c(colnames(metadata))
Alcorlo_et_al_2001 <-  metadata[metadata$author=="Alcorlo_et_al_2001",]

for (i in 1:sample_times) {
  
  Alcorlo_et_al_2001_random_samples[i,] <- as.matrix(sample_n(Alcorlo_et_al_2001, 1))
  
}

#Cohen_et_al_2003
Cohen_et_al_2003_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Cohen_et_al_2003_random_samples) <- c(colnames(metadata))
Cohen_et_al_2003 <-  metadata[metadata$author=="Cohen_et_al_2003",]

for (i in 1:sample_times) {
  
  Cohen_et_al_2003_random_samples[i,] <- as.matrix(sample_n(Cohen_et_al_2003, 1))
  
}

#Fryer_1959
Fryer_1959_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Fryer_1959_random_samples) <- c(colnames(metadata))
Fryer_1959 <-  metadata[metadata$author=="Fryer_1959",]

for (i in 1:sample_times) {
  
  Fryer_1959_random_samples[i,] <- as.matrix(sample_n(Fryer_1959, 1))
  
}

#Havens_1992
Havens_1992_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Havens_1992_random_samples) <- c(colnames(metadata))
Havens_1992 <-  metadata[metadata$author=="Havens_1992",]
Havens_1992 <-  Havens_1992[Havens_1992$Primary_type!="Aquatic and terrestrial",]

for (i in 1:sample_times) {
  
  Havens_1992_random_samples[i,] <- as.matrix(sample_n(Havens_1992, 1))
  
}

#for food webs from a publication that produced only a single network and are "lake"
single_aquatic_Lake_food_webs <- metadata[metadata$author=="One_network_per_publication",]
single_aquatic_Lake_food_webs <- single_aquatic_Lake_food_webs[single_aquatic_Lake_food_webs$ Specialized_aquatic_type=="Lake",]

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_lake <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_lake) <- c("realization","mean_distance","number_of_networks")

for (i in 1:sample_times){
  
  #obtaining each food web from publications that were randomly sampled for a specific realization
  realization <- as.data.frame(rbind(single_aquatic_Lake_food_webs,Angelini_et_al_2013_random_samples[i,],Stewart_and_Sprules_2011_random_samples[i,],
                                     Alcorlo_et_al_2001_random_samples[i,],Cohen_et_al_2003_random_samples[i,],Fryer_1959_random_samples[i,],
                                     Havens_1992_random_samples[i,]))
  
  #list of rows and columns to delete which are not
  #of interest right now
  rows_columns_delete <- c()
  
  #looping through the GCD-11 matrix
  for (j in 1:ncol(ecological)) {
    
    #the name of the row we are looking through
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      #if we exceed the number of rows in the realization, this food web should be removed
      if (rowNumber>nrow(realization)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
      #if we have not exceed the number of rows in the realization
      if (rowNumber<=nrow(realization)) {
        
        rowOfInterest <- as.data.frame(realization[rowNumber,])
        
        #if the row name in the matrix is the food web we are looking for
        #we want to keep this food web
        if (as.character(rowOfInterest$name) == network_name) {
          
          found = TRUE
          
        #otherwise we have not found the food web in the metadata list
        } else {
          
          rowNumber = rowNumber + 1
          
        }
        
      }
      
      
    }
    
  }
  
  #deleting food webs that are not the grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #realization number
  average_distance_lake[i,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_lake[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the realization
  average_distance_lake[i,3] <- ncol(ecological_intermediate)
  
}


##
##  Mean pairwise GCD-11 between "marine" food webs 
##  with publication removed
##

#Baeta_et_al_2011
Baeta_et_al_2011_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Baeta_et_al_2011_random_samples) <- c(colnames(metadata))
Baeta_et_al_2011 <-  metadata[metadata$author=="Baeta_et_al_2011",]

for (i in 1:sample_times) {
  
  Baeta_et_al_2011_random_samples[i,] <- as.matrix(sample_n(Baeta_et_al_2011, 1))
  
}

#Menge_and_Sutherland_1976
Menge_and_Sutherland_1976_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Menge_and_Sutherland_1976_random_samples) <- c(colnames(metadata))
Menge_and_Sutherland_1976 <-  metadata[metadata$author=="Menge_and_Sutherland_1976",]

for (i in 1:sample_times) {
  
  Menge_and_Sutherland_1976_random_samples[i,] <- as.matrix(sample_n(Menge_and_Sutherland_1976, 1))
  
}

#for food webs from publications that produced only a single network and are "marine"
single_aquatic_Marine_food_webs <- metadata[metadata$author=="One_network_per_publication",]
single_aquatic_Marine_food_webs <- single_aquatic_Marine_food_webs[single_aquatic_Marine_food_webs$ Specialized_aquatic_type=="Marine",]

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_marine <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_marine) <- c("realization","mean_distance","number_of_networks")

for (i in 1:sample_times){
  
  #obtaining each food web from publications that were randomly sampled for a specific realization
  realization <- as.data.frame(rbind(single_aquatic_Marine_food_webs,Baeta_et_al_2011_random_samples[i,],
                                     Menge_and_Sutherland_1976_random_samples[i,]))
  
  #list of rows and columns to delete which are not
  #of interest right now
  rows_columns_delete <- c()
  
  #looping through the GCD-11 matrix
  for (j in 1:ncol(ecological)) {
    
    #the name of the row we are looking through
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      #if we exceed the number of rows in the realization, this food web should be removed
      if (rowNumber>nrow(realization)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
      #if we have not exceed the number of rows in the realization
      if (rowNumber<=nrow(realization)) {
        
        rowOfInterest <- as.data.frame(realization[rowNumber,])
        
        #if the row name in the matrix is the food web we are looking for
        #we want to keep this food web
        if (as.character(rowOfInterest$name) == network_name) {
          
          found = TRUE
          
        #otherwise we have not found the food web in the metadata list
        } else {
          
          rowNumber = rowNumber + 1
          
        }
        
      }
      
      
    }
    
  }
  
  #deleting food webs that are not the grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #realization number
  average_distance_marine[i,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_marine[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the realization
  average_distance_marine[i,3] <- ncol(ecological_intermediate)
  
}

##
##  Mean pairwise GCD-11 between "river" food webs
##  with publication removed
##

#Angelini_et_al_2006
Angelini_et_al_2006_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Angelini_et_al_2006_random_samples) <- c(colnames(metadata))
Angelini_et_al_2006 <-  metadata[metadata$author=="Angelini_et_al_2006",]

for (i in 1:sample_times) {
  
  Angelini_et_al_2006_random_samples[i,] <- as.matrix(sample_n(Angelini_et_al_2006, 1))
  
}

#for food webs from publications that produced only a single network and are "river"
single_aquatic_River_food_webs <- metadata[metadata$author=="One_network_per_publication",]
single_aquatic_River_food_webs <- single_aquatic_River_food_webs[single_aquatic_River_food_webs$Specialized_aquatic_type=="River",]

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_river <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_river) <- c("realization","mean_distance","number_of_networks")

for (i in 1:sample_times){
  
  #obtaining each food web from publications that were randomly sampled for a specific realization
  realization <- as.data.frame(rbind(single_aquatic_River_food_webs,Angelini_et_al_2006_random_samples[i,]))
  
  #list of rows and columns to delete which are not
  #of interest right now
  rows_columns_delete <- c()
  
  #looping through the GCD-11 matrix
  for (j in 1:ncol(ecological)) {
    
    #the name of the row we are looking through
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      #if we exceed the number of rows in the realization, this food web should be removed
      if (rowNumber>nrow(realization)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
      #if we have not exceed the number of rows in the realization
      if (rowNumber<=nrow(realization)) {
        
        rowOfInterest <- as.data.frame(realization[rowNumber,])
        
        #if the row name in the matrix is the food web we are looking for
        #we want to keep this food web
        if (as.character(rowOfInterest$name) == network_name) {
          
          found = TRUE
          
        #otherwise we have not found the food web in the metadata list
        } else {
          
          rowNumber = rowNumber + 1
          
        }
        
      }
      
      
    }
    
  }
  
  #deleting food webs that are not the grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #realization number
  average_distance_river[i,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_river[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the realization
  average_distance_river[i,3] <- ncol(ecological_intermediate)
  
}

##
##  Mean pairwise GCD-11 between "stream" food webs
##  with publication removed
##

#Closs_and_Lake_1994
Closs_and_Lake_1994_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Closs_and_Lake_1994_random_samples) <- c(colnames(metadata))
Closs_and_Lake_1994 <-  metadata[metadata$author=="Closs_and_Lake_1994",]

for (i in 1:sample_times) {
  
  Closs_and_Lake_1994_random_samples[i,] <- as.matrix(sample_n(Closs_and_Lake_1994, 1))
  
}

#Layer_et_al_2010
Layer_et_al_2010_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Layer_et_al_2010_random_samples) <- c(colnames(metadata))
Layer_et_al_2010 <-  metadata[metadata$author=="Layer_et_al_2010",]

for (i in 1:sample_times) {
  
  Layer_et_al_2010_random_samples[i,] <- as.matrix(sample_n(Layer_et_al_2010, 1))
  
}

#O_Gorman_et_al_2019
O_Gorman_et_al_2019_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(O_Gorman_et_al_2019_random_samples) <- c(colnames(metadata))
O_Gorman_et_al_2019 <-  metadata[metadata$author=="O_Gorman_et_al_2019",]

for (i in 1:sample_times) {
  
  O_Gorman_et_al_2019_random_samples[i,] <- as.matrix(sample_n(O_Gorman_et_al_2019, 1))
  
}

#Parker_and_Huryn_2006
Parker_and_Huryn_2006_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Parker_and_Huryn_2006_random_samples) <- c(colnames(metadata))
Parker_and_Huryn_2006 <-  metadata[metadata$author=="Parker_and_Huryn_2006",]

for (i in 1:sample_times) {
  
  Parker_and_Huryn_2006_random_samples[i,] <- as.matrix(sample_n(Parker_and_Huryn_2006, 1))
  
}

#Tavares-Cromar_and_Williams_1996
Tavares_Cromar_and_Williams_1996_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Tavares_Cromar_and_Williams_1996_random_samples) <- c(colnames(metadata))
Tavares_Cromar_and_Williams_1996 <-  metadata[metadata$author=="Tavares-Cromar_and_Williams_1996",]

for (i in 1:sample_times) {
  
  Tavares_Cromar_and_Williams_1996_random_samples[i,] <- as.matrix(sample_n(Tavares_Cromar_and_Williams_1996, 1))
  
}

#Thompson_and_Townsend_2003
Thompson_and_Townsend_2003_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Thompson_and_Townsend_2003_random_samples) <- c(colnames(metadata))
Thompson_and_Townsend_2003 <-  metadata[metadata$author=="Thompson_and_Townsend_2003",]

for (i in 1:sample_times) {
  
  Thompson_and_Townsend_2003_random_samples[i,] <- as.matrix(sample_n(Thompson_and_Townsend_2003, 1))
  
}

#Thompson_and_Townsend_2004
Thompson_and_Townsend_2004_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Thompson_and_Townsend_2004_random_samples) <- c(colnames(metadata))
Thompson_and_Townsend_2004 <-  metadata[metadata$author=="Thompson_and_Townsend_2004",]

for (i in 1:sample_times) {
  
  Thompson_and_Townsend_2004_random_samples[i,] <- as.matrix(sample_n(Thompson_and_Townsend_2004, 1))
  
}

#for food webs from publications that produced only a single network and are "stream"
single_aquatic_Stream_food_webs <- metadata[metadata$author=="One_network_per_publication",]
single_aquatic_Stream_food_webs <- single_aquatic_Stream_food_webs[single_aquatic_Stream_food_webs$Specialized_aquatic_type=="Stream",]

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_stream <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_stream) <- c("realization","mean_distance","number_of_networks")

for (i in 1:sample_times){
  
  #obtaining each food web from publications that were randomly sampled for a specific realization
  realization <- as.data.frame(rbind(single_aquatic_Stream_food_webs,Closs_and_Lake_1994_random_samples[i,],
                                     Layer_et_al_2010_random_samples[i,],O_Gorman_et_al_2019_random_samples[i,],
                                     Parker_and_Huryn_2006_random_samples[i,],Tavares_Cromar_and_Williams_1996_random_samples[i,],
                                     Thompson_and_Townsend_2003_random_samples[i,],Thompson_and_Townsend_2004_random_samples[i,]))
  
  #list of rows and columns to delete which are not
  #of interest right now
  rows_columns_delete <- c()
  
  #looping through the GCD-11 matrix
  for (j in 1:ncol(ecological)) {
    
    #the name of the row we are looking through
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      #if we exceed the number of rows in the realization, this food web should be removed
      if (rowNumber>nrow(realization)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
      #if we have not exceed the number of rows in the realization
      if (rowNumber<=nrow(realization)) {
        
        rowOfInterest <- as.data.frame(realization[rowNumber,])
        
        #if the row name in the matrix is the food web we are looking for
        #we want to keep this food web
        if (as.character(rowOfInterest$name) == network_name) {
          
          found = TRUE
          
        #otherwise we have not found the food web in the metadata list
        } else {
          
          rowNumber = rowNumber + 1
          
        }
        
      }
      
      
    }
    
  }
  
  #deleting food webs that are not the grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #realization number
  average_distance_stream[i,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_stream[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the realization
  average_distance_stream[i,3] <- ncol(ecological_intermediate)
  
}

##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs from different ecological systems 
##  with publication removed
##

##
##  Mean pairwise GCD-11 between "lake" 
##  and "marine" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_lake_v_marine <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_lake_v_marine) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

for (i in 1:sample_times){
  
  #obtaining each "lake" food web that was randomly sampled for a specific realization
  realization_lake <- as.data.frame(rbind(single_aquatic_Lake_food_webs,Angelini_et_al_2013_random_samples[i,],
                                          Stewart_and_Sprules_2011_random_samples[i,],Alcorlo_et_al_2001_random_samples[i,],
                                          Cohen_et_al_2003_random_samples[i,],Fryer_1959_random_samples[i,],
                                          Havens_1992_random_samples[i,]))
  
  #obtaining each "marine" food web that was randomly sampled for a specific realization
  realization_marine <- as.data.frame(rbind(single_aquatic_Marine_food_webs,Baeta_et_al_2011_random_samples[i,],
                                            Menge_and_Sutherland_1976_random_samples[i,]))
  
  #list of rows and columns to delete which are not of the 
  #type we are looking for (i.e., "river" and "stream",
  #and unchosen "marine" and "lake" food webs)
  rows_columns_delete <- c()
  
  #columns to delete which are the food web type we are 
  #trying to compare too (i.e., "marine") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
  columns_delete <- c()
  
  #rows to delete which are the food web type we are 
  #trying to compare too (i.e., "lake") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
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
        
        #if the row is either "river" or "stream" food web, we will need 
        #to delete it (also deleting "spring" and "non-aquatic" food webs)
        if (rowOfInterest$Specialized_aquatic_type == "River" | rowOfInterest$Specialized_aquatic_type == "Stream" |
            rowOfInterest$Specialized_aquatic_type == "" | rowOfInterest$Specialized_aquatic_type == "Spring") {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        #if food web is of type "marine"
        } else if (rowOfInterest$Specialized_aquatic_type == "Marine") {
          
          marine_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_marine)) {
            
            #if the "marine" food web is in the realization, only delete from 
            #column
            if (realization_marine[n,]$name==network_name) {
              
              marine_realization <- TRUE
              
              #delete in the column 
              columns_delete <- c(columns_delete,j)
              
            }
            
          }
          
          #if "marine" food web not in the realization, delete from
          #both row and column
          if (marine_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
        #if food web is of type "lake"
        } else if (rowOfInterest$Specialized_aquatic_type == "Lake") {
          
          
          lake_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_lake)) {
            
            #if the "lake" food web is in the realization, only delete from 
            #row
            if (realization_lake[n,]$name==network_name) {
              
              lake_realization <- TRUE
              
              #delete in the row
              rows_delete <- c(rows_delete,j)
              
            }
            
          }
          
          #if "lake" food web not in the realization, delete from
          #both row and column
          if (lake_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
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
  
  #deleting food webs that are not the grouping of interest now
  ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]
  
  #realization number
  average_distance_lake_v_marine[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "lake" food webs are columns, and "marine" food webs are rows
  average_distance_lake_v_marine[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "lake" and in the realization
  average_distance_lake_v_marine[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "marine" and in the realization
  average_distance_lake_v_marine[i,4] <- nrow(ecological_intermediate)
  
}

##
##  Mean pairwise GCD-11 between "lake" 
##  and "river" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_lake_v_river <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_lake_v_river) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

for (i in 1:sample_times){
  
  #obtaining each "lake" food web that was randomly sampled for a specific realization
  realization_lake <- as.data.frame(rbind(single_aquatic_Lake_food_webs,Angelini_et_al_2013_random_samples[i,],
                                          Stewart_and_Sprules_2011_random_samples[i,],Alcorlo_et_al_2001_random_samples[i,],
                                          Cohen_et_al_2003_random_samples[i,],Fryer_1959_random_samples[i,],
                                          Havens_1992_random_samples[i,]))
  
  #obtaining each "river" food web that was randomly sampled for a specific realization
  realization_river <- as.data.frame(rbind(single_aquatic_River_food_webs,Angelini_et_al_2006_random_samples[i,]))
  
  #list of rows and columns to delete which are not of the
  #type we are looking for (i.e., "marine" and "stream",
  #and unchosen "lake" and "river" food webs)
  rows_columns_delete <- c()
  
  #columns to delete which are the food web type we are 
  #trying to compare too (i.e., "river") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
  columns_delete <- c()
  
  #rows to delete which are the food web type we are 
  #trying to compare too (i.e., "lake") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
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
        
        #if the row is either "marine" or "stream" food web, we will need 
        #to delete it (also deleting "spring" and "non-aquatic" networks)
        if (rowOfInterest$Specialized_aquatic_type == "Marine" | rowOfInterest$Specialized_aquatic_type == "Stream" |
            rowOfInterest$Specialized_aquatic_type == "" | rowOfInterest$Specialized_aquatic_type == "Spring") {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        #if food web is of type "river"
        } else if (rowOfInterest$Specialized_aquatic_type == "River") {
          
          river_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_river)) {
            
            #if the "river" food web is in the realization, only delete from 
            #column
            if (realization_river[n,]$name==network_name) {
              
              river_realization <- TRUE
              
              #delete in the column 
              columns_delete <- c(columns_delete,j)
              
            }
            
          }
          
          #if "river" food web not in the realization, delete from
          #both row and column
          if (river_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
        #if food web is of type "lake"
        } else if (rowOfInterest$Specialized_aquatic_type == "Lake") {
          
          
          lake_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_lake)) {
            
            #if the "lake" food web is in the realization, only delete from 
            #row
            if (realization_lake[n,]$name==network_name) {
              
              lake_realization <- TRUE
              
              #delete in the row
              rows_delete <- c(rows_delete,j)
              
            }
            
          }
          
          #if "lake" food web not in the realization, delete from
          #both row and column
          if (lake_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
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
  
  #deleting food webs that are not the grouping of interest now
  ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]
  
  #realization number
  average_distance_lake_v_river[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "lake" food webs are columns, and "river" food webs are rows
  average_distance_lake_v_river[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "lake" and in the realization
  average_distance_lake_v_river[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "river" and in the realization
  average_distance_lake_v_river[i,4] <- nrow(ecological_intermediate)
  
}

##
##  Mean pairwise GCD-11 between "lake" 
##  and "stream" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_lake_v_stream <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_lake_v_stream) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

for (i in 1:sample_times){
  
  #obtaining each "lake" food web that was randomly sampled for a specific realization
  realization_lake <- as.data.frame(rbind(single_aquatic_Lake_food_webs,Angelini_et_al_2013_random_samples[i,],
                                          Stewart_and_Sprules_2011_random_samples[i,],Alcorlo_et_al_2001_random_samples[i,],
                                          Cohen_et_al_2003_random_samples[i,],Fryer_1959_random_samples[i,],
                                          Havens_1992_random_samples[i,]))
  
  #obtaining each "stream" food web that was randomly sampled for a specific realization
  realization_stream <- as.data.frame(rbind(single_aquatic_Stream_food_webs,Closs_and_Lake_1994_random_samples[i,],
                                            Layer_et_al_2010_random_samples[i,],O_Gorman_et_al_2019_random_samples[i,],
                                            Parker_and_Huryn_2006_random_samples[i,],Tavares_Cromar_and_Williams_1996_random_samples[i,],
                                            Thompson_and_Townsend_2003_random_samples[i,],Thompson_and_Townsend_2004_random_samples[i,]))
  
  #list of rows and columns to delete which are not of the
  #type we are looking for (i.e., "marine" and "river",
  #and unchosen "lake" and "stream" food webs)
  rows_columns_delete <- c()
  
  #columns to delete which are the food web type we are 
  #trying to compare too (i.e., "stream") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
  columns_delete <- c()
  
  #rows to delete which are the food web type we are 
  #trying to compare too (i.e., "lake") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
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
        
        #if the row is either "marine" or "river" food web, we will need 
        #to delete it (also deleting "spring" and "non-aquatic" food webs)
        if (rowOfInterest$Specialized_aquatic_type == "Marine" | rowOfInterest$Specialized_aquatic_type == "River" |
            rowOfInterest$Specialized_aquatic_type == "" | rowOfInterest$Specialized_aquatic_type == "Spring") {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        #if food web is of type "stream"
        } else if (rowOfInterest$Specialized_aquatic_type == "Stream") {
          
          stream_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_stream)) {
            
            #if the "stream" food web is in the realization, only delete from 
            #column
            if (realization_stream[n,]$name==network_name) {
              
              stream_realization <- TRUE
              
              #delete in the column 
              columns_delete <- c(columns_delete,j)
              
            }
            
          }
          
          #if "stream" food web not in the realization, delete from
          #both row and column
          if (stream_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
        #if food web is of type "lake"
        } else if (rowOfInterest$Specialized_aquatic_type == "Lake") {
          
          
          lake_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_lake)) {
            
            #if the "lake" food web is in the realization, only delete from 
            #row
            if (realization_lake[n,]$name==network_name) {
              
              lake_realization <- TRUE
              
              #delete in the row
              rows_delete <- c(rows_delete,j)
              
            }
            
          }
          
          #if "lake" food web not in the realization, delete from
          #both row and column
          if (lake_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
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
  
  #deleting food webs that are not the grouping of interest now
  ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]
  
  #realization number
  average_distance_lake_v_stream[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "lake" food webs are columns, and "stream" food webs are rows
  average_distance_lake_v_stream[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "lake" and in the realization
  average_distance_lake_v_stream[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "stream" and in the realization
  average_distance_lake_v_stream[i,4] <- nrow(ecological_intermediate)
  
}

##
##  Mean pairwise GCD-11 between "marine" 
##  and "river" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_marine_v_river <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_marine_v_river) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

for (i in 1:sample_times){
  
  #obtaining each "marine" food web that was randomly sampled for a specific realization
  realization_marine <- as.data.frame(rbind(single_aquatic_Marine_food_webs,Baeta_et_al_2011_random_samples[i,],
                                            Menge_and_Sutherland_1976_random_samples[i,]))
  
  #obtaining each "river" food web that was randomly sampled for a specific realization
  realization_river <- as.data.frame(rbind(single_aquatic_River_food_webs,Angelini_et_al_2006_random_samples[i,]))
  
  #list of rows and columns to delete which are not of the
  #type we are looking for (i.e., "lake" and "stream",
  #and unchosen "marine" and "river" food webs)
  rows_columns_delete <- c()
  
  #columns to delete which are the food web type we are 
  #trying to compare to (i.e., "river") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
  columns_delete <- c()
  
  #rows to delete which are the food web type we are 
  #trying to compare to (i.e., "marine") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
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
        
        #if the row is either "lake" or "stream" food web, we will need 
        #to delete it (also deleting "spring" and "non-aquatic" food webs)
        if (rowOfInterest$Specialized_aquatic_type == "Stream" | rowOfInterest$Specialized_aquatic_type == "Lake" |
            rowOfInterest$Specialized_aquatic_type == "" | rowOfInterest$Specialized_aquatic_type == "Spring") {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        #if food is of type "river"
        } else if (rowOfInterest$Specialized_aquatic_type == "River") {
          
          river_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_river)) {
            
            #if the "river" food web is in the realization, only delete from 
            #column
            if (realization_river[n,]$name==network_name) {
              
              river_realization <- TRUE
              
              #delete in the column 
              columns_delete <- c(columns_delete,j)
              
            }
            
          }
          
          #if "river" food web not in the realization, delete from
          #both row and column
          if (river_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
        #if food web is of type "marine"
        } else if (rowOfInterest$Specialized_aquatic_type == "Marine") {
          
          
          marine_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_marine)) {
            
            #if the "marine" food web is in the realization, only delete from 
            #row
            if (realization_marine[n,]$name==network_name) {
              
              marine_realization <- TRUE
              
              #delete in the row
              rows_delete <- c(rows_delete,j)
              
            }
            
          }
          
          #if "marine" food web not in the realization, delete from
          #both row and column
          if (marine_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
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
  
  #deleting food webs that are not the grouping of interest now
  ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]
  
  #realization number
  average_distance_marine_v_river[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "marine" food webs are columns, and "river" food webs are rows
  average_distance_marine_v_river[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "marine" and in the realization
  average_distance_marine_v_river[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "river" and in the realization
  average_distance_marine_v_river[i,4] <- nrow(ecological_intermediate)
  
}

##
##  Mean pairwise GCD-11 between "marine" 
##  and "stream" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_marine_v_stream <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_marine_v_stream) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

for (i in 1:sample_times){
  
  #obtaining each "marine" food web that was randomly sampled for a specific realization
  realization_marine <- as.data.frame(rbind(single_aquatic_Marine_food_webs,Baeta_et_al_2011_random_samples[i,],
                                            Menge_and_Sutherland_1976_random_samples[i,]))
  
  #obtaining each "stream" food web that was randomly sampled for a specific realization
  realization_stream <- as.data.frame(rbind(single_aquatic_Stream_food_webs,Closs_and_Lake_1994_random_samples[i,],
                                            Layer_et_al_2010_random_samples[i,],O_Gorman_et_al_2019_random_samples[i,],
                                            Parker_and_Huryn_2006_random_samples[i,],Tavares_Cromar_and_Williams_1996_random_samples[i,],
                                            Thompson_and_Townsend_2003_random_samples[i,],Thompson_and_Townsend_2004_random_samples[i,]))
  
  #list of rows and columns to delete which are not of the
  #type we are looking for (i.e., "lake" and "river",
  #and unchosen "marine" and "stream" food webs)
  rows_columns_delete <- c()
  
  #columns to delete which are the food web type we are 
  #trying to compare too (i.e., "stream") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
  columns_delete <- c()
  
  #rows to delete which are the food web type we are 
  #trying to compare too (i.e., "marine") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
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
        
        #if the row is either "lake" or "river" food web, we will need 
        #to delete it (also deleting "spring" and "non-aquatic" networks)
        if (rowOfInterest$Specialized_aquatic_type == "River" | rowOfInterest$Specialized_aquatic_type == "Lake" |
            rowOfInterest$Specialized_aquatic_type == "" | rowOfInterest$Specialized_aquatic_type == "Spring") {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        #if food web is of type "stream"
        } else if (rowOfInterest$Specialized_aquatic_type == "Stream") {
          
          stream_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_stream)) {
            
            #if the "stream" food web is in the realization, only delete from 
            #column
            if (realization_stream[n,]$name==network_name) {
              
              stream_realization <- TRUE
              
              #delete in the column 
              columns_delete <- c(columns_delete,j)
              
            }
            
          }
          
          #if "stream" food web not in the realization, delete from
          #both row and column
          if (stream_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
        #if food web is of type "marine"
        } else if (rowOfInterest$Specialized_aquatic_type == "Marine") {
          
          
          marine_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_marine)) {
            
            #if the "marine" food web is in the realization, only delete from 
            #row
            if (realization_marine[n,]$name==network_name) {
              
              marine_realization <- TRUE
              
              #delete in the row
              rows_delete <- c(rows_delete,j)
              
            }
            
          }
          
          #if "marine" food web not in the realization, delete from
          #both row and column
          if (marine_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
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
  
  #deleting food webs that are not the grouping of interest now
  ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]
  
  #realization number
  average_distance_marine_v_stream[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "marine" food webs are columns, and "stream" food webs are rows
  average_distance_marine_v_stream[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "marine" and in the realization
  average_distance_marine_v_stream[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "stream" and in the realization
  average_distance_marine_v_stream[i,4] <- nrow(ecological_intermediate)
  
}

##
##  Mean pairwise GCD-11 between "river" 
##  and "stream" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_river_v_stream <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_river_v_stream) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

for (i in 1:sample_times){
  
  #obtaining each "river" food web that was randomly sampled for a specific realization
  realization_river <- as.data.frame(rbind(single_aquatic_River_food_webs,Angelini_et_al_2006_random_samples[i,]))
  
  #obtaining each "stream" food web that was randomly sampled for a specific realization
  realization_stream <- as.data.frame(rbind(single_aquatic_Stream_food_webs,Closs_and_Lake_1994_random_samples[i,],
                                            Layer_et_al_2010_random_samples[i,],O_Gorman_et_al_2019_random_samples[i,],
                                            Parker_and_Huryn_2006_random_samples[i,],Tavares_Cromar_and_Williams_1996_random_samples[i,],
                                            Thompson_and_Townsend_2003_random_samples[i,],Thompson_and_Townsend_2004_random_samples[i,]))
  
  #list of rows and columns to delete which are not of the
  #type we are looking for (i.e., "lake" and "marine",
  #and unchosen "river" and "stream" food webs)
  rows_columns_delete <- c()
  
  #columns to delete which are the food web type we are 
  #trying to compare too (i.e., "stream") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
  columns_delete <- c()
  
  #rows to delete which are the food web type we are 
  #trying to compare too (i.e., "river") and are in 
  #the realization so as not to compare food webs 
  #twice (i.e, they are in both the rows and columns)
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
        
        #if the row is either "lake" or "marine" food web, we will need 
        #to delete it (also deleting "spring" and "non-aquatic" food webs)
        if (rowOfInterest$Specialized_aquatic_type == "Lake" | rowOfInterest$Specialized_aquatic_type == "Marine" |
            rowOfInterest$Specialized_aquatic_type == "" | rowOfInterest$Specialized_aquatic_type == "Spring") {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        #if food web is of type "stream"
        } else if (rowOfInterest$Specialized_aquatic_type == "Stream") {
          
          stream_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_stream)) {
            
            #if the "stream" food web is in the realization, only delete from 
            #column
            if (realization_stream[n,]$name==network_name) {
              
              stream_realization <- TRUE
              
              #delete in the column 
              columns_delete <- c(columns_delete,j)
              
            }
            
          }
          
          #if "stream" food web not in the realization, delete from
          #both row and column
          if (stream_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
        #if food web is of type "river"
        } else if (rowOfInterest$Specialized_aquatic_type == "River") {
          
          
          river_realization <- FALSE
          
          #looping through to see if food web is in the realization
          for (n in 1:nrow(realization_river)) {
            
            #if the "river" food web is in the realization, only delete from 
            #row
            if (realization_river[n,]$name==network_name) {
              
              river_realization <- TRUE
              
              #delete in the row
              rows_delete <- c(rows_delete,j)
              
            }
            
          }
          
          #if "river" food web not in the realization, delete from
          #both row and column
          if (river_realization==FALSE) {
            
            rows_columns_delete <- c(rows_columns_delete,j)
            
          }
          
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
  
  #deleting food webs that are not the grouping of interest now
  ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]
  
  #realization number
  average_distance_river_v_stream[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "river" food webs are columns, and "stream" food webs are rows
  average_distance_river_v_stream[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "river" and in the realization
  average_distance_river_v_stream[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "stream" and in the realization
  average_distance_river_v_stream[i,4] <- nrow(ecological_intermediate)
  
}

##
##  All information for table S2
##

#mean pairwise GCD-11 of "lake" food webs
round(mean(average_distance_lake[,2]),2)

#mean pairwise GCD-11 of "marine" food webs 
round(mean(average_distance_marine[,2]),2)

#mean pairwise GCD-11 of "river" food webs
round(mean(average_distance_river[,2]),2)

#mean pairwise GCD-11 of "stream" food webs
round(mean(average_distance_stream[,2]),2)

#mean pairwise GCD-11 between "lake" and "marine" food webs
round(mean(average_distance_lake_v_marine[,2]),2)

#mean pairwise GCD-11 between "lake" and "river" food webs
round(mean(average_distance_lake_v_river[,2]),2)

#mean pairwise GCD-11 between "lake" and "stream" food webs
round(mean(average_distance_lake_v_stream[,2]),2)

#mean pairwise GCD-11 between "marine" and "river" food webs
round(mean(average_distance_marine_v_river[,2]),2)

#mean pairwise GCD-11 between "marine" and "stream" food webs
round(mean(average_distance_marine_v_stream[,2]),2)

#mean pairwise GCD-11 between "river" and "stream" food webs
round(mean(average_distance_river_v_stream[,2]),2)