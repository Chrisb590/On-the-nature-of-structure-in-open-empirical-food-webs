########################################################
##
##
##
##  For recreating Appendix S1: Section S1.8 of 
##  Brimacombe et al. (2024): On the nature of structure 
##  in collections of freely available food webs
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
sample_times <- 50

########################################################

##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs that are "aquatic" and not created using 
##  Ecopath (Appendix S1: Section S1.8.1)
##

#Randomly obtaining one food web per publication grouping

#Closs_and_Lake_1994 non-Ecopath "aquatic" food webs
Closs_and_Lake_1994_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Closs_and_Lake_1994_random_samples) <- c(colnames(metadata))
Closs_and_Lake_1994 <-  metadata[metadata$author=="Closs_and_Lake_1994",]

for (i in 1:sample_times) {
  
  Closs_and_Lake_1994_random_samples[i,] <- as.matrix(sample_n(Closs_and_Lake_1994, 1))
  
}

#Thompson_and_Townsend_2003 non-Ecopath "aquatic" food webs
Thompson_and_Townsend_2003_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Thompson_and_Townsend_2003_random_samples) <- c(colnames(metadata))
Thompson_and_Townsend_2003 <-  metadata[metadata$author=="Thompson_and_Townsend_2003",]

for (i in 1:sample_times) {
  
  Thompson_and_Townsend_2003_random_samples[i,] <- as.matrix(sample_n(Thompson_and_Townsend_2003, 1))
  
}

#Alcorlo_et_al_2001 non-Ecopath "aquatic" food webs
Alcorlo_et_al_2001_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Alcorlo_et_al_2001_random_samples) <- c(colnames(metadata))
Alcorlo_et_al_2001 <-  metadata[metadata$author=="Alcorlo_et_al_2001",]

for (i in 1:sample_times) {
  
  Alcorlo_et_al_2001_random_samples[i,] <- as.matrix(sample_n(Alcorlo_et_al_2001, 1))
  
}

#Menge_and_Sutherland_1976 non-Ecopath "aquatic" food webs
Menge_and_Sutherland_1976_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Menge_and_Sutherland_1976_random_samples) <- c(colnames(metadata))
Menge_and_Sutherland_1976 <-  metadata[metadata$author=="Menge_and_Sutherland_1976",]

for (i in 1:sample_times) {
  
  Menge_and_Sutherland_1976_random_samples[i,] <- as.matrix(sample_n(Menge_and_Sutherland_1976, 1))
  
}

#Tavares_Cromar_and_Williams_1996 non-Ecopath "aquatic" food webs
Tavares_Cromar_and_Williams_1996_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Tavares_Cromar_and_Williams_1996_random_samples) <- c(colnames(metadata))
Tavares_Cromar_and_Williams_1996 <-  metadata[metadata$author=="Tavares-Cromar_and_Williams_1996",]

for (i in 1:sample_times) {
  
  Tavares_Cromar_and_Williams_1996_random_samples[i,] <- as.matrix(sample_n(Tavares_Cromar_and_Williams_1996, 1))
  
}

#Parker_and_Huryn_2006 non-Ecopath "aquatic" food webs
Parker_and_Huryn_2006_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Parker_and_Huryn_2006_random_samples) <- c(colnames(metadata))
Parker_and_Huryn_2006 <-  metadata[metadata$author=="Parker_and_Huryn_2006",]

for (i in 1:sample_times) {
  
  Parker_and_Huryn_2006_random_samples[i,] <- as.matrix(sample_n(Parker_and_Huryn_2006, 1))
  
}

#Thompson_and_Townsend_2004 non-Ecopath "aquatic" food webs
Thompson_and_Townsend_2004_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Thompson_and_Townsend_2004_random_samples) <- c(colnames(metadata))
Thompson_and_Townsend_2004 <-  metadata[metadata$author=="Thompson_and_Townsend_2004",]

for (i in 1:sample_times) {
  
  Thompson_and_Townsend_2004_random_samples[i,] <- as.matrix(sample_n(Thompson_and_Townsend_2004, 1))
  
}

#Cohen_et_al_2003 non-Ecopath "aquatic" food webs
Cohen_et_al_2003_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Cohen_et_al_2003_random_samples) <- c(colnames(metadata))
Cohen_et_al_2003 <-  metadata[metadata$author=="Cohen_et_al_2003",]

for (i in 1:sample_times) {
  
  Cohen_et_al_2003_random_samples[i,] <- as.matrix(sample_n(Cohen_et_al_2003, 1))
  
}

#Fryer_1959 non-Ecopath "aquatic" food webs
Fryer_1959_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Fryer_1959_random_samples) <- c(colnames(metadata))
Fryer_1959 <-  metadata[metadata$author=="Fryer_1959",]

for (i in 1:sample_times) {
  
  Fryer_1959_random_samples[i,] <- as.matrix(sample_n(Fryer_1959, 1))
  
}

#Havens_1992 non-Ecopath "aquatic" food webs
Havens_1992_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Havens_1992_random_samples) <- c(colnames(metadata))
Havens_1992 <-  metadata[metadata$author=="Havens_1992",]
Havens_1992 <-  Havens_1992[Havens_1992$Primary_type!="Aquatic and terrestrial",]

for (i in 1:sample_times) {
  
  Havens_1992_random_samples[i,] <- as.matrix(sample_n(Havens_1992, 1))
  
}

#Layer_et_al_2010 non-Ecopath "aquatic" food webs
Layer_et_al_2010_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Havens_1992_random_samples) <- c(colnames(metadata))
Layer_et_al_2010 <-  metadata[metadata$author=="Layer_et_al_2010",]

for (i in 1:sample_times) {
  
  Layer_et_al_2010_random_samples[i,] <- as.matrix(sample_n(Layer_et_al_2010, 1))
  
}

#O_Gorman_et_al_2019 non-Ecopath "aquatic" food webs
O_Gorman_et_al_2019_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(O_Gorman_et_al_2019_random_samples) <- c(colnames(metadata))
O_Gorman_et_al_2019 <-  metadata[metadata$author=="O_Gorman_et_al_2019",]

for (i in 1:sample_times) {
  
  O_Gorman_et_al_2019_random_samples[i,] <- as.matrix(sample_n(O_Gorman_et_al_2019, 1))
  
}

#for non-Ecopath "aquatic" food webs sourced from publications that provided only a single network
single_aquatic_not_ecopath_food_webs <- metadata[metadata$author=="One_network_per_publication",]
single_aquatic_not_ecopath_food_webs <- single_aquatic_not_ecopath_food_webs[single_aquatic_not_ecopath_food_webs$Primary_type=="Aquatic",]
single_aquatic_not_ecopath_food_webs <- single_aquatic_not_ecopath_food_webs[single_aquatic_not_ecopath_food_webs$Ecopath=="Not_Ecopath_model",]

#storage for mean pairwise GCD-11 between food webs that are non-Ecopath "aquatic" from particular realization
average_distance_non_ecopath <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_non_ecopath) <- c("realization","mean_distance","number_of_networks")

#row counter for average_distance_non_ecopath
row_counter_for_average_distance_non_ecopath <- 1

for (i in 1:sample_times){
  
  #obtaining each food web from publications that were randomly sampled for a specific realization
  realization <- as.data.frame(rbind(single_aquatic_not_ecopath_food_webs,Closs_and_Lake_1994_random_samples[i,], Thompson_and_Townsend_2003_random_samples[i,],
                                     Alcorlo_et_al_2001_random_samples[i,], Menge_and_Sutherland_1976_random_samples[i,], Tavares_Cromar_and_Williams_1996_random_samples[i,],
                                     Parker_and_Huryn_2006_random_samples[i,], Thompson_and_Townsend_2004_random_samples[i,],Cohen_et_al_2003_random_samples[i,],
                                     Fryer_1959_random_samples[i,],Havens_1992_random_samples[i,],Layer_et_al_2010_random_samples[i,],
                                     O_Gorman_et_al_2019_random_samples[i,]))
  
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
  average_distance_non_ecopath[row_counter_for_average_distance_non_ecopath,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_non_ecopath[row_counter_for_average_distance_non_ecopath,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the realization
  average_distance_non_ecopath[row_counter_for_average_distance_non_ecopath,3] <- ncol(ecological_intermediate)
  
  row_counter_for_average_distance_non_ecopath = row_counter_for_average_distance_non_ecopath + 1
  
}

##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs that are "aquatic" and created using 
##  Ecopath (Appendix S1: Section S1.8.2)
##

#Obtaining all 72 combinations of Ecopath "aquatic" food webs
#(i.e., removing for publication effect)

#Angelini_et_al_2006 Ecopath "aquatic" food webs
Angelini_et_al_2006_storage <- matrix(,nrow=0,ncol(metadata))
colnames(Angelini_et_al_2006_storage) <- c(colnames(metadata))
Angelini_et_al_2006 <-  metadata[metadata$author=="Angelini_et_al_2006",]

for (i in 1:36) {
  
  Angelini_et_al_2006_storage <- rbind(Angelini_et_al_2006_storage,Angelini_et_al_2006)
  
}

#Stewart_&_Sprules_2011 Ecopath "aquatic" food webs
Stewart_and_Sprules_storage <- matrix(,nrow=0,ncol(metadata))
colnames(Stewart_and_Sprules_storage) <- c(colnames(metadata))
Stewart_and_Sprules_2011 <-  metadata[metadata$author=="Stewart_and_Sprules_2011",]

for (i in 1:72) {
  
  if (i%%4==1 | i%%4==2) {
    
    Stewart_and_Sprules_storage <- rbind(Stewart_and_Sprules_storage,Stewart_and_Sprules_2011[1,])
    
  } else {
    
    Stewart_and_Sprules_storage <- rbind(Stewart_and_Sprules_storage,Stewart_and_Sprules_2011[2,])
    
  }
  
  
  
}

#Angelini_et_al_2013 Ecopath "aquatic" food webs
Angelini_et_al_2013_storage <- matrix(,nrow=0,ncol(metadata))
colnames(Angelini_et_al_2013_storage) <- c(colnames(metadata))
Angelini_et_al_2013 <-  metadata[metadata$author=="Angelini_et_al_2013",]

for (i in 1:72) {
  
  if (i%%12==1 | i%%12==2 | i%%12==3 | i%%12==4) {
    
    Angelini_et_al_2013_storage <- rbind(Angelini_et_al_2013_storage,Angelini_et_al_2013[1,])
    
  } else if (i%%12==5 | i%%12==6 | i%%12==7 | i%%12==8) {
    
    Angelini_et_al_2013_storage <- rbind(Angelini_et_al_2013_storage,Angelini_et_al_2013[2,])
    
  } else {
    
    Angelini_et_al_2013_storage <- rbind(Angelini_et_al_2013_storage,Angelini_et_al_2013[3,])
    
  }
  
  
}

#Baeta_et_al_2011 Ecopath "aquatic" food webs
Baeta_et_al_2011_storage <- matrix(0,nrow=0,ncol(metadata))
colnames(Baeta_et_al_2011_storage) <- c(colnames(metadata))
Baeta_et_al_2011 <-  metadata[metadata$author=="Baeta_et_al_2011",]

for (i in 1:72) {
  
  if (i%%72==1 | i%%72==2 | i%%72==3 | i%%72==4 | i%%72==5 | i%%72==6 | i%%72==7 | i%%72==8 | i%%72==9 | i%%72==10 | 
      i%%72==11 | i%%72==12) {
    
    Baeta_et_al_2011_storage <- rbind(Baeta_et_al_2011_storage,Baeta_et_al_2011[1,])
    
  } else if (i%%72==13 | i%%72==14 | i%%72==15 | i%%72==16 | i%%72==17 | i%%72==18 | i%%72==19 | i%%72==20 | i%%72==21 | i%%72==22 | 
             i%%72==23 | i%%72==24) {
    
    Baeta_et_al_2011_storage <- rbind(Baeta_et_al_2011_storage,Baeta_et_al_2011[2,])
    
  } else if (i%%72==25 | i%%72==26 | i%%72==27 | i%%72==28 | i%%72==29 | i%%72==30 | i%%72==31 | i%%72==32 | i%%72==33 | i%%72==34 | 
             i%%72==35 | i%%72==36) {
    
    Baeta_et_al_2011_storage <- rbind(Baeta_et_al_2011_storage,Baeta_et_al_2011[3,])
    
  } else if (i%%72==37 | i%%72==38 | i%%72==39 | i%%72==40 | i%%72==41 | i%%72==42 | i%%72==43 | i%%72==44 | i%%72==45 | i%%72==46 | 
             i%%72==47 | i%%72==48) {
    
    Baeta_et_al_2011_storage <- rbind(Baeta_et_al_2011_storage,Baeta_et_al_2011[4,])
    
  } else if (i%%72==49 | i%%72==50 | i%%72==51 | i%%72==52 | i%%72==53 | i%%72==54 | i%%72==55 | i%%72==56 | i%%72==57 | i%%72==58 | 
             i%%72==59 | i%%72==60) {
    
    Baeta_et_al_2011_storage <- rbind(Baeta_et_al_2011_storage,Baeta_et_al_2011[5,])
    
  } else {
    
    Baeta_et_al_2011_storage <- rbind(Baeta_et_al_2011_storage,Baeta_et_al_2011[6,])
    
  }
  
}

#food webs that are Ecopath "aquatic" but from publications that produced only a single network
single_food_webs_ecopath <- metadata[metadata$author=="One_network_per_publication",]
single_food_webs_ecopath <- single_food_webs_ecopath[single_food_webs_ecopath$Ecopath=="Ecopath_model",]

#storage for all 72 possible combinations of mean pairwise GCD-11 between 
#food webs that are Ecopath "aquatic"
average_distance_ecopath <- matrix(,nrow=72,ncol=3)
colnames(average_distance_ecopath) <- c("realization","mean_distance","number_of_networks")

#row counter for average_distance_ecopath
row_counter_for_average_distance_ecopath <- 1

for (i in 1:72){

  #obtaining each publication's sampled food web for a specific realization, and all food webs
  #sourced from publications that each provided only a single network
  realization <- as.data.frame(rbind(single_food_webs_ecopath,Angelini_et_al_2006_storage[i,],Angelini_et_al_2013_storage[i,],
                                     Baeta_et_al_2011_storage[i,],Stewart_and_Sprules_storage[i,]))
  
  #list of rows and columns to delete which are not of the 
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
  average_distance_ecopath[row_counter_for_average_distance_ecopath,1] <- i
  
  #mean pairwise GCD-11 between food webs that are Ecopath "aquatic"
  average_distance_ecopath[row_counter_for_average_distance_ecopath,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are Ecopath "aquatic"
  average_distance_ecopath[row_counter_for_average_distance_ecopath,3] <- ncol(ecological_intermediate)
  
  row_counter_for_average_distance_ecopath = row_counter_for_average_distance_ecopath + 1
  
}

##
##  All information for Appendix S1: Section S1.8
##

#value for mean pairwise GCD-11 between non-Ecopath "aquatic" food webs
#with publication effect removed
round(mean(average_distance_non_ecopath[,2]),2)

#value for mean pairwise GCD-11 between Ecopath "aquatic" food webs
#with publication effect removed
round(mean(average_distance_ecopath[,2]),2)
