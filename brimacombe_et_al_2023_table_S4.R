########################################################
##
##
##
##  For recreating table S4 of Brimacombe et al. (2024):
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

#storage for median pairwise GCD-11 between food webs from the same publication grouping
each_median_distance_publication_grouping <- matrix(,nrow=(length(levels(metadata$author))),ncol=3)
colnames(each_median_distance_publication_grouping) <- c("network_number","median_distance","number_of_networks")

#the different publication groupings
for (i in 1:nrow(each_median_distance_publication_grouping)) {
  
  each_median_distance_publication_grouping[i,1] <- levels(metadata$author)[i]
  
}

#looping through the pairwise GCD-11 matrix to find 
#food webs source from specific publication groupings
for (i in 1:nrow(each_median_distance_publication_grouping)) {
  
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
        if (rowOfInterest$author != each_median_distance_publication_grouping[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
    }
    
  }
  
  #deleting food webs that are not sourced from the publication grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #median pairwise GCD-11 between specific publication's food webs
  each_median_distance_publication_grouping[i,2] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific publication grouping of interest
  each_median_distance_publication_grouping[i,3] <- ncol(ecological_intermediate)
  
}

each_median_distance_publication_grouping <- as.data.frame(each_median_distance_publication_grouping)

##
##  Used to evaluate the median pairwise GCD-11 between
##  food webs sourced from publications that produced 
##  multiple networks, before 1990 and after 1990  
##

#only for food webs sourced from publications that each produced only
#a single network 
metadata_one_network_per_publication <- metadata[complete.cases(metadata),] 

#storage for median pairwise GCD-11 between food webs sourced from publications
#that each produced only a single network, and binned by decade of publication
each_median_distance <- matrix(,nrow=(length(levels(metadata_one_network_per_publication$decade_published))),ncol=4)
colnames(each_median_distance) <- c("publication_year","median_distance","number_of_networks","Identity")

for (i in 1:nrow(each_median_distance)) {
  
  each_median_distance[i,1] <- levels(metadata_one_network_per_publication$decade_published)[i]
  
}

#looping through the pairwise GCD-11 matrix for food webs sourced from
#publications that each produced only a single network, and binned by 
#specific decade of publication
for (i in 1:nrow(each_median_distance)) {
  
  #list of rows and columns to delete which are not of the 
  #decade of interest
  rows_columns_delete <- c()
  
  for (j in 1:ncol(ecological)) {
    
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata_one_network_per_publication[rowNumber,])
      
      #if the row name in the matrix is the food web we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not a decade we are interested in, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$decade_published != each_median_distance[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
      #if food web not in metadata, we do not want it
      if (rowNumber > nrow(metadata_one_network_per_publication)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
    }
    
  }
  
  #deleting food webs that are not the decade of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #median pairwise GCD-11 between food webs of specific decade
  each_median_distance[i,2] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs of the specific decade
  each_median_distance[i,3] <- ncol(ecological_intermediate)
  
  #identifier
  each_median_distance[i,4] <- "Decade Published"
  
}

each_median_distance <- as.data.frame(each_median_distance)
each_median_distance$number_of_networks <- as.numeric(each_median_distance$number_of_networks)
each_median_distance$publication_year <- as.numeric(each_median_distance$publication_year)
each_median_distance$median_distance <- as.numeric(each_median_distance$median_distance)

##
##  Same as above but for the median GCD-11 between food webs from publications
##  that provided multiple networks
##

#storage for median GCD-11 between food webs sourced from the same publication that produced multiple networks
each_median_distance_for_publication <- matrix(,nrow=(length(levels(metadata$author))),ncol=5)
colnames(each_median_distance_for_publication) <- c("publication_name","publication_year","median_distance","number_of_networks","Identity")

#obtaining information about all publications that provided multiple networks
for (i in 1:nrow(each_median_distance_for_publication)) {
  
  each_median_distance_for_publication[i,1] <- levels(metadata$author)[i]
  each_median_distance_for_publication[i,2] <- metadata[which(metadata$author==levels(metadata$author)[i])[1],]$year_published
  each_median_distance_for_publication[i,5] <- "Publication"
  
}

each_median_distance_for_publication <- as.data.frame(each_median_distance_for_publication)
each_median_distance_for_publication <- each_median_distance_for_publication[!(each_median_distance_for_publication$publication_name=="One_network_per_publication"),]

#looping through the pairwise GCD-11 matrix to find food webs 
#associated with specific publications that provided multiple
#networks
for (i in 1:nrow(each_median_distance_for_publication)) {
  
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
        
        #if the row is not food web we want, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$author != each_median_distance_for_publication[i,1]) {
          
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
  
  #median pairwise GCD-11 between specific publication's food webs
  each_median_distance_for_publication[i,3] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific publication grouping of interest
  each_median_distance_for_publication[i,4] <- ncol(ecological_intermediate)
  
}

each_median_distance_for_publication <- as.data.frame(each_median_distance_for_publication)
each_median_distance_for_publication$median_distance <- as.numeric(each_median_distance_for_publication$median_distance)
each_median_distance_for_publication$publication_year <- as.numeric(each_median_distance_for_publication$publication_year)
each_median_distance_for_publication$number_of_networks <- as.numeric(each_median_distance_for_publication$number_of_networks)

decades <- c(1950,1970,1980,1990,2000,2010)

#need to weight interpolation that is to be used in the plot for averaging across decade 
new_for_interpolation <- matrix(,nrow=length(unique(c(1950,1970,1980,1990,2000,2010))),ncol=4)

#putting in the years in the decade column
for (i in 1:nrow(new_for_interpolation)) {
  
  new_for_interpolation[i,1] <- decades[i] 
  
}

colnames(new_for_interpolation) <- c("publication_year","median_distance","number_of_networks","Identity")

new_for_interpolation <- as.data.frame(new_for_interpolation)

for (i in 1:length(unique(new_for_interpolation$publication_year))) {
  
  #for year 1950
  if (new_for_interpolation$publication_year[i]==1950) {
    
    weighted_average_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year<1960 &
                                                                            each_median_distance_for_publication$publication_year>=1950,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$median_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- "Publication"
    
    #for year 1970
  } else if (new_for_interpolation$publication_year[i]==1970) {
    
    weighted_average_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year<1980 &
                                                                            each_median_distance_for_publication$publication_year>=1970,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$median_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- "Publication"
    
    #for year 1980
  } else if (new_for_interpolation$publication_year[i]==1980) {
    
    weighted_average_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year<1990 &
                                                                            each_median_distance_for_publication$publication_year>=1980,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$median_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- "Publication"
    
    #for year 1990
  } else if (new_for_interpolation$publication_year[i]==1990) {
    
    weighted_average_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year<2000 &
                                                                            each_median_distance_for_publication$publication_year>=1990,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$median_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- "Publication"
    
    #for year 2000
  } else if (new_for_interpolation$publication_year[i]==2000) {
    
    weighted_average_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year<2010 &
                                                                            each_median_distance_for_publication$publication_year>=2000,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$median_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- "Publication"
    
    #for 2010
  } else if (new_for_interpolation$publication_year[i]==2010) {
    
    weighted_average_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year<2020 &
                                                                            each_median_distance_for_publication$publication_year>=2010,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$median_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- "Publication"
    
  }
  
  
}

##
##  All information for table S4
##

#median pairwise GCD-11 of food webs from publications that only produced a single network
webs_from_pubs_with_single_network <- each_median_distance_publication_grouping[each_median_distance_publication_grouping$network_number=="One_network_per_publication",]
webs_from_pubs_with_single_network

#median pairwise GCD-11 of food webs from the same publication that produced multiple networks
webs_from_pubs_with_multiple_networks <- each_median_distance_publication_grouping[!(each_median_distance_publication_grouping$network_number=="One_network_per_publication"),]
round(median(as.numeric(webs_from_pubs_with_multiple_networks$median_distance)),2)

#median of the median pairwise GCD-11 from publications that produced multiple networks before 2000
each_median_distance_for_publication_before_1990 <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year<2000,]
round(median(each_median_distance_for_publication_before_1990$median_distance),2)

#median of the median pairwise GCD-11 from publications that produced multiple networks after 2000
each_median_distance_for_publication_after_1990 <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year>=2000,]
round(median(each_median_distance_for_publication_after_1990$median_distance),2)
