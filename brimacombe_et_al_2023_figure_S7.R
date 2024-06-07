########################################################
##
##
##
##  For recreating figure S7 of Brimacombe et al. (2024):
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

#all pairwise GCD-11s (only lower triangular matrix since
#it is symmetric [don't want to count twice])
ecological_intermediate <- ecological[lower.tri(ecological)]

#storage for number of counts of each pairwise GCD-11
storage_for_entries <- as.matrix(matrix(0,nrow=3,ncol=15))
colnames(storage_for_entries) <- c("0.0<=x<=0.5","0.5<x<=1.0","1.0<x<=1.5","1.5<x<=2.0",
                                   "2.0<x<=2.5","2.5<x<=3.0","3.0<x<=3.5","3.5<x<=4.0",
                                   "4.0<x<=4.5","4.5<x<=5.0","5.0<x<=5.5","5.5<x<=6.0",
                                   "6.0<x<=6.5","6.5<x<=7.0","7.0<x<=7.5")
rownames(storage_for_entries) <- c("Only between food webs from publications that produced a single network (n=3403)",
                                   "Only between food webs from publications that produced multiple networks (n=2487)",
                                   "Between all other possible food web pairwise comparisons (n=31511)")

row_of_storage <- 3

#
# All pairwise GCD-11s (regardless of grouping)
# 

#determining pairwise GCD-11 bin
for (k in 1:length(ecological_intermediate)) {
  
  if (ecological_intermediate[k] <= 0.5) {
    
    storage_for_entries[row_of_storage,1] <- storage_for_entries[row_of_storage,1] + 1
    
  } else if (ecological_intermediate[k] > 0.5 & ecological_intermediate[k] <= 1.0) {
    
    storage_for_entries[row_of_storage,2] <- storage_for_entries[row_of_storage,2] + 1
    
  } else if (ecological_intermediate[k] > 1.0 & ecological_intermediate[k] <= 1.5) {
    
    storage_for_entries[row_of_storage,3] <- storage_for_entries[row_of_storage,3] + 1
    
  } else if (ecological_intermediate[k] > 1.5 & ecological_intermediate[k] <= 2.0) {
    
    storage_for_entries[row_of_storage,4] <- storage_for_entries[row_of_storage,4] + 1
    
  } else if (ecological_intermediate[k] > 2.0 & ecological_intermediate[k] <= 2.5) {
    
    storage_for_entries[row_of_storage,5] <- storage_for_entries[row_of_storage,5] + 1
    
  } else if (ecological_intermediate[k] > 2.5 & ecological_intermediate[k] <= 3.0) {
    
    storage_for_entries[row_of_storage,6] <- storage_for_entries[row_of_storage,6] + 1
    
  } else if (ecological_intermediate[k] > 3.0 & ecological_intermediate[k] <= 3.5) {
    
    storage_for_entries[row_of_storage,7] <- storage_for_entries[row_of_storage,7] + 1
    
  } else if (ecological_intermediate[k] > 3.5 & ecological_intermediate[k] <= 4.0) {
    
    storage_for_entries[row_of_storage,8] <- storage_for_entries[row_of_storage,8] + 1
    
  } else if (ecological_intermediate[k] > 4.0 & ecological_intermediate[k] <= 4.5) {
    
    storage_for_entries[row_of_storage,9] <- storage_for_entries[row_of_storage,9] + 1
    
  } else if (ecological_intermediate[k] > 4.5 & ecological_intermediate[k] <= 5.0) {
    
    storage_for_entries[row_of_storage,10] <- storage_for_entries[row_of_storage,10] + 1
    
  } else if (ecological_intermediate[k] > 5.0 & ecological_intermediate[k] <= 5.5) {
    
    storage_for_entries[row_of_storage,11] <- storage_for_entries[row_of_storage,11] + 1
    
  } else if (ecological_intermediate[k] > 5.5 & ecological_intermediate[k] <= 6.0) {
    
    storage_for_entries[row_of_storage,12] <- storage_for_entries[row_of_storage,12] + 1
    
  } else if (ecological_intermediate[k] > 6.0 & ecological_intermediate[k] <= 6.5) {
    
    storage_for_entries[row_of_storage,13] <- storage_for_entries[row_of_storage,13] + 1
    
  } else if (ecological_intermediate[k] > 6.5 & ecological_intermediate[k] <= 7.0) {
    
    storage_for_entries[row_of_storage,14] <- storage_for_entries[row_of_storage,14] + 1
    
  } else {
    
    storage_for_entries[row_of_storage,15] <- storage_for_entries[row_of_storage,15] + 1
    
  }
  
}

#
# All pairwise GCD-11s between (i) publications that provided multiple networks or
# (2) pariwise GCD-11s from publications that provided only a single network 
# 

#storage for pairwise GCD-11s between food webs from the same publication grouping
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
  
  #deleting food webs that are not the publication grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  ecological_intermediate <- ecological_intermediate[lower.tri(ecological_intermediate)]
  
  #if searching for "one_network_per_publication", then store in first row
  if (each_average_distance[i]=="One_network_per_publication") {
    
    row_of_storage <- 1
    
  } else {
    
    row_of_storage <- 2
    
  }
  
  for (k in 1:length(ecological_intermediate)) {
    
    if (ecological_intermediate[k] <= 0.5) {
      
      storage_for_entries[row_of_storage,1] <- storage_for_entries[row_of_storage,1] + 1
      
    } else if (ecological_intermediate[k] > 0.5 & ecological_intermediate[k] <= 1.0) {
      
      storage_for_entries[row_of_storage,2] <- storage_for_entries[row_of_storage,2] + 1
      
    } else if (ecological_intermediate[k] > 1.0 & ecological_intermediate[k] <= 1.5) {
      
      storage_for_entries[row_of_storage,3] <- storage_for_entries[row_of_storage,3] + 1
      
    } else if (ecological_intermediate[k] > 1.5 & ecological_intermediate[k] <= 2.0) {
      
      storage_for_entries[row_of_storage,4] <- storage_for_entries[row_of_storage,4] + 1
      
    } else if (ecological_intermediate[k] > 2.0 & ecological_intermediate[k] <= 2.5) {
      
      storage_for_entries[row_of_storage,5] <- storage_for_entries[row_of_storage,5] + 1
      
    } else if (ecological_intermediate[k] > 2.5 & ecological_intermediate[k] <= 3.0) {
      
      storage_for_entries[row_of_storage,6] <- storage_for_entries[row_of_storage,6] + 1
      
    } else if (ecological_intermediate[k] > 3.0 & ecological_intermediate[k] <= 3.5) {
      
      storage_for_entries[row_of_storage,7] <- storage_for_entries[row_of_storage,7] + 1
      
    } else if (ecological_intermediate[k] > 3.5 & ecological_intermediate[k] <= 4.0) {
      
      storage_for_entries[row_of_storage,8] <- storage_for_entries[row_of_storage,8] + 1
      
    } else if (ecological_intermediate[k] > 4.0 & ecological_intermediate[k] <= 4.5) {
      
      storage_for_entries[row_of_storage,9] <- storage_for_entries[row_of_storage,9] + 1
      
    } else if (ecological_intermediate[k] > 4.5 & ecological_intermediate[k] <= 5.0) {
      
      storage_for_entries[row_of_storage,10] <- storage_for_entries[row_of_storage,10] + 1
      
    } else if (ecological_intermediate[k] > 5.0 & ecological_intermediate[k] <= 5.5) {
      
      storage_for_entries[row_of_storage,11] <- storage_for_entries[row_of_storage,11] + 1
      
    } else if (ecological_intermediate[k] > 5.5 & ecological_intermediate[k] <= 6.0) {
      
      storage_for_entries[row_of_storage,12] <- storage_for_entries[row_of_storage,12] + 1
      
    } else if (ecological_intermediate[k] > 6.0 & ecological_intermediate[k] <= 6.5) {
      
      storage_for_entries[row_of_storage,13] <- storage_for_entries[row_of_storage,13] + 1
      
    } else if (ecological_intermediate[k] > 6.5 & ecological_intermediate[k] <= 7.0) {
      
      storage_for_entries[row_of_storage,14] <- storage_for_entries[row_of_storage,14] + 1
      
    } else {
      
      storage_for_entries[row_of_storage,15] <- storage_for_entries[row_of_storage,15] + 1
      
    }
    
  }
  
}

storage_for_entries_2 <- storage_for_entries

#all pairwise GCD-11s "Between all other possible food web pairwise comparisons (n=31511)"
storage_for_entries_2[3,] <- storage_for_entries_2[3,] - storage_for_entries_2[1,] - storage_for_entries_2[2,]

intermediate <- as_tibble(as.data.frame(as.table(storage_for_entries_2)))
plotting_data <- as.data.frame(intermediate)

#
# Figure S7
# 

ggplot(plotting_data, aes(fill=Var1, x=Var2, y=Freq)) + 
  geom_bar(position="stack",stat = "identity") +
  xlab("Pairwise GCD-11") +
  ylab("Frequency") + 
  scale_fill_manual(values=c("#85D4E3","#F4B5BD","#FAD77B")) +
  theme(axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=15, colour ="black"),
        axis.title.y = element_text(size=15, colour ="black"),
        text = element_text(size = 15,color="black"),
        axis.text.x = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color="black"))
