########################################################
##
##
##
##  For recreating figure 4 of Brimacombe et al. (2024):
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

#only for food webs sourced from publications that 
#each produced only a single network 
metadata_one_network_per_publication <- metadata[complete.cases(metadata),] 

#storage for mean pairwise GCD-11 between food webs sourced from publications
#that each provided only a single network, and binned by decade of publication
each_average_distance <- matrix(,nrow=(length(levels(metadata_one_network_per_publication$decade_published))),ncol=5)
colnames(each_average_distance) <- c("publication_year","mean_distance","SD_distance","number_of_networks","Identity")

for (i in 1:nrow(each_average_distance)) {
  
  each_average_distance[i,1] <- levels(metadata_one_network_per_publication$decade_published)[i]
  
}

#looping through the pairwise GCD-11 matrix for food webs sourced from
#publications that each provided only a single network, and binned by 
#specific decade of publication
for (i in 1:nrow(each_average_distance)) {
  
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
        if (rowOfInterest$decade_published != each_average_distance[i,1]) {
          
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
  
  #mean pairwise GCD-11 between food webs of specific decade
  each_average_distance[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #standard deviation in pairwise GCD-11 between food webs of specific decade
  each_average_distance[i,3] <- round(sd(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs of the specific decade
  each_average_distance[i,4] <- ncol(ecological_intermediate)
  
  #identifier
  each_average_distance[i,5] <- "Decade Published"
  
}

each_average_distance <- as.data.frame(each_average_distance)
each_average_distance$number_of_networks <- as.numeric(each_average_distance$number_of_networks)
each_average_distance$publication_year <- as.numeric(each_average_distance$publication_year)
each_average_distance$mean_distance <- as.numeric(each_average_distance$mean_distance)
each_average_distance$SD_distance <- as.numeric(each_average_distance$SD_distance)

##
##  Same as above but for the mean pairwise GCD-11 between food webs sourced from 
##  publications that provided multiple networks
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
        #food web in the GCD-11 matrix)
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

#replacing with zeros
each_average_distance_for_publication[is.na(each_average_distance_for_publication)] <- 0

decades <- c(1950,1970,1980,1990,2000,2010)

#need to weight interpolation that is to be used in the plot for averaging across decade 
new_for_interpolation <- matrix(,nrow=length(unique(c(1950,1970,1980,1990,2000,2010))),ncol=5)

#putting in the years in the decade column
for (i in 1:nrow(new_for_interpolation)) {
  
  new_for_interpolation[i,1] <- decades[i] 
  
}

colnames(new_for_interpolation) <- c("publication_year","mean_distance","SD_distance","number_of_networks","Identity")

new_for_interpolation <- as.data.frame(new_for_interpolation)

for (i in 1:length(unique(new_for_interpolation$publication_year))) {
  
  #for year 1950
  if (new_for_interpolation$publication_year[i]==1950) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<1960 &
                                                                             each_average_distance_for_publication$publication_year>=1950,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- weighted.mean(weighted_average_intermediate$SD_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,5] <- "Publication"
    
    #for year 1970
  } else if (new_for_interpolation$publication_year[i]==1970) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<1980 &
                                                                             each_average_distance_for_publication$publication_year>=1970,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- weighted.mean(weighted_average_intermediate$SD_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,5] <- "Publication"
    
    #for year 1980
  } else if (new_for_interpolation$publication_year[i]==1980) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<1990 &
                                                                             each_average_distance_for_publication$publication_year>=1980,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- weighted.mean(weighted_average_intermediate$SD_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,5] <- "Publication"
    
    #for year 1990
  } else if (new_for_interpolation$publication_year[i]==1990) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<2000 &
                                                                             each_average_distance_for_publication$publication_year>=1990,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- weighted.mean(weighted_average_intermediate$SD_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,5] <- "Publication"
    
    #for year 2000
  } else if (new_for_interpolation$publication_year[i]==2000) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<2010 &
                                                                             each_average_distance_for_publication$publication_year>=2000,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- weighted.mean(weighted_average_intermediate$SD_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,5] <- "Publication"
    
    #for 2010
  } else if (new_for_interpolation$publication_year[i]==2010) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<2020 &
                                                                             each_average_distance_for_publication$publication_year>=2010,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,3] <- weighted.mean(weighted_average_intermediate$SD_distance,weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,4] <- sum(weighted_average_intermediate$number_of_networks)
    new_for_interpolation[i,5] <- "Publication"
    
  }
  
}

new_for_interpolation <- as.data.frame(new_for_interpolation)

#data needed for plotting
data_for_plotting <- rbind(new_for_interpolation,each_average_distance)

##
##  Figure 4
##

ggplot(data=data_for_plotting,aes(x=publication_year,y=mean_distance,color=Identity,size=number_of_networks)) +
  geom_errorbar(data=data_for_plotting,aes(ymin = mean_distance-SD_distance, ymax = mean_distance+SD_distance),width = 0.9,size=1)+
  geom_point(alpha = 0.5) +
  scale_x_continuous(breaks = seq(1910,2020,by=10)) + 
  scale_color_manual(values=c("#029386", "#011288")) +
  scale_fill_manual(values=c("black", "black")) + 
  scale_size_continuous(range=c(0,15),limits=c(0,80),breaks=c(2,5,10,20,40,60,80)) +
  geom_line(data=each_average_distance,aes(x=publication_year,y=mean_distance),color="#029386",size=0.7) +
  geom_line(data=new_for_interpolation,aes(x=publication_year,y=mean_distance),color="#011288",size=0.7) + 
  theme_bw() +
  xlab("Decade published") +
  ylab("Mean pairwise GCD-11") +
  theme(axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=20, colour ="black"),
        axis.title.y = element_text(size=20, colour ="black"),
        text = element_text(size = 20,color="black"),
        axis.text.x = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color="black"))
