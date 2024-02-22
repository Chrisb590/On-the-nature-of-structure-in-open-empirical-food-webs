########################################################
#
# For recreating the results of Brimacombe et al. (2024):
# On the nature of structure in collections of freely 
# available food webs
#
# This script takes both the pairwise GCD11 matrix 
# (gcd11.csv) between all networks and information about 
# each network (metaData.csv) and calculates the median 
# pairwise GCD-11 between networks
#
########################################################

#clearing environment
rm(list = ls())

#loading necessary packages
library(ggplot2)
library(igraph)
library(dplyr)
library(spatstat)

#set working directory where gcd11.csv file is
setwd(".")

#reading in the pairwise GCD-11 matrix
ecological <- as.matrix(read.csv("gcd11.csv",row.names = 1))

#the associated metadata of all networks
metadata <- as.matrix(read.csv("metaData.csv"))
metadata <- as.data.frame(metadata)
metadata$author <- as.factor(metadata$author)
metadata$number_of_nodes <- as.numeric(metadata$number_of_nodes)
metadata$year_published <- as.numeric(metadata$year_published)
metadata$decade_published <- as.factor(metadata$decade_published)
metadata$Ecopath <- as.factor(metadata$Ecopath)

########################################################
##
##
##
##  Used to evaluate the median pairwise GCD-11 between 
##  networks from the same publication grouping as
##  shown in parts of Table S3
##
##
##
########################################################

#storage for median pairwise GCD-11 between networks from the same publication grouping
each_median_distance <- matrix(,nrow=(length(levels(metadata$author))),ncol=3)
colnames(each_median_distance) <- c("network_number","median_distance","number_of_networks")

#the different publication groupings
for (i in 1:nrow(each_median_distance)) {
  
  each_median_distance[i,1] <- levels(metadata$author)[i]
  
}

#looping through the pairwise GCD-11 matrix to find 
#networks associated with specific publication groupings
for (i in 1:nrow(each_median_distance)) {
  
  #list of rows and columns to delete which are not of the 
  #publication grouping of interest right now
  rows_columns_delete <- c()
  
  for (j in 1:ncol(ecological)) {
    
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the network is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata[rowNumber,])
      
      #if the row name in the matrix is the network we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not a publication network we want, we will need 
        #to delete it (here, we just keep the location of this network
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$author != each_median_distance[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
        #otherwise we have not found the network in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
    }
    
  }
  
  #deleting networks that are not the publication grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #median pairwise GCD-11 between specific publication networks
  each_median_distance[i,2] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the specific publication grouping of interest
  each_median_distance[i,3] <- ncol(ecological_intermediate)
  
}

each_median_distance <- as.data.frame(each_median_distance)

#median pairwise GCD-11 of networks from the same publication (in Table S3)
webs_from_pubs_with_multiple_networks <- each_median_distance[!(each_median_distance$network_number=="One_network_per_publication"),]
round(median(as.numeric(webs_from_pubs_with_multiple_networks$median_distance)),2)

#median pairwise GCD-11 of networks from publications that only produced a single network (in Table S3)
webs_from_pubs_with_single_network <- each_median_distance[each_median_distance$network_number=="One_network_per_publication",]
webs_from_pubs_with_single_network

########################################################
##
##
##
##  Used to evaluate the median pairwise GCD-11 between
##  networks from publications that produced multiple
##  networks, before 1990 and after 1990 (Figure S7 and
##  parts of Table S3)
##
##
##
########################################################

#only for networks that are sourced from publications that each produced only
#a single network 
metadata_one_network_per_publication <- metadata[complete.cases(metadata),] 

#storage for median pairwise GCD-11 between networks sourced from publications
#that each provided only a single network, and binned by decade of publication
each_median_distance <- matrix(,nrow=(length(levels(metadata_one_network_per_publication$decade_published))),ncol=4)
colnames(each_median_distance) <- c("publication_year","median_distance","number_of_networks","Identity")

for (i in 1:nrow(each_median_distance)) {
  
  each_median_distance[i,1] <- levels(metadata_one_network_per_publication$decade_published)[i]
  
}

#looping through the pairwise GCD-11 matrix for networks sourced from
#publications that each provided only a single network, and binned by 
#specific decade of publication
for (i in 1:nrow(each_median_distance)) {
  
  #list of rows and columns to delete which are not of the 
  #decade of interest
  rows_columns_delete <- c()
  
  for (j in 1:ncol(ecological)) {
    
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the network is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata_one_network_per_publication[rowNumber,])
      
      #if the row name in the matrix is the network we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not a decade we are interested in, we will need 
        #to delete it (here, we just keep the location of this network
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$decade_published != each_median_distance[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
        #otherwise we have not found the network in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
      #if network not in metadata, we do not want it
      if (rowNumber > nrow(metadata_one_network_per_publication)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
    }
    
  }
  
  #deleting networks that are not the decade of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #median pairwise GCD-11 between specific decade
  each_median_distance[i,2] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks of the specific decade
  each_median_distance[i,3] <- ncol(ecological_intermediate)
  
  #identifier
  each_median_distance[i,4] <- "Decade Published"
  
}

each_median_distance <- as.data.frame(each_median_distance)
each_median_distance$number_of_networks <- as.numeric(each_median_distance$number_of_networks)
each_median_distance$publication_year <- as.numeric(each_median_distance$publication_year)
each_median_distance$median_distance <- as.numeric(each_median_distance$median_distance)

##
##  Same as above but for the median GCD-11 between networks from publications
##  that provided multiple networks
##

#storage for median GCD-11 between networks from the same publication
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

#looping through the pairwise GCD-11 matrix to find networks 
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
    
    #looping to find what number the network is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata[rowNumber,])
      
      #if the row name in the matrix is the network we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not network we want, we will need 
        #to delete it (here, we just keep the location of this network
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$author != each_median_distance_for_publication[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
        #otherwise we have not found the network in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
    }
    
  }
  
  #deleting networks that are not the publication grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #median pairwise GCD-11 between specific publication's networks
  each_median_distance_for_publication[i,3] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the specific publication grouping of interest
  each_median_distance_for_publication[i,4] <- ncol(ecological_intermediate)
  
}

each_median_distance_for_publication <- as.data.frame(each_median_distance_for_publication)
each_median_distance_for_publication$median_distance <- as.numeric(each_median_distance_for_publication$median_distance)
each_median_distance_for_publication$publication_year <- as.numeric(each_median_distance_for_publication$publication_year)
each_median_distance_for_publication$number_of_networks <- as.numeric(each_median_distance_for_publication$number_of_networks)

#data needed for plotting
data_for_plotting <- rbind(each_median_distance_for_publication[,2:5],each_median_distance)

#need to weight interpolation that is to be used in the plot for 2003, 2006, 2011
new_for_interpolation <- matrix(,nrow=length(unique(each_median_distance_for_publication$publication_year)),ncol=2)
colnames(new_for_interpolation) <- c("year","weighted_median")

for (i in 1:length(unique(each_median_distance_for_publication$publication_year))) {
  
  new_for_interpolation[i,1] <- unique(each_median_distance_for_publication$publication_year)[i]
  
  #for year 2006
  if (unique(each_median_distance_for_publication$publication_year)[i] == 2006) {
    
    weighted_median_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year==2006,]
    new_for_interpolation[i,2] <- weighted.median(weighted_median_intermediate$median_distance,weighted_median_intermediate$number_of_networks)
    
    #for year 2003
  } else if (unique(each_median_distance_for_publication$publication_year)[i] == 2003) {
    
    weighted_median_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year==2003,]
    new_for_interpolation[i,2] <- weighted.median(weighted_median_intermediate$median_distance,weighted_median_intermediate$number_of_networks)
    
    #for year 2011
  } else if (unique(each_median_distance_for_publication$publication_year)[i] == 2011) {
    
    weighted_median_intermediate <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year==2011,]
    new_for_interpolation[i,2] <- weighted.median(weighted_median_intermediate$median_distance,weighted_median_intermediate$number_of_networks)
    
    #otherwise no need to interpolate year
  } else {
    
    new_for_interpolation[i,2] <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year==new_for_interpolation[i,1],]$median_distance
    
  }
  
  
}

new_for_interpolation <- as.data.frame(new_for_interpolation)

#Figure S7 in the manuscript
ggplot(data=data_for_plotting,aes(x=publication_year,y=median_distance,color=Identity,size=number_of_networks)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values=c("#029386", "#011288")) +
  scale_fill_manual(values=c("black", "black")) + 
  scale_size_continuous(range=c(0,10),limits=c(0,30),breaks=c(2,3,4,5,10,20,30)) +
  geom_line(data=each_median_distance,aes(x=publication_year,y=median_distance),color="#029386",size=0.7) +
  geom_line(data=new_for_interpolation,aes(x=year,y=weighted_median),color="#011288",size=0.7) + 
  theme_bw() +
  xlab("Publication year") +
  ylab("Median GCD-11") +
  theme(axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=20, colour ="black"),
        axis.title.y = element_text(size=20, colour ="black"),
        text = element_text(size = 20,color="black"),
        axis.text.x = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color="black"))

#median of the median pairwise GCD-11 from publications before 1990 (parts of Table S3)
each_median_distance_for_publication_before_1990 <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year<1990,]
round(median(each_median_distance_for_publication_before_1990$median_distance),2)

#median of the median pairwise GCD-11 from publications after 1990 (parts of Table S3)
each_median_distance_for_publication_after_1990 <- each_median_distance_for_publication[each_median_distance_for_publication$publication_year>1990,]
round(median(each_median_distance_for_publication_after_1990$median_distance),2)

########################################################
##
##
##
##  Used to evaluate the median pairwise GCD-11 between
##  networks from the same ecological system (Parts
##  of Table S2)
##
##
##
########################################################

#storage for median pairwise GCD-11 between networks from the same ecosystem
each_median_distance <- matrix(,nrow=(length(levels(as.factor(metadata$Primary_type)))),ncol=3)
colnames(each_median_distance) <- c("Ecosystem_type","median_distance","number_of_networks")

for (i in 1:nrow(each_median_distance)) {
  
  each_median_distance[i,1] <- levels(as.factor(metadata$Primary_type))[i]
  
}

#looping through the pairwise GCD-11 matrix to find networks associated 
#with specific ecosystems
for (i in 1:nrow(each_median_distance)) {
  
  #list of rows and columns to delete which are not of the 
  #ecosystem grouping of interest right now
  rows_columns_delete <- c()
  
  for (j in 1:ncol(ecological)) {
    
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the network is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata[rowNumber,])
      
      #if the row name in the matrix is the network we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not an ecosystem network we want, we will need 
        #to delete it (here, we just keep the location of this network
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$Primary_type != each_median_distance[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        #otherwise we have not found the network in the metadata list
        
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
      #if network not in metadata, we do not want it
      if (rowNumber > nrow(metadata)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
    }
    
  }
  
  #deleting networks that are not the ecosystem grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #median pairwise GCD-11 between specific ecosystem networks
  each_median_distance[i,2] <- round(median(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the specific ecosystem grouping of interest
  each_median_distance[i,3] <- ncol(ecological_intermediate)
  
}

each_median_distance <- as.data.frame(each_median_distance)
each_median_distance$median_distance <- as.numeric(each_median_distance$median_distance)
each_median_distance$number_of_networks <- as.numeric(each_median_distance$number_of_networks)

#parts of Table S2
each_median_distance

########################################################
##
##
##
##  Used to evaluate the median pairwise GCD-11 between
##  networks from different ecological systems (Parts 
##  of Table S2)
##
##
##
########################################################

##
##  Median pairwise GCD-11 between "Aquatic" and 
##  "Aquatic and terrestrial" networks
##

#list of rows and columns to delete which are not of the 
#type we are looking for (i.e., "Terrestrial")
rows_columns_delete <- c()

#columns to delete which are the network type we are 
#trying to compare too (i.e., "Aquatic") so as not to
#compare networks twice (i.e, they are in both the rows
#and columns)
columns_delete <- c()

#rows to delete which are the network type we are 
#trying to compare too (i.e., "Aquatic and terrestrial") 
#so as not to compare networks twice (i.e, they are in 
#both the rows and columns)
rows_delete <- c()

for (j in 1:ncol(ecological)) {
  
  network_name <- row.names(ecological)[j]
  
  found = FALSE
  rowNumber = 1
  
  #looping to find what number the network is
  while (found != TRUE) {
    
    rowOfInterest <- as.data.frame(metadata[rowNumber,])
    
    #if the row name in the matrix is the network we are looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is a terrestrial network, we will need 
      #to delete it
      if (rowOfInterest$Primary_type == "Terrestrial") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Aquatic") {
        
        columns_delete <- c(columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Aquatic and terrestrial") {
        
        rows_delete <- c(rows_delete,j)
        
      }
      
      found = TRUE
      
      #otherwise we have not found the network in the metadata list
    } else {
      
      rowNumber = rowNumber + 1
      
    }
    
    #if network not in metadata, we do not want it
    if (rowNumber > nrow(metadata)) {
      
      rows_columns_delete <- c(rows_columns_delete,j)
      
      found = TRUE
      
    }
    
  }
  
}

#deleting networks that are not of interest
ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]

#median pairwise GCD-11 between "aquatic" and "aquatic and terrestrial" (parts of Table S2)
round(median(ecological_intermediate),2)

##
##  Median pairwise GCD-11 between "Aquatic" and 
##  "Terrestrial" networks
##

#list of rows and columns to delete which are not of the 
#type we are looking for (i.e., "Aquatic and terrestrial")
rows_columns_delete <- c()

#columns to delete which are the network type we are 
#trying to compare too (i.e., "Aquatic") so as not to
#compare networks twice (i.e, they are in both the rows
#and columns)
columns_delete <- c()

#rows to delete which are the network type we are 
#trying to compare too (i.e., "Terrestrial") so as 
#not to compare networks twice (i.e, they are in 
#both the rows and columns)
rows_delete <- c()

for (j in 1:ncol(ecological)) {
  
  network_name <- row.names(ecological)[j]
  
  found = FALSE
  rowNumber = 1
  
  #looping to find what number the network is
  while (found != TRUE) {
    
    rowOfInterest <- as.data.frame(metadata[rowNumber,])
    
    #if the row name in the matrix is the network we are looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is a Aquatic and terrestrial network, we will need 
      #to delete it
      if (rowOfInterest$Primary_type == "Aquatic and terrestrial") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Aquatic") {
        
        columns_delete <- c(columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Terrestrial") {
        
        rows_delete <- c(rows_delete,j)
        
      }
      
      found = TRUE
      
      #otherwise we have not found the network in the metadata list
    } else {
      
      rowNumber = rowNumber + 1
      
    }
    
    #if network not in metadata, we do not want it
    if (rowNumber > nrow(metadata)) {
      
      rows_columns_delete <- c(rows_columns_delete,j)
      
      found = TRUE
      
    }
    
  }
  
}

#deleting networks that are not of interest
ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]

#median pairwise GCD-11 between "aquatic" and "terrestrial" (parts of Table S2)
round(median(ecological_intermediate),2)

##
##  Median pairwise GCD-11 between "Aquatic and terrestrial" 
##  and "Terrestrial" networks
##

#list of rows and columns to delete which are not of the 
#type we are looking for (i.e., "Aquatic")
rows_columns_delete <- c()

#columns to delete which are the network type we are 
#trying to compare too (i.e., "Aquatic and terrestrial") 
#so as not to compare networks twice (i.e, they are in 
#both the rows and columns)
columns_delete <- c()

#rows to delete which are the network type we are 
#trying to compare too (i.e., "Terrestrial") so as 
#not to compare networks twice (i.e, they are in 
#both the rows and columns)
rows_delete <- c()

for (j in 1:ncol(ecological)) {
  
  network_name <- row.names(ecological)[j]
  
  found = FALSE
  rowNumber = 1
  
  #looping to find what number the network is
  while (found != TRUE) {
    
    rowOfInterest <- as.data.frame(metadata[rowNumber,])
    
    #if the row name in the matrix is the network we are looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is a Aquatic network, we will need 
      #to delete it
      if (rowOfInterest$Primary_type == "Aquatic") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Aquatic and terrestrial") {
        
        columns_delete <- c(columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Terrestrial") {
        
        rows_delete <- c(rows_delete,j)
        
      }
      
      found = TRUE
      
      #otherwise we have not found the network in the metadata list
    } else {
      
      rowNumber = rowNumber + 1
      
    }
    
    #if network not in metadata, we do not want it
    if (rowNumber > nrow(metadata)) {
      
      rows_columns_delete <- c(rows_columns_delete,j)
      
      found = TRUE
      
    }
    
  }
  
}

#deleting networks that are not of interest
ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]

#median pairwise distance between "aquatic and terrestrial" and "terrestrial" (parts of Table S2)
round(median(ecological_intermediate),2)
