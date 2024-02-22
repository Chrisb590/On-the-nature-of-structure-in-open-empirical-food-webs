########################################################
#
# For recreating the results of Brimacombe et al. (2024):
# On the nature of structure in collections of freely 
# available food webs 
#
# This script takes both the pairwise GCD11 matrix 
# (gcd11.csv) between all networks and information about 
# each network (metaData.csv) and calculates the mean 
# pairwise GCD-11 between networks
#
########################################################

#clearing environment
rm(list = ls())

#loading necessary packages
library(ggplot2)
library(igraph)
library(dplyr)

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
##  Used to evaluate the mean pairwise GCD-11 between 
##  networks from the same publication grouping as
##  shown in parts of Table 2, Table S5
##
##
##
########################################################

#storage for mean pairwise GCD-11 between networks from the same publication grouping
each_average_distance <- matrix(,nrow=(length(levels(metadata$author))),ncol=3)
colnames(each_average_distance) <- c("network_number","mean_distance","number_of_networks")

#the different publication groupings
for (i in 1:nrow(each_average_distance)) {
  
  each_average_distance[i,1] <- levels(metadata$author)[i]
  
}

#looping through the pairwise GCD-11 matrix to find 
#networks associated with specific publication groupings
for (i in 1:nrow(each_average_distance)) {
  
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
        if (rowOfInterest$author != each_average_distance[i,1]) {
          
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
  
  #mean pairwise GCD-11 between specific publication networks
  each_average_distance[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the specific publication grouping of interest
  each_average_distance[i,3] <- ncol(ecological_intermediate)
  
}

#mean pairwise GCD-11 (in Table S5)
each_average_distance <- as.data.frame(each_average_distance)
each_average_distance$mean_distance <- as.numeric(each_average_distance$mean_distance)
each_average_distance$number_of_networks <- as.numeric(each_average_distance$number_of_networks)
each_average_distance

#weighted mean pairwise GCD-11 of networks from the same publication (in Table 2)
webs_from_pubs_with_multiple_networks <- each_average_distance[!(each_average_distance$network_number=="One_network_per_publication"),]
round(weighted.mean(webs_from_pubs_with_multiple_networks$mean_distance,webs_from_pubs_with_multiple_networks$number_of_networks),2)

#mean pairwise GCD-11 of networks from publications that only produced a single network (in Table 2)
webs_from_pubs_with_single_network <- each_average_distance[each_average_distance$network_number=="One_network_per_publication",]
webs_from_pubs_with_single_network

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 as a
##  function of the mean number of nodes per network
##  and standard deviation in the number of nodes per
##  network for publications that provided multiple
##  networks (in Figure S8)
##
##
##
########################################################

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

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  networks from publications that produced multiple
##  networks, before 1990 and after 1990 (Figure 4 and
##  parts of Table 2)
##
##
##
########################################################

#only for networks that are sourced from publications that each produced only
#a single network 
metadata_one_network_per_publication <- metadata[complete.cases(metadata),] 

#storage for mean pairwise GCD-11 between networks sourced from publications
#that each provided only a single network, and binned by decade of publication
each_average_distance <- matrix(,nrow=(length(levels(metadata_one_network_per_publication$decade_published))),ncol=4)
colnames(each_average_distance) <- c("publication_year","mean_distance","number_of_networks","Identity")

for (i in 1:nrow(each_average_distance)) {
  
  each_average_distance[i,1] <- levels(metadata_one_network_per_publication$decade_published)[i]
  
}

#looping through the pairwise GCD-11 matrix for networks sourced from
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
    
    #looping to find what number the network is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata_one_network_per_publication[rowNumber,])
      
      #if the row name in the matrix is the network we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not a decade we are interested in, we will need 
        #to delete it (here, we just keep the location of this network
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$decade_published != each_average_distance[i,1]) {
          
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
  
  #mean pairwise GCD-11 between specific decade
  each_average_distance[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks of the specific decade
  each_average_distance[i,3] <- ncol(ecological_intermediate)
  
  #identifier
  each_average_distance[i,4] <- "Decade Published"
  
}

each_average_distance <- as.data.frame(each_average_distance)
each_average_distance$number_of_networks <- as.numeric(each_average_distance$number_of_networks)
each_average_distance$publication_year <- as.numeric(each_average_distance$publication_year)
each_average_distance$mean_distance <- as.numeric(each_average_distance$mean_distance)


##
##  Same as above but for the mean GCD-11 between networks from publications
##  that provided multiple networks
##

#storage for mean GCD-11 between networks from the same publication
each_average_distance_for_publication <- matrix(,nrow=(length(levels(metadata$author))),ncol=5)
colnames(each_average_distance_for_publication) <- c("publication_name","publication_year","mean_distance","number_of_networks","Identity")

#obtaining information about all publications that provided multiple networks
for (i in 1:nrow(each_average_distance_for_publication)) {
  
  each_average_distance_for_publication[i,1] <- levels(metadata$author)[i]
  each_average_distance_for_publication[i,2] <- metadata[which(metadata$author==levels(metadata$author)[i])[1],]$year_published
  each_average_distance_for_publication[i,5] <- "Publication"
  
}

each_average_distance_for_publication <- as.data.frame(each_average_distance_for_publication)
each_average_distance_for_publication <- each_average_distance_for_publication[!(each_average_distance_for_publication$publication_name=="One_network_per_publication"),]

#looping through the pairwise GCD-11 matrix to find networks 
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
    
    #looping to find what number the network is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(metadata[rowNumber,])
      
      #if the row name in the matrix is the network we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not network we want, we will need 
        #to delete it (here, we just keep the location of this network
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$author != each_average_distance_for_publication[i,1]) {
          
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
  
  #mean pairwise GCD-11 between specific publication's networks
  each_average_distance_for_publication[i,3] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the specific publication grouping of interest
  each_average_distance_for_publication[i,4] <- ncol(ecological_intermediate)
  
}

each_average_distance_for_publication <- as.data.frame(each_average_distance_for_publication)
each_average_distance_for_publication$mean_distance <- as.numeric(each_average_distance_for_publication$mean_distance)
each_average_distance_for_publication$publication_year <- as.numeric(each_average_distance_for_publication$publication_year)
each_average_distance_for_publication$number_of_networks <- as.numeric(each_average_distance_for_publication$number_of_networks)

#data needed for plotting
data_for_plotting <- rbind(each_average_distance_for_publication[,2:5],each_average_distance)

#need to weight interpolation that is to be used in the plot for 2003, 2006, 2011
new_for_interpolation <- matrix(,nrow=length(unique(each_average_distance_for_publication$publication_year)),ncol=2)
colnames(new_for_interpolation) <- c("year","weighted_average")

for (i in 1:length(unique(each_average_distance_for_publication$publication_year))) {
  
  new_for_interpolation[i,1] <- unique(each_average_distance_for_publication$publication_year)[i]
  
  #for year 2006
  if (unique(each_average_distance_for_publication$publication_year)[i] == 2006) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year==2006,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
    
  #for year 2003
  } else if (unique(each_average_distance_for_publication$publication_year)[i] == 2003) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year==2003,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
    
  #for year 2011
  } else if (unique(each_average_distance_for_publication$publication_year)[i] == 2011) {
    
    weighted_average_intermediate <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year==2011,]
    new_for_interpolation[i,2] <- weighted.mean(weighted_average_intermediate$mean_distance,weighted_average_intermediate$number_of_networks)
  
  #otherwise no need to interpolate year
  } else {
    
    new_for_interpolation[i,2] <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year==new_for_interpolation[i,1],]$mean_distance
    
  }
  
  
}

new_for_interpolation <- as.data.frame(new_for_interpolation)


#Figure 4 in the manuscript
ggplot(data=data_for_plotting,aes(x=publication_year,y=mean_distance,color=Identity,size=number_of_networks)) +
  geom_point(alpha = 0.5) +
  scale_color_manual(values=c("#029386", "#011288")) +
  scale_fill_manual(values=c("black", "black")) + 
  scale_size_continuous(range=c(0,10),limits=c(0,30),breaks=c(2,3,4,5,10,20,30)) +
  geom_line(data=each_average_distance,aes(x=publication_year,y=mean_distance),color="#029386",size=0.7) +
  geom_line(data=new_for_interpolation,aes(x=year,y=weighted_average),color="#011288",size=0.7) + 
  theme_bw() +
  xlab("Publication year") +
  ylab("Mean GCD-11") +
  theme(axis.line = element_line(colour = "black"),
        axis.title.x = element_text(size=20, colour ="black"),
        axis.title.y = element_text(size=20, colour ="black"),
        text = element_text(size = 20,color="black"),
        axis.text.x = element_text(color="black"),
        axis.ticks = element_line(color = "black"),
        axis.text.y = element_text(color="black"))

#weighted mean of mean pairwise GCD-11 from publications before 1990 (parts of Table 2)
each_average_distance_for_publication_before_1990 <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<1990,]
round(weighted.mean(each_average_distance_for_publication_before_1990$mean_distance,each_average_distance_for_publication_before_1990$number_of_networks),2)

#weighted mean of mean pairwise GCD-11 from publications after 1990 (parts of Table 2)
each_average_distance_for_publication_after_1990 <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year>1990,]
round(weighted.mean(each_average_distance_for_publication_after_1990$mean_distance,each_average_distance_for_publication_after_1990$number_of_networks),2)

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  networks from the same ecological system (Parts
##  of Table 1)
##
##
##
########################################################

#storage for mean pairwise GCD-11 between networks from the same ecosystem
each_average_distance <- matrix(,nrow=(length(levels(as.factor(metadata$Primary_type)))),ncol=3)
colnames(each_average_distance) <- c("Ecosystem_type","mean_distance","number_of_networks")

for (i in 1:nrow(each_average_distance)) {
  
  each_average_distance[i,1] <- levels(as.factor(metadata$Primary_type))[i]
  
}

#looping through the pairwise GCD-11 matrix to find networks associated 
#with specific ecosystems
for (i in 1:nrow(each_average_distance)) {
  
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
        if (rowOfInterest$Primary_type != each_average_distance[i,1]) {
          
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
  
  #mean pairwise GCD-11 between specific ecosystem networks
  each_average_distance[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the specific ecosystem grouping of interest
  each_average_distance[i,3] <- ncol(ecological_intermediate)
  
}

each_average_distance <- as.data.frame(each_average_distance)
each_average_distance$mean_distance <- as.numeric(each_average_distance$mean_distance)
each_average_distance$number_of_networks <- as.numeric(each_average_distance$number_of_networks)

#parts of Table 1
each_average_distance

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  networks from different ecological systems (Parts 
##  of Table 1)
##
##
##
########################################################

##
##  Mean pairwise GCD-11 between "Aquatic" and 
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

#mean pairwise GCD-11 between "aquatic" and "aquatic and terrestrial" (parts of Table 1)
round(mean(ecological_intermediate),2)

##
##  Mean pairwise GCD-11 between "Aquatic" and 
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

#mean pairwise GCD-11 between "aquatic" and "terrestrial" (parts of Table 1)
round(mean(ecological_intermediate),2)

##
##  Mean pairwise GCD-11 between "Aquatic and terrestrial" 
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
 
#mean pairwise distance between "aquatic and terrestrial" and "terrestrial" (parts of Table 1)
round(mean(ecological_intermediate),2)

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  networks that are aquatic and not created using 
##  Ecopath (Appendix S1.6.1)
##
##
##
########################################################

##
##  Randomly obtaining one network per publication grouping
##

#number of times you want to randomly sample/number of realizations
sample_times <- 200

#Closs_and_Lake_1994 non-Ecopath aquatic networks
Closs_and_Lake_1994_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Closs_and_Lake_1994_random_samples) <- c(colnames(metadata))
Closs_and_Lake_1994 <-  metadata[metadata$author=="Closs_and_Lake_1994",]

for (i in 1:sample_times) {
  
  Closs_and_Lake_1994_random_samples[i,] <- as.matrix(sample_n(Closs_and_Lake_1994, 1))
  
}

#Thompson_and_Townsend_2003 non-Ecopath aquatic networks
Thompson_and_Townsend_2003_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Thompson_and_Townsend_2003_random_samples) <- c(colnames(metadata))
Thompson_and_Townsend_2003 <-  metadata[metadata$author=="Thompson_and_Townsend_2003",]

for (i in 1:sample_times) {
  
  Thompson_and_Townsend_2003_random_samples[i,] <- as.matrix(sample_n(Thompson_and_Townsend_2003, 1))
  
}

#Alcorlo_et_al_2001 non-Ecopath aquatic networks
Alcorlo_et_al_2001_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Alcorlo_et_al_2001_random_samples) <- c(colnames(metadata))
Alcorlo_et_al_2001 <-  metadata[metadata$author=="Alcorlo_et_al_2001",]

for (i in 1:sample_times) {
  
  Alcorlo_et_al_2001_random_samples[i,] <- as.matrix(sample_n(Alcorlo_et_al_2001, 1))
  
}

#Menge_and_Sutherland_1976 non-Ecopath aquatic networks
Menge_and_Sutherland_1976_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Menge_and_Sutherland_1976_random_samples) <- c(colnames(metadata))
Menge_and_Sutherland_1976 <-  metadata[metadata$author=="Menge_and_Sutherland_1976",]

for (i in 1:sample_times) {
  
  Menge_and_Sutherland_1976_random_samples[i,] <- as.matrix(sample_n(Menge_and_Sutherland_1976, 1))
  
}

#Tavares_Cromar_and_Williams_1996 non-Ecopath aquatic networks
Tavares_Cromar_and_Williams_1996_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Tavares_Cromar_and_Williams_1996_random_samples) <- c(colnames(metadata))
Tavares_Cromar_and_Williams_1996 <-  metadata[metadata$author=="Tavares-Cromar_and_Williams_1996",]

for (i in 1:sample_times) {
  
  Tavares_Cromar_and_Williams_1996_random_samples[i,] <- as.matrix(sample_n(Tavares_Cromar_and_Williams_1996, 1))
  
}

#Parker_and_Huryn_2006 non-Ecopath aquatic networks
Parker_and_Huryn_2006_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Parker_and_Huryn_2006_random_samples) <- c(colnames(metadata))
Parker_and_Huryn_2006 <-  metadata[metadata$author=="Parker_and_Huryn_2006",]

for (i in 1:sample_times) {
  
  Parker_and_Huryn_2006_random_samples[i,] <- as.matrix(sample_n(Parker_and_Huryn_2006, 1))
  
}

#Thompson_and_Townsend_2004 non-Ecopath aquatic networks
Thompson_and_Townsend_2004_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Thompson_and_Townsend_2004_random_samples) <- c(colnames(metadata))
Thompson_and_Townsend_2004 <-  metadata[metadata$author=="Thompson_and_Townsend_2004",]

for (i in 1:sample_times) {
  
  Thompson_and_Townsend_2004_random_samples[i,] <- as.matrix(sample_n(Thompson_and_Townsend_2004, 1))
  
}

#Cohen_et_al_2003 non-Ecopath aquatic networks
Cohen_et_al_2003_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Cohen_et_al_2003_random_samples) <- c(colnames(metadata))
Cohen_et_al_2003 <-  metadata[metadata$author=="Cohen_et_al_2003",]

for (i in 1:sample_times) {
  
  Cohen_et_al_2003_random_samples[i,] <- as.matrix(sample_n(Cohen_et_al_2003, 1))
  
}

#Fryer_1959 non-Ecopath aquatic networks
Fryer_1959_random_samples <- matrix(0,sample_times,ncol(metadata))
colnames(Fryer_1959_random_samples) <- c(colnames(metadata))
Fryer_1959 <-  metadata[metadata$author=="Fryer_1959",]

for (i in 1:sample_times) {
  
  Fryer_1959_random_samples[i,] <- as.matrix(sample_n(Fryer_1959, 1))
  
}


#for aquatic non-Ecopath networks sourced from publications 
#that provided only a single network each
single_aquatic_not_ecopath_food_webs <- metadata[metadata$author=="One_network_per_publication",]
single_aquatic_not_ecopath_food_webs <- single_aquatic_not_ecopath_food_webs[single_aquatic_not_ecopath_food_webs$Primary_type=="Aquatic",]
single_aquatic_not_ecopath_food_webs <- single_aquatic_not_ecopath_food_webs[single_aquatic_not_ecopath_food_webs$Ecopath=="Not_Ecopath_model",]

#storage for mean pairwise GCD-11 between networks that are aquatic non-Ecopath 
#from particular realization
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_networks")

#row counter for average_distance_per_realization
row_counter_for_average_distance_per_realization <- 1

for (i in 1:sample_times){
  
  #obtaining each networks from publications that were randomly sampled for a specific realization
  realization <- as.data.frame(rbind(single_aquatic_not_ecopath_food_webs,Closs_and_Lake_1994_random_samples[i,], Thompson_and_Townsend_2003_random_samples[i,],
                                     Alcorlo_et_al_2001_random_samples[i,], Menge_and_Sutherland_1976_random_samples[i,], Tavares_Cromar_and_Williams_1996_random_samples[i,],
                                     Parker_and_Huryn_2006_random_samples[i,], Thompson_and_Townsend_2004_random_samples[i,],Cohen_et_al_2003_random_samples[i,],
                                     Fryer_1959_random_samples[i,]))
  
  #list of rows and columns to delete which are not
  #of interest right now
  rows_columns_delete <- c()
  
  #looping through the GCD-11 matrix
  for (j in 1:ncol(ecological)) {
    
    #the name of the row we are looking through
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the network is
    while (found != TRUE) {
      
      #if we exceed the number of rows in the realization, this network should be removed
      if (rowNumber>nrow(realization)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
      #if we have not exceed the number of rows in the realization
      if (rowNumber<=nrow(realization)) {
        
        rowOfInterest <- as.data.frame(realization[rowNumber,])
        
        #if the row name in the matrix is the network we are looking for
        #we want to keep this network
        if (as.character(rowOfInterest$name) == network_name) {
          
          found = TRUE
          
          #otherwise we have not found the network in the metadata list
        } else {
          
          rowNumber = rowNumber + 1
          
        }
        
      }
      
      
    }
    
  }
  
  #deleting networks that are not the grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #realization number
  average_distance_per_realization[row_counter_for_average_distance_per_realization,1] <- i
  
  #mean pairwise GCD-11 between networks each from a different publication
  average_distance_per_realization[row_counter_for_average_distance_per_realization,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the realization
  average_distance_per_realization[row_counter_for_average_distance_per_realization,3] <- ncol(ecological_intermediate)
  
  row_counter_for_average_distance_per_realization = row_counter_for_average_distance_per_realization + 1
  
}

#value for mean pairwise GCD-11 between aquatic non-Ecopath networks
#with publication effect removed (S1 Appendix Section  S1.6.1)
round(mean(average_distance_per_realization[,2]),2)

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  networks that are aquatic and created using 
##  Ecopath (Appendix S1.6.2)
##
##
##
########################################################

##
##  Obtaining all 72 combinations of Ecopath networks
##  (i.e., removing for publication effect)
##

#Angelini_et_al_2006 Ecopath networks
Angelini_et_al_2006_storage <- matrix(,nrow=0,ncol(metadata))
colnames(Angelini_et_al_2006_storage) <- c(colnames(metadata))
Angelini_et_al_2006 <-  metadata[metadata$author=="Angelini_et_al_2006",]

for (i in 1:36) {
  
  Angelini_et_al_2006_storage <- rbind(Angelini_et_al_2006_storage,Angelini_et_al_2006)
  
}

#Stewart_&_Sprules_2011 Ecopath networks
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

#Angelini_et_al_2013 Ecopath networks
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

#Baeta_et_al_2011 Ecopath networks
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


#networks that are Ecopath models but from publications that produced only a single network
single_food_webs_ecopath <- metadata[metadata$author=="One_network_per_publication",]
single_food_webs_ecopath <- single_food_webs_ecopath[single_food_webs_ecopath$Ecopath=="Ecopath_model",]

#storage for all 72 possible combinations of mean pairwise GCD-11 between 
#networks that are Ecopath
average_distance_per_realization <- matrix(,nrow=72,ncol=3)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_networks")

#row counter for average_distance_per_realization
row_counter_for_average_distance_per_realization <- 1

for (i in 1:72){
  
  #obtaining each publications sampled network for a specific realization, and all networks
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
      
    #looping to find what number the network is
    while (found != TRUE) {
      
      #if we exceed the number of rows in the realization, this network should be removed
      if (rowNumber>nrow(realization)) {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
        found = TRUE
        
      }
      
      #if we have not exceed the number of rows in the realization
      if (rowNumber<=nrow(realization)) {
        
        rowOfInterest <- as.data.frame(realization[rowNumber,])
        
        #if the row name in the matrix is the network we are looking for
        #we want to keep this network
        if (as.character(rowOfInterest$name) == network_name) {
          
          found = TRUE
          
        #otherwise we have not found the network in the metadata list
        } else {
          
          rowNumber = rowNumber + 1
          
        }
        
      }

        
    }
      
  }
    
  #deleting networks that are not the grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
    
  #realization number
  average_distance_per_realization[row_counter_for_average_distance_per_realization,1] <- i
    
  #mean pairwise GCD-11 between networks that are Ecopath
  average_distance_per_realization[row_counter_for_average_distance_per_realization,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
    
  #number of networks that are Ecopath
  average_distance_per_realization[row_counter_for_average_distance_per_realization,3] <- ncol(ecological_intermediate)
    
  row_counter_for_average_distance_per_realization = row_counter_for_average_distance_per_realization + 1
  
}

#value for mean pairwise GCD-11 between Ecopath networks
#with publication effect removed (S1 Appendix Section  S1.6.2)
round(mean(average_distance_per_realization[,2]),2)
