########################################################
##
##
##
##  For recreating the results of Brimacombe et al. 
##  (2024): On the nature of structure in collections  
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
##
##
##  Used to evaluate the mean pairwise GCD-11 between 
##  food webs from the same publication grouping as
##  shown in parts of table 2 and table S6
##
##
##
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

#mean pairwise GCD-11 (in table S6)
each_average_distance <- as.data.frame(each_average_distance)
each_average_distance$mean_distance <- as.numeric(each_average_distance$mean_distance)
each_average_distance$number_of_networks <- as.numeric(each_average_distance$number_of_networks)
each_average_distance

#weighted mean pairwise GCD-11 of food webs sourced from the same publication that produced multiple networks (in table 2)
webs_from_pubs_with_multiple_networks <- each_average_distance[!(each_average_distance$network_number=="One_network_per_publication"),]
round(weighted.mean(webs_from_pubs_with_multiple_networks$mean_distance,webs_from_pubs_with_multiple_networks$number_of_networks),2)

#mean pairwise GCD-11 of food webs sourced from publications that only produced a single network (in table 2)
webs_from_pubs_with_single_network <- each_average_distance[each_average_distance$network_number=="One_network_per_publication",]
webs_from_pubs_with_single_network

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs sourced from publications that produced 
##  multiple networks, before 1990 and after 1990 
##  (figure 4 and parts of table 2)
##
##
##
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
        #food web in the pairwise GCD-11 matrix)
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

#figure 4
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

#weighted mean of mean pairwise GCD-11 from publications before or during 1990s (parts of table 2)
each_average_distance_for_publication_before_1990 <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year<2000,]
round(weighted.mean(each_average_distance_for_publication_before_1990$mean_distance,each_average_distance_for_publication_before_1990$number_of_networks),2)

#weighted mean of mean pairwise GCD-11 from publications after 1990s (parts of table 2)
each_average_distance_for_publication_after_1990 <- each_average_distance_for_publication[each_average_distance_for_publication$publication_year>=2000,]
round(weighted.mean(each_average_distance_for_publication_after_1990$mean_distance,each_average_distance_for_publication_after_1990$number_of_networks),2)

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs from the same ecological system (parts
##  of table 1)
##
##
##
########################################################

#storage for mean pairwise GCD-11 between food webs from the same ecosystem
each_average_distance <- matrix(,nrow=(length(levels(as.factor(metadata$Primary_type)))),ncol=3)
colnames(each_average_distance) <- c("Ecosystem_type","mean_distance","number_of_networks")

for (i in 1:nrow(each_average_distance)) {
  
  each_average_distance[i,1] <- levels(as.factor(metadata$Primary_type))[i]
  
}

#looping through the pairwise GCD-11 matrix to find food webs associated 
#with specific ecosystems
for (i in 1:nrow(each_average_distance)) {
  
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
        
        #if the row is not an ecosystem food web we want, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the GCD-11 matrix)
        if (rowOfInterest$Primary_type != each_average_distance[i,1]) {
          
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
  
  #mean pairwise GCD-11 between specific ecosystem's food webs
  each_average_distance[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific ecosystem grouping of interest
  each_average_distance[i,3] <- ncol(ecological_intermediate)
  
}

each_average_distance <- as.data.frame(each_average_distance)
each_average_distance$mean_distance <- as.numeric(each_average_distance$mean_distance)
each_average_distance$number_of_networks <- as.numeric(each_average_distance$number_of_networks)

#parts of table 1
each_average_distance

##
##  Mean pairwise GCD-11 between "terrestrial" food webs, 
##  with "Digel et al. (2014)" food webs removed (parts
##  of table 1)
##

#storage for mean pairwise GCD-11 between food webs from the same ecosystem
each_average_distance <- matrix(,nrow=1,ncol=3)
colnames(each_average_distance) <- c("Ecosystem_type","mean_distance","number_of_networks")
each_average_distance[1,1] <- "Terrestrial"
  

#looping through the pairwise GCD-11 matrix to find food webs associated 
#with specific ecosystems
for (i in 1:nrow(each_average_distance)) {
  
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
        
        #if the row is not an ecosystem food web we want, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the pairwise GCD-11 matrix)
        if (rowOfInterest$Primary_type != each_average_distance[i,1]) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        #if the row is Digel_et_al_2014, remove
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
  
  #mean pairwise GCD-11 between specific ecosystem's food webs
  each_average_distance[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific ecosystem grouping of interest
  each_average_distance[i,3] <- ncol(ecological_intermediate)
  
}

#parts of table 1
each_average_distance

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs from different ecological systems (parts 
##  of table 1)
##
##
##
########################################################

##
##  Mean pairwise GCD-11 between "aquatic" and 
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
    
    #if the row name in the matrix is the food web we 
    #are looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is a "terrestrial" food web, we will need 
      #to delete it
      if (rowOfInterest$Primary_type == "Terrestrial") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Aquatic") {
        
        columns_delete <- c(columns_delete,j)
        
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
ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]

#mean pairwise GCD-11 between "aquatic" and "aquatic and terrestrial" food webs (parts of table 1)
round(mean(ecological_intermediate),2)

##
##  Mean pairwise GCD-11 between "aquatic" and 
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
    
    #if the row name in the matrix is the food web we are 
    #looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is an "aquatic and terrestrial" food web, we will 
      #need to delete it
      if (rowOfInterest$Primary_type == "Aquatic and terrestrial") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Aquatic") {
        
        columns_delete <- c(columns_delete,j)
        
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
ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]

#mean pairwise GCD-11 between "aquatic" and "terrestrial" food webs (parts of table 1)
round(mean(ecological_intermediate),2)

##
##  Mean pairwise GCD-11 between "aquatic and terrestrial" 
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
    
    #if the row name in the matrix is the food web we are 
    #looking for
    if (as.character(rowOfInterest$name) == network_name) {
      
      #if the row is an "aquatic" food web, we will need 
      #to delete it
      if (rowOfInterest$Primary_type == "Aquatic") {
        
        rows_columns_delete <- c(rows_columns_delete,j)
        
      } else if (rowOfInterest$Primary_type == "Aquatic and terrestrial") {
        
        columns_delete <- c(columns_delete,j)
        
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
ecological_intermediate <- ecological[-c(rows_columns_delete,rows_delete),-c(rows_columns_delete,columns_delete)]
 
#mean pairwise GCD-11 between "aquatic and terrestrial" and "terrestrial" food webs (parts of table 1)
round(mean(ecological_intermediate),2)

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs that are "lake", "marine", "river", and 
##  "stream" (Appendix S1: Section S1.5)
##
##
##
########################################################

##
##  Mean pairwise GCD-11 between "lake" food webs 
##  with publication removed
##

#number of times you want to randomly sample/number of realizations
sample_times <- 200

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
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
  
}

#parts of table S2
round(mean(average_distance_per_realization[,2]),2)

##
##  Mean pairwise GCD-11 between "marine" food webs 
##  with publication removed
##

#number of times you want to randomly sample/number of realizations
sample_times <- 200

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
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of networks that are of the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
}

#parts of table S2
round(mean(average_distance_per_realization[,2]),2)

##
##  Mean pairwise GCD-11 between "river" food webs
##  with publication removed
##

#number of times you want to randomly sample/number of realizations
sample_times <- 200

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
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
}

#parts of table S2
round(mean(average_distance_per_realization[,2]),2)

##
##  Mean pairwise GCD-11 between "stream" food webs
##  with publication removed
##

#number of times you want to randomly sample/number of realizations
sample_times <- 200

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
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
}

#parts of table S2
round(mean(average_distance_per_realization[,2]),2)

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs from different ecological systems (parts 
##  of table S2) with publication removed
##
##
##
########################################################

##
##  Mean pairwise GCD-11 between "lake" 
##  and "marine" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "lake" food webs are columns, and "marine" food webs are rows
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "lake" and in the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "marine" and in the realization
  average_distance_per_realization[i,4] <- nrow(ecological_intermediate)
  
}

#mean pairwise GCD-11 between "lake" and "marine" food webs with publication effect removed (part of table S2)
round(mean(average_distance_per_realization[,2]),2)

##
##  Mean pairwise GCD-11 between "lake" 
##  and "river" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "lake" food webs are columns, and "river" food webs are rows
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "lake" and in the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "river" and in the realization
  average_distance_per_realization[i,4] <- nrow(ecological_intermediate)
  
}

#mean pairwise GCD-11 between "lake" and "river" food webs with publication effect removed (part of table S2)
round(mean(average_distance_per_realization[,2]),2)

##
##  Mean pairwise GCD-11 between "lake" 
##  and "stream" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "lake" food webs are columns, and "stream" food webs are rows
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "lake" and in the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "stream" and in the realization
  average_distance_per_realization[i,4] <- nrow(ecological_intermediate)
  
}

#mean pairwise GCD-11 between "lake" and "stream" food webs with publication effect removed (part of table S2)
round(mean(average_distance_per_realization[,2]),2)

##
##  Mean pairwise GCD-11 between "marine" 
##  and "river" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "marine" food webs are columns, and "river" food webs are rows
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "marine" and in the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "river" and in the realization
  average_distance_per_realization[i,4] <- nrow(ecological_intermediate)
  
}

#mean pairwise GCD-11 between "marine" and "river" food webs with publication effect removed (part of table S2)
round(mean(average_distance_per_realization[,2]),2)

##
##  Mean pairwise GCD-11 between "marine" 
##  and "stream" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "marine" food webs are columns, and "stream" food webs are rows
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "marine" and in the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "stream" and in the realization
  average_distance_per_realization[i,4] <- nrow(ecological_intermediate)
  
}

#mean pairwise GCD-11 between "marine" and "stream" food webs with publication effect removed (part of table S2)
round(mean(average_distance_per_realization[,2]),2)

##
##  Mean pairwise GCD-11 between "river" 
##  and "stream" food webs with publication
##  removed
##

#storage for mean pairwise GCD-11 between food webs from particular realization
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=4)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_column_networks","number_of_row_networks")

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
  average_distance_per_realization[i,1] <- i
  
  #mean pairwise GCD-11 between food webs 
  #where for ecological_intermediate: the "river" food webs are columns, and "stream" food webs are rows
  average_distance_per_realization[i,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are both "river" and in the realization
  average_distance_per_realization[i,3] <- ncol(ecological_intermediate)
  
  #number of food webs that are both "stream" and in the realization
  average_distance_per_realization[i,4] <- nrow(ecological_intermediate)
  
}

#mean pairwise GCD-11 between "river" and "stream" food webs with publication effect removed (part of table S2)
round(mean(average_distance_per_realization[,2]),2)

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs sourced from a single publication 
##  (parts of Appendix S1: Section S1.7)
##
##
##
########################################################

#storage for mean and standard deviation in pairwise GCD-11
average_distance_per_realization_network_size<- matrix(,nrow=4,ncol=5)
colnames(average_distance_per_realization_network_size) <- c("network_size_quantile_number","mean_network_size","SD_network_size","mean_distance","number_of_networks")

#obtaining all food webs sourced from a publication that produced only a single network
realization <- metadata[metadata$author=="One_network_per_publication",]

realization$number_of_nodes <- as.numeric(realization$number_of_nodes)

#identifying the quartile that each food web belongs to (based on number of nodes)
Quartile<-cut(realization$number_of_nodes,quantile(realization$number_of_nodes,type=8),include.lowest=TRUE,labels=FALSE)

#adding quartile information as its own column
realization <- cbind(realization,Quartile)

#obtaining only quartile and food web name
network_name_and_quartile <- as.data.frame(cbind(realization$name,realization$Quartile))

#naming the columns
colnames(network_name_and_quartile) <- c("name","network_size_quartile")

network_name_and_quartile$network_size_quartile <- as.numeric(network_name_and_quartile$network_size_quartile)

all_networks_and_quartiles <- merge(metadata,network_name_and_quartile,by="name",all=T)

#looping through the pairwise GCD-11 matrix to find food webs associated 
#with specific author/publication (here "4" is the number of quartiles)
for (k in 1:4) {
  
  #list of rows and columns to delete which are not of the 
  #author/publication grouping of interest right now
  rows_columns_delete <- c()
  
  #looping through the GCD-11 matrix
  for (j in 1:ncol(ecological)) {
    
    #the name of the row we are looking through
    network_name <- row.names(ecological)[j]
    
    found = FALSE
    rowNumber = 1
    
    #looping to find what number the food web is
    while (found != TRUE) {
      
      rowOfInterest <- as.data.frame(all_networks_and_quartiles[rowNumber,])
      
      #if the row name in the matrix is the food web we are looking for
      if (as.character(rowOfInterest$name) == network_name) {
        
        #if the row is not an author's food web we want, or is NA, we will need 
        #to delete it (here, we just keep the location of this food web
        #in the GCD-11 matrix)
        if (is.na(rowOfInterest$network_size_quartile)) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        } else if (rowOfInterest$network_size_quartile != k) {
          
          rows_columns_delete <- c(rows_columns_delete,j)
          
        }
        
        found = TRUE
        
      #otherwise we have not found the food web in the metadata list
      } else {
        
        rowNumber = rowNumber + 1
        
      }
      
    }
    
  }
  
  #deleting food webs that are not the author/publication grouping of interest now
  #(pairwise GCD-11 matrix is symmetric so we can delete the same 
  #rows and columns)
  ecological_intermediate <- ecological[-rows_columns_delete,-rows_columns_delete]
  
  #quartile number
  average_distance_per_realization_network_size[k,1] <- k
  
  #mean food web size between specific author/publication's food webs
  average_distance_per_realization_network_size[k,2] <- round(mean(realization[realization$Quartile==k,]$number_of_nodes),2)
  
  #standard deviation in food web size between specific author/publication's food webs
  average_distance_per_realization_network_size[k,3] <- round(sd(realization[realization$Quartile==k,]$number_of_nodes),2)
  
  #mean pairwise GCD-11 between specific author/publication's food webs
  average_distance_per_realization_network_size[k,4] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the specific author/publication grouping of interest
  average_distance_per_realization_network_size[k,5] <- ncol(ecological_intermediate)
  
}

#table S5
average_distance_per_realization_network_size

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs that are "aquatic" and not created using 
##  Ecopath (Appendix S1: Section S1.8.1)
##
##
##
########################################################

##
##  Randomly obtaining one food web per publication grouping
##

#number of times you want to randomly sample/number of realizations
sample_times <- 200

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

#for "aquatic" non-Ecopath food webs sourced from publications 
#that provided only a single network
single_aquatic_not_ecopath_food_webs <- metadata[metadata$author=="One_network_per_publication",]
single_aquatic_not_ecopath_food_webs <- single_aquatic_not_ecopath_food_webs[single_aquatic_not_ecopath_food_webs$Primary_type=="Aquatic",]
single_aquatic_not_ecopath_food_webs <- single_aquatic_not_ecopath_food_webs[single_aquatic_not_ecopath_food_webs$Ecopath=="Not_Ecopath_model",]

#storage for mean pairwise GCD-11 between food webs that are "aquatic" non-Ecopath 
#from particular realization
average_distance_per_realization <- matrix(,nrow=sample_times,ncol=3)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_networks")

#row counter for average_distance_per_realization
row_counter_for_average_distance_per_realization <- 1

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
  average_distance_per_realization[row_counter_for_average_distance_per_realization,1] <- i
  
  #mean pairwise GCD-11 between food webs each from a different publication
  average_distance_per_realization[row_counter_for_average_distance_per_realization,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
  
  #number of food webs that are of the realization
  average_distance_per_realization[row_counter_for_average_distance_per_realization,3] <- ncol(ecological_intermediate)
  
  row_counter_for_average_distance_per_realization = row_counter_for_average_distance_per_realization + 1
  
}

#value for mean pairwise GCD-11 between "aquatic" non-Ecopath food webs
#with publication effect removed (Appendix S1: Section  S1.8.1)
round(mean(average_distance_per_realization[,2]),2)

########################################################
##
##
##
##  Used to evaluate the mean pairwise GCD-11 between
##  food webs that are "aquatic" and created using 
##  Ecopath (Appendix S1: Section S1.8.2)
##
##
##
########################################################

##
##  Obtaining all 72 combinations of Ecopath "aquatic" food webs
##  (i.e., removing for publication effect)
##

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
average_distance_per_realization <- matrix(,nrow=72,ncol=3)
colnames(average_distance_per_realization) <- c("realization","mean_distance","number_of_networks")

#row counter for average_distance_per_realization
row_counter_for_average_distance_per_realization <- 1

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
  average_distance_per_realization[row_counter_for_average_distance_per_realization,1] <- i
    
  #mean pairwise GCD-11 between food webs that are Ecopath "aquatic"
  average_distance_per_realization[row_counter_for_average_distance_per_realization,2] <- round(mean(ecological_intermediate[lower.tri(ecological_intermediate)]),2)
    
  #number of food webs that are Ecopath "aquatic"
  average_distance_per_realization[row_counter_for_average_distance_per_realization,3] <- ncol(ecological_intermediate)
    
  row_counter_for_average_distance_per_realization = row_counter_for_average_distance_per_realization + 1
  
}

#value for mean pairwise GCD-11 between Ecopath "aquatic" food webs
#with publication effect removed (Appendix S1: Section  S1.8.2)
round(mean(average_distance_per_realization[,2]),2)

########################################################
##
##
##
##  Used to evaluate the pairwise GCD-11 between
##  food webs (Appendix S1: Section S1.4)
##
##
##
########################################################

#reading in the pairwise GCD-11 matrix
ecological <- as.matrix(read.csv("gcd11.csv",row.names = 1))
ecological_intermediate <- ecological[lower.tri(ecological)]

storage_for_entries <- as.matrix(matrix(0,nrow=3,ncol=15))
colnames(storage_for_entries) <- c("0.0<=x<=0.5","0.5<x<=1.0","1.0<x<=1.5","1.5<x<=2.0",
                                   "2.0<x<=2.5","2.5<x<=3.0","3.0<x<=3.5","3.5<x<=4.0",
                                   "4.0<x<=4.5","4.5<x<=5.0","5.0<x<=5.5","5.5<x<=6.0",
                                   "6.0<x<=6.5","6.5<x<=7.0","7.0<x<=7.5")
rownames(storage_for_entries) <- c("Only between food webs from publications that produced a single network (n=3403)",
                                   "Only between food webs from publications that produced multiple networks (n=2487)",
                                   "Between all other possible food web pairwise comparisons (n=31511)")

row_of_storage <- 3

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
# All publications that provided multiple networks
# 

#storage for mean pairwise GCD-11 between food webs from the same publication grouping
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

#
# For plotting figure S7
# 

storage_for_entries_2 <- storage_for_entries

#so that we are not double counting
storage_for_entries_2[3,] <- storage_for_entries_2[3,] - storage_for_entries_2[1,] - storage_for_entries_2[2,]

intermediate <- as_tibble(as.data.frame(as.table(storage_for_entries_2)))
plotting_data <- as.data.frame(intermediate)

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

########################################################
##
##
##
##  Used to evaluate all pairwise GCD-11s based on
##  absolute difference in food web size and standard 
##  deviation in food web size (figure S10)
##
##
##
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
  
  #looping through the column of the GCD-11 matrix
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

#absolute food web size and pairwise GCD-11
summary(lm(distances$absolute_network_size~distances$pairwise_distance))

#standard deviation of food web size and pairwise GCD-11
summary(lm(distances$SD_network_size~distances$pairwise_distance))

#figure S10
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

#figure S10
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

