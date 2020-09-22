##
## This function requires two inputs a filepath indicating your input data in the format of the example data provided in this
## repository, and a list of the species you want to include from those recorded in the data file. The species list should be
## written using the english names used in your data file.
## This function is based on code provided in the book by Kéry M & Royle JA (2016) (Applied Hierarchical Modeling in Ecology: 
## Analysis of distribution, abundance and species richness in R and BUGS: Volume 1).
##
## The code returns a list containing different data you will need for the analysis. We used date and temperature as covariates
## in the analysis, so these are also returned in the list for use in the MSOM. There are parameters in here that will need 
## adjusting for use in other datasets.
##


dataprep <- function(filepath, specieslist) {
  # Load packages
  require(ggplot2)
  require(tidyverse)
  require(data.table)
  
  # Read in data
  df <- read.csv(filepath, row.names=1, header = T, sep = ",", dec = ".") #Read dataset
  dfilter <- df %>% dplyr:: filter(engname %in% specieslist) # Filter the dataset by your species list
  dfilter <- droplevels(dfilter) # Remove the factor levels from the species no longer included in the dataset
  dfilter <- dfilter[order(dfilter$engname),] # Reorder alphabetically
  
  # Create various species lists
  species.list <- levels(dfilter$engname)  # create alphabetic species list
  spec.name.list <- tapply(dfilter$specid, dfilter$engname, mean) #species ID
  spec.id.list <- unique(dfilter$specid) #ID list
  ordered.spec.name.list <- spec.name.list[order(spec.name.list)]
  
  COUNTS <- cbind(dfilter$count1, dfilter$count2, dfilter$count3, dfilter$count4, dfilter$count5)  # Put count columns into a new dataset
  DET <- COUNTS
  DET[DET > 1] <- 1   # Turn the counts into detection/ nondetection data
  
  # Put detection data into 3D array
  nsite <- length(unique(dfilter$site))  # Number of sites in dataset
  nrep <- 5     # Number of replicate surveys per season
  nspec <- length(species.list) # How many species occur in the dataset
  
  # Create a 3D species matrix
  y <- array(NA, dim=c(nsite, nrep, nspec))   
  for(i in 1:nspec){
    y[,,i] <- DET[((i-1)*nsite + 1):(i*nsite),]
  }
  dimnames(y) <- list (NULL, NULL, names(spec.name.list))
  class(y)<-"numeric"
  y[is.na(y)]<-0
  (nspec <- dim(y)[3]) # Redefine nspec as the observed number of species
  
  # Standardise the date and temperature data
  temp <- as.matrix(dfilter %>% select(maxtemp1, maxtemp2, maxtemp3, maxtemp4, maxtemp5))
  stdtemp <- (temp - mean(temp)) / sd(temp)
  
  date <- as.matrix(dfilter %>% select(date1, date2, date3, date4, date5))
  stddate <- (date - mean(date)) / sd(date)
  
  
  # Output some simple data checks using the 3D species matrix
  tmp <- apply(y, c(1,3), max, na.rm=T) 
  
  sort(obs.occ <- apply(tmp, 2, sum, na.rm=T)) # Observed number of occupied sites
  ## Plot species 'occurrence frequency' distribution
  occurfreq <- plot(sort(obs.occ), xlab="Species number", ylab="Number of sites*days with detections") 

  ## Get observed number of species per site
  tmp <- apply(y, c(1,3), max, na.rm=T)
  tmp[tmp == "-Inf"] <- NA
  sort(C <- apply(tmp, 1, sum)) #Compute and print sorted species counts
  
  ## Plot the observed species richness
  obssprich <- plot(table(C), xlim=c(0,10), xlab="Observed Number of Species", ylab = "Number of Sites", frame=F)
  abline(v=mean(C, na.rm=T), col="blue", lwd=3)
  
  
  ## Condense data to look at detections of each species by site
  ysum <- as.data.frame(apply(y, c(1,3), sum))
  ysum <- setDT(ysum, keep.rownames="site")
  df_x <- pivot_longer(ysum, cols = -1, names_to = "Species", values_to = "Detections")
  df_x$site <- as.numeric(df_x$site)
  
  spdections <- ggplot(df_x, aes(x=site, y=Species)) +
    geom_tile(aes(fill = Detections), colour = "white") +
    scale_fill_gradient(low = "white", high = "black") +
    labs(x = "Site Number", y = "Species") +
    theme(panel.background = element_rect("white"),
          panel.border = element_rect(color="black", fill="transparent"),
          axis.text=element_text(size=9),
          axis.title=element_text(size=10),
          legend.title=element_text(size=10),
          legend.text=element_text(size=9))
  
  return(list(y=y, ysum=ysum, nsite=nsite, nrep=nrep, nspec=nspec, species.list=species.list, spec.name.list=spec.name.list,
              spec.id.list=spec.id.list, ordered.spec.name.list=ordered.spec.name.list, C=C,obs.occ=obs.occ,
              stddate=stddate, stdtemp=stdtemp, date=date, temp=temp, plotspdet = spdections))
}
