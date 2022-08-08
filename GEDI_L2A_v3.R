#######   Description   ######
#' Author: Savannah Cooley
#' Description: This script runs a GEDI data processing function that filters the data following Potapov et al. (2021),
#' then spatially subsets the data over a user-defined area of interest and saves the data as .csv files.
#######   Install packages ######
library(data.table)
library(matrixStats)
library(dplyr)
library(ggplot2)
library(bit64) # bit64 to address Warning message: In require_bit64_if_needed(x) : Some columns are type 'integer64' but package bit64 is not installed. Those columns will print as strange looking floating point data. There is no need to reload the data. Simply install.packages('bit64') to obtain the integer64 print method and print the data again.
library(htm2txt) # for file download part of script
library(magrittr) # for file download part of script
library(rhdf5) # for reading in .hdf5 files
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("rhdf5")

#######   Prep work: define AOI and input variables    ######
root<-"/Volumes/Cooley_Sgate/SERVIR_Am/Data/Peru/L2A"
download_dir<-paste(root,"/","Orig_hdf5/Hdfs_round2",sep="")
outPath<-paste(root, "/_Subset_QA/r2_processed_Peru_Am",sep="")

ls_files<- list.files(download_dir, full.names = F, pattern = "*.h5", recursive = F)
ls_files_long<- list.files(download_dir, full.names = T, pattern = "*.h5", recursive = F)
done=list.files(outPath, pattern = "*.csv")

# Only use data from full power beams (i.e. not split beams: "BEAM0000","BEAM0001","BEAM0010","BEAM0011")
beam_list<-as.vector(c("BEAM0101","BEAM0110","BEAM1000","BEAM1011")) 

# AOI boundary -- cover all Peruvian Amazon # Already clipped to AOI by EarthExplorer
#x_max= -68.5
#x_min= -78.3
#y_max= -1.6
#y_min= -14.2

#######   Define functions   ######

outliers_to_NA <- function(v) {

lower_bound <- quantile(v, 0.025, na.rm=T)
upper_bound <- quantile(v, 0.975, na.rm=T)
#lower_bound <- boxplot.stats(v)$conf[1] # Alternate approach
#upper_bound <- boxplot.stats(v)$conf[2]

# Remove the farthest outlier among the 6 algorithms
v[v<lower_bound]<-NA
v[v>upper_bound]<-NA
return(v)
}

Pre_process_GEDI_L2A <- function(f, f_short, beams, outdir) {
  setwd(outdir)
  dt_all <- data.table()
  #yyyyddd<-substr(f_short, 10,16) 
  #hhmmss<-substr(f_short, 17,22) 
  tryCatch({ # This is the 'try' part
    for (beam_num in beams){ # iterate over all the beam numbers
      #######   Read in data     ######
      dt <- data.table() # reset
      
      # For QA filtering
      dt$L2A_QA <- h5read(f, paste("/", beam_num,"/rx_assess/quality_flag", sep="")) # Flag simpilfying selection of most useful data. # A quality_flag value of 1 indicates the laser shot meets criteria based on energy, sensitivity, amplitude, real-time surface tracking quality and difference to a DEM. 
      dt$L2A_QA <-as.integer(dt$L2A_QA) 
      #dt$surf_flag <- h5read(f, paste("/", beam_num,"/surface_flag", sep="")) # Redundant if using quality_flag. Indicates elev_lowestmode is within 300m of DEM or MSS
      #dt$surf_flag <- as.integer(dt$surf_flag)
      
      if(length(dt[L2A_QA==1])>0){ # Ensure some quality, non-noise data exist 
        
        #######   Retrieve most important information    ######
        dt$shot_number <- h5read(f, paste("/", beam_num,"/shot_number", sep=""), bit64conversion='bit64') # Shot ID
        dt$delta_time <- h5read(f, paste("/",beam_num,"/delta_time", sep="")) # Time delta since Jan 1 00:00 2018.
        dt$lat <- h5read(f, paste("/",beam_num,"/lat_lowestmode", sep="")) # Latitude of center of lowest mode
        dt$lon <- h5read(f, paste("/",beam_num,"/lon_lowestmode", sep="")) # Longitude of center of lowest mode
        #dt$lat_a1 <- h5read(f, paste("/",beam_num,"/geolocation/lat_lowestmode_a1", sep="")) # lat for a1 algorithm
        #dt$lon_a1 <- h5read(f, paste("/",beam_num,"/geolocation/lon_lowestmode_a1", sep="")) # lon for a1 algorithm
        #dt$lat_h <- h5read(f, paste("/",beam_num,"/geolocation/lat_highestreturn", sep="")) # lat for highest detected return
        #dt$lon_h <- h5read(f, paste("/",beam_num,"/geolocation/lat_highestreturn", sep="")) # lon for highest detected return
        
        dt$solar_elevation<- h5read(f, paste("/",beam_num,"/solar_elevation", sep="")) # Solar_elevation for filtering
        
        #rx_modeamps <- h5read(f, paste("/",beam_num,"/rx_processing_a1/rx_modeamps", sep="")) # Amplitudes of each detected mode within waveform
        dt$rx_energy  <- h5read(f, paste("/",beam_num,"/rx_assess/rx_energy", sep="")) # Integrated energy in RX waveform after subtracting the noise mean.
        dt$rx_gamplitude <- h5read(f, paste("/",beam_num,"/rx_1gaussfit/rx_gamplitude", sep="")) # Amplitude of single gaussian fit to the rxwaveform
        dt$rx_gamplitude_error <- h5read(f, paste("/",beam_num,"/rx_1gaussfit/rx_gamplitude_error", sep="")) # Error in ampltiude estimate for single gaussian fit to the rxwaveform
        
        # Elevation
        dt$DEM <- h5read(f, paste("/",beam_num,"/digital_elevation_model", sep="")) # TanDEM-X elevation at GEDI footprint location 
        dt$elev_lowestmode <- h5read(f, paste("/",beam_num,"/elev_lowestmode", sep="")) # Elevation of lowest return from selected algorithm
        #dt$elev_highestreturn <- h5read(f, paste("/",beam_num,"/elev_highestreturn", sep="")) # Elevation of highest return
        
        #######   Retrieve relative height metrics for selected best algorithm   ######
        rh_m <- h5read(f, paste("/",beam_num,"/rh", sep="")) # Relative height metrics at 1% interval (in cm)
        rh_m <- t(rh_m) # transpose matrix so each of the 100 RH percentiles is a column
        rh_m <- data.table(rh_m) # typecast to data.table
        rh_m <- subset(x=rh_m, select=c("V10", "V20", "V30", "V40", "V50", "V60", "V70", "V80", "V90", "V100","V25", "V75", "V95", "V98"))
        
        dt<-cbind(dt, rh_m)
      
        #######   Noise removal    ######
        
        # Elevation
        elev_a1 <- h5read(f, paste("/",beam_num,"/geolocation/elev_lowestmode_a1", sep="")) # Elevation of lowest return for a1 algorithm
        elev_a2 <- h5read(f, paste("/",beam_num,"/geolocation/elev_lowestmode_a2", sep="")) # Elevation of lowest return for a2 algorithm
        elev_a3 <- h5read(f, paste("/",beam_num,"/geolocation/elev_lowestmode_a3", sep="")) # Elevation of lowest return for a3 algorithm
        elev_a4 <- h5read(f, paste("/",beam_num,"/geolocation/elev_lowestmode_a4", sep="")) # Elevation of lowest return for a4 algorithm
        elev_a5 <- h5read(f, paste("/",beam_num,"/geolocation/elev_lowestmode_a5", sep="")) # Elevation of lowest return for a5 algorithm
        elev_a6 <- h5read(f, paste("/",beam_num,"/geolocation/elev_lowestmode_a6", sep="")) # Elevation of lowest return for a6 algorithm
        
        # Remove outlier elevation estimates among the 6 algorithms
        t<-cbind(elev_a1,elev_a2,elev_a3,elev_a4,elev_a5,elev_a6)
        
        for(k in 1:dim(t)[1]){
          # Update row
          t[k,]<-outliers_to_NA(t[k,])
        }
        
        # Calculate the range among the remaining non-NA elevation values for each row
        dt$elev_min<- rowMins(t, na.rm = T)
        dt$elev_max<- rowMaxs(t, na.rm = T)
        dt$elev_range<- dt$elev_max-dt$elev_min # For filtering later
        
        # Subset data based on range in elevation estimates
        dt<- subset(dt, dt$elev_range<=2) # Filter such that the range of predicted ground elevations among the non-outlier rows from the elevation algorithms is ≤2m
        
        # Subset all data based on sensitivity (only include data with strong enough signal to penetrate dense canopy)
        dt$sensitivity <- h5read(f, paste("/",beam_num,"/sensitivity", sep="")) 
        dt<- subset(dt, sensitivity>=0.9)
        
        # Additional filtering 
        #dt<- subset(dt, solar_elevation < 0) # Limit the background noise effects of reflected solar radiation)
        #dt<- subset(dt, surf_flag>0) # only keep quality data 
        
        dt<- subset(dt, V98>0) # Quality check
        
        #dt<- subset(dt, L2A_QA>0) # A quality_flag value of 1 indicates the laser shot meets criteria based on energy, sensitivity, amplitude, real-time surface tracking quality and difference to a DEM. 
        #dt$L2A_QA<-NULL # Drop column since L2A_QA value is 1 for all remaining data
        
        dt$beam<- beam_num
        dt_all<-rbind(dt, dt_all)
      } # end QA if
      
    } # end j loop 
    #######   Write data to file     ######
    #dt_all$yyyyddd<-yyyyddd
    #dt_all$hhmmss<-hhmmss
    
    # Write clipped data to file
    fname<-paste("sub_",substr(f_short,1,nchar(f_short)-3), ".csv", sep="")
    fwrite(dt_all, fname)
    message(paste("Processed file:", f))
  },
  error=function(cond) {
    message(paste(i, "--> Issue opening the data:", f_short))
    message("Here's the original error message:")
    message(cond)
  },
  warning=function(cond) {
    message(paste("Data caused a warning:", f_short))
    message("Here's the original warning message:")
    message(cond)
  }) # end Try/Catch  
}

#######   Run processing function for GEDI L2A data     ######
setwd(outPath) 

for(i in 964:length(ls_files)){ 
  fn<-substr(ls_files[i],1,nchar(ls_files[i])-3) 
  
  if(any(done %like% fn)){
    print(paste("Already processed: ", ls_files[i]))
  } else {
    
    # Read in, filter, spatially subset & write GEDI L2A data
    Pre_process_GEDI_L2A(ls_files_long[i], ls_files[i], beam_list, outPath)
  }
  
  #file.remove(ls_files_long[i])
  #print(paste("Deleted file: ", ls_files[i]))
}