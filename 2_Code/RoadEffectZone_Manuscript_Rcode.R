# R code for Manuscript entitled: 
  # Influence of traffic volume on mammal beta diversity within the road effect zone
# Manuscript authors: 
  # Thomas J. Yamashita*
  # David B. Wester
  # Zachary M. Wardle
  # Daniel G. Scognamillo
  # Landon R. Schofield 
  # Michael E. Tewes
  # John H. Young Jr.
  # Jason V. Lombardi

# Primary code writer: 
  # Thomas J. Yamashita


#####################################################################################################

# Install packages ####
install.packages("devtools")
devtools::install_github("tomyamashita/cameraTrapping")

#####################################################################################################

#####################################################################################################################
############################################## Part 1: Loading the Data ############################################# 
#####################################################################################################################

rm(list = ls()); gc()

# Step 0: Get a list of the files in the data folder ####
files <- list.files(path = file.path(getwd(), "Data"), full.names = T)
files


# Step 1: Load the raw data ####
## DataOrganize Files
### Road Effect Zone Study
AP1_REZ <- do.call(rbind, pbapply::pblapply(files[grep("DataOrganize_REZ", files)], read.table, header = T))
### Ranch-Level Comparison
AP1_Ranch <- do.call(rbind, pbapply::pblapply(files[grep("DataOrganize_Grid", files)], read.table, header = T))

## Timelapse Files
timelapse <- pbapply::pblapply(files[grep("timelapse_out_", files)], read.csv, header = T)
names(timelapse) <- sub(".csv", "", sub("timelapse_out_", "", basename(files[grep("timelapse", files)])))
### Road Effect Zone Study
timelapse_REZ <- timelapse[!grepl("Grid", names(timelapse))]
### If using timelapse files instead of dataorganize files, create AP1
if(!exists("AP1_REZ")){AP1_REZ <- do.call(rbind, pbapply::pblapply(timelapse_REZ, cameraTrapping::doTimelapse))}
### Ranch-Level Comparison
timelapse_Ranch <- timelapse[grepl("Grid", names(timelapse))]

## Environmental Data
envdata_REZ <- openxlsx::read.xlsx(files[grep("EnvData_REZ", files)])
envdata_Ranch <- openxlsx::read.xlsx(files[grep("EnvData_Ranch", files)])

## CT Tables
### Load the CT Tables
cttable_REZ <- openxlsx::read.xlsx(files[grep("CTtable_REZ", files)], sheet = 1, detectDates = T)
cttable_REZ_all <- openxlsx::read.xlsx(files[grep("CTtable_REZ", files)], sheet = 2, detectDates = T)
cttable_Ranch <- openxlsx::read.xlsx(files[grep("CTtable_Ranch", files)], sheet = 1, detectDates = T)
cttable_Ranch_all <- openxlsx::read.xlsx(files[grep("CTtable_Ranch", files)], sheet = 2, detectDates = T)
### We may need to convert cttable dates into characters, because camtrapR is annoying and keeps changing things...
cttable_REZ <- cameraTrapping::ctDates(cttable_REZ)
cttable_REZ_all <- cameraTrapping::ctDates(cttable_REZ_all)
cttable_Ranch <- cameraTrapping::ctDates(cttable_Ranch)
cttable_Ranch_all <- cameraTrapping::ctDates(cttable_Ranch_all)
### Need to remove some non-relevant data from cttable_X_all because it gets used to ensure no issues with half-days
cttable_REZ_all <- cttable_REZ_all[cttable_REZ_all$Location != "FM1847",]
cttable_Ranch_all <- cttable_Ranch_all[cttable_Ranch_all$Ranch != "Yturria",]

## Species List
species_list <- openxlsx::read.xlsx(files[grep("Species_List", files)])


# Step 2: Check image and trap rates ####
## Road Effect Zone
### Use Timelapse files instead of dataorganize files because easier to access as separate items
effort_Image_REZ <- cameraTrapping::imageEffort(timelapse = timelapse_REZ)
effort_Trap_REZ <- cameraTrapping::trapEffort(cttable_REZ, camOP = list(stationCol = "Camera", setupCol = "Setup_date", retrievalCol = "Retrieval_date", hasProblems = T, cameraCol = "Camera", byCamera = F, allCamsOn = F, camerasIndependent = F))
### Display the results
effort_Image_REZ
effort_Trap_REZ
### Save the output
#openxlsx::write.xlsx(list("ImageEffort" = effort_Image_REZ, "TrapEffort" = effort_Trap_REZ), file = file.path("Results", paste("CameraEffort_REZ_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))

## Ranch-Level Cameras
### Calculate trap and image effort
effort_Image_Ranch <- cameraTrapping::imageEffort(timelapse = timelapse_Ranch)
effort_Trap_Ranch <- cameraTrapping::trapEffort(cttable_Ranch, camOP = list(stationCol = "Site", setupCol = "Setup_date", retrievalCol = "Retrieval_date", hasProblems = T, cameraCol = "Camera", byCamera = F, allCamsOn = F, camerasIndependent = F))
### Display the results
effort_Image_Ranch
effort_Trap_Ranch
### Save the Output
#openxlsx::write.xlsx(list("ImageEffort" = effort_Image_Ranch, "TrapEffort" = effort_Trap_Ranch), file = file.path("Results", paste("CameraEffort_Ranch_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))


# Step 3: Calculate number of independent events at 2 different time intervals ####
## Road Effect Zone Study
### First we need the start and end dates of the study
startdate_REZ <- "2022-05-01"
enddate_REZ <- "2023-04-30"
### Calculate number of independent events at 30 min interval
AP2_30min <- cameraTrapping::calculateEvents(AP1_REZ, envdata_REZ, start_date = startdate_REZ, end_date = enddate_REZ, interval = "30 min")
### Calculate number of independent events at 1 min interval
AP2_1min <- cameraTrapping::calculateEvents(AP1_REZ, envdata_REZ, start_date = startdate_REZ, end_date = enddate_REZ, interval = "1 min")

## Ranch-Level Comparison
### Start and End dates of the study
startdate_Ranch <- "2022-02-01"
enddate_Ranch <- "2022-05-31"
### Calculate number of independent events at 30 min interval
#### No diel activity analysis for this analysis so no distinguisher is needed for time interval
AP2_Ranch <- cameraTrapping::calculateEvents(do = AP1_Ranch, envdata = envdata_Ranch, sort.col = "Site", start_date = startdate_Ranch, end_date = enddate_Ranch, interval = "30 min")


# Step 4: Save the outputs for easy loading later ####
## Save the results of this process to use later
#save.image(file = file.path("Results", paste("Events_Data_", format(Sys.Date(), "%Y%m%d"), ".RData", sep = "")))

## Load the RData file for the calculated events data
load(file = file.path("Results", "Events_Data_20231206.RData"))

#####################################################################################################################
################################## Part 2: Prepping REZ data for PERMANOVA Analysis ################################# 
#####################################################################################################################
rm(list = ls()); gc()


# Step 0: Load the RData file for the calculated events data ####
load(file = file.path("Results", "Events_Data_20231206.RData"))


# Step 1: Check the data for errors and number of detections ####
## Number of unique species
sort(unique(AP2_30min$species))
## Number of events per species
with(AP2_30min, table(species))
## Number of events per species and location
with(AP2_30min, table(species, Location))
## Number of events by Location and distance
REZ_table <- as.array(with(AP2_30min, table(Location, Distance, species)))
REZ_table[,,c("badger", "bobcat", "coyote", "ocelot", "raccoon", "striped_skunk", "weasel")]
## Which species should be included in the analysis?
include_specs <- species_list$Timelapse[species_list$Process==1]


# Step 2: Calculate Detections/Month ####
## Running my function for calculating number of events/month (uses max number of individuals detected in the event sequence)
camop <- list(stationCol = "Camera", setupCol = "Setup_date", retrievalCol = "Retrieval_date", hasProblems = T, cameraCol = "Camera", byCamera = F, allCamsOn = F, camerasIndependent = F)
AP3_month <- cameraTrapping::summarizeEvents(x = AP2_30min, ct = cttable_REZ_all, unit = "1 month", 
                                             include = include_specs, 
                                             camOP = camop, out_form = "long", out_data = c("AB"), out_correction = "all")
AP3_month_PD <- lapply(AP3_month, function(x){x$l_AB_pd})
AP3_month_PU <- lapply(AP3_month, function(x){x$l_AB_pu})
AP3_month_RM <- lapply(AP3_month, function(x){x$l_AB_rm})
AP3_month_raw <- lapply(AP3_month, function(x){x$l_AB_raw})

# Step 3: Format data for export ####
## Function for doing so
convertFun <- function(x){
  #x <- AP3_month_AB
  
  # Add relevant predictors and covariates
  x1 <- lapply(1:length(x), function(i){
    y1 <- merge.data.frame(envdata_REZ, x[[i]], by.x = "Camera", by.y = "Site", all.y = T)
    y1$month <- paste(lubridate::year(y1$Interval), formatC(lubridate::month(y1$Interval), width = 2, flag = "0"), sep = "")
    y1$LocXDist <- with(y1, paste(Location, Distance, sep = "_"))
    y1$LocXMo <- with(y1, paste(Location, month, sep = "_"))
    y1$DistXMo <- with(y1, paste(Distance, month, sep = "_"))
    y1$LocXDistXMo <- with(y1, paste(Location, Distance, month, sep = "_"))
    y1$wholeplot <- with(y1, paste(Location, Transect, sep = "_"))
    y1$subplot <- with(y1, paste(Location, Transect, Distance, sep = "_"))
    y1$species <- names(AP3_month)[i]
    return(y1)
  })
  # Combine each species into a single object
  x2 <- do.call(rbind, x1)
  # Create a common name for the values field
  colnames(x2)[grep("AB", colnames(x2))] <- "values"
  # Convert to wide format where each species is its own column
  x3 <- tidyr::pivot_wider(x2, names_from = species, values_from = values)
  x3$Label <- with(x3, paste(Location, Transect, Distance, month, sep = "_"))
  
  # Rearrange the data and select only relevant columns
  factors_design <- c("Site", "Location", "Transect", "Distance", "month")
  factors_diag <- c("Good", "Interval", "totald", "actived", "activet")
  factors_perm <- c("LocXDist", "LocXMo", "DistXMo", "LocXDistXMo", "wholeplot", "subplot")
  factors_cov <- c("CanCover_Percent", "CanHeight_Mn", "CanHeight_Med")
  factors_spec <- unique(x2$species)
  
  x4 <- data.frame(x3[,c("Label", factors_design, factors_cov, factors_diag, factors_perm, factors_spec)])
  
  # Remove months that are not part of the study period (study period = May 2022 through April 2023)
  x5 <- x4[!(x4$month %in% c("202203", "202204", "202305")),]
  
  return(x5)
  
  # Remove those experimental units that have no data (0 active trap nights)
  #x6 <- x5[x5$actived != 0,]
  rm(x1, x2, x3, x4, x5)
}

## Run the function
AP4_month_PD <- convertFun(AP3_month_PD)
AP4_month_PU <- convertFun(AP3_month_PU)
AP4_month_rm <- convertFun(AP3_month_RM)
AP4_month_raw <- convertFun(AP3_month_raw)

## Remove those experimental units that have no data (0 active trap nights)
AP5_month_PD <- AP4_month_PD[AP4_month_PD$actived != 0,]
AP5_month_PU <- AP4_month_PU[AP4_month_PU$actived != 0,]
AP5_month_rm <- AP4_month_rm[AP4_month_rm$actived != 0,]
AP5_month_raw <- AP4_month_raw[AP4_month_raw$actived != 0,]

## Define groups of fields by their properties
factors_design <- c("Site", "Location", "Transect", "Distance", "month")
factors_diag <- c("Good", "Interval", "totald", "actived", "activet")
factors_perm <- c("LocXDist", "LocXMo", "DistXMo", "LocXDistXMo", "wholeplot", "subplot")
factors_cov <- c("CanCover_Percent", "CanHeight_Mn", "CanHeight_Med")
factors_spec <- sub("-", ".", names(AP3_month))

## Create an object with just the factors of interest
factors_all <- AP4_month_PU[,c("Label", factors_design, factors_diag, factors_perm)]
factors_rm <- AP5_month_PU[,c("Label", factors_design, factors_diag, factors_perm)]

## Create an indicators object for use in PERMANOVA
indicators <- merge.data.frame(data.frame(Timelapse = names(AP3_month)[which(names(AP3_month) %in% include_specs)]), species_list[,c("Timelapse", "Name", "Class", "Order", "Tag", "MammalAnalysis", "NativePlusAnalysis", "NativeAnalysis", "CarnivoreAnalysis", "UngulateAnalysis", "FelidAnalysis", "IndividualAnalyses")], by = "Timelapse", all.x = T)
indicators$Timelapse <- sub("-", ".", indicators$Timelapse)

## Create a covariates file for the subplot (level of variation for covariates = location x transect x distance)
### Need to be able to summarize covariates at the subplot and wholeplot levels for analysis
covariates_all <- AP4_month_PU[,c("Label", factors_cov)]
covariates_rm <- AP5_month_PU[,c("Label", factors_cov)]
#covariates <- unique(AP4_month_AB[,c("subplot", factors_cov)])

## Bring all the data together into a large output file for export and analysis
out <- list(data_PD_all = AP4_month_PD[,c("Label", factors_spec)], data_PU_all = AP4_month_PU[,c("Label", factors_spec)], data_RM_all = AP4_month_rm[,c("Label", factors_spec)], data_raw_all = AP4_month_raw[,c("Label", factors_spec)], 
            data_PD_rm = AP5_month_PD[,c("Label", factors_spec)], data_PU_rm = AP5_month_PU[,c("Label", factors_spec)], data_RM_rm = AP5_month_rm[,c("Label", factors_spec)], data_raw_rm = AP5_month_raw[,c("Label", factors_spec)], 
            factors_all = factors_all, factors_rm = factors_rm, 
            indicators = indicators, 
            covariates_all = covariates_all, covariates_rm = covariates_rm, 
            all_PD_all = AP4_month_PD, all_PU_all = AP4_month_PU, all_RM_all = AP4_month_rm, all_raw_all = AP4_month_raw, 
            all_PD_rm = AP5_month_PD, all_PU_rm = AP5_month_PU, all_RM_rm = AP5_month_rm, all_raw_rm = AP5_month_raw)

## Export the data for use in PRIMER and SAS
openxlsx::write.xlsx(out, file = file.path("Outputs_to_Other_Programs", paste("PERMANOVA_data_REZ_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))


#####################################################################################################################
####################################### Part 3: PERMANOVA and LSmeans Analyses ###################################### 
#####################################################################################################################

# PERMANOVA analyses performed in Primer v7

# LSmeans (least squares means) and associated standard errors calculated using SAS v9.4

#####################################################################################################################
##################################### Part 4: Graphs from the PERMANOVA Analyses #################################### 
#####################################################################################################################
rm(list = ls()); gc()

# Step 1: Road Effect Zone Multivariate plots from PERMANOVA results ####
## Find the data
files <- list.files(path = file.path(getwd(), "Data_from_Others"), full.names = T)
files

## Separate REZ from Ranch files
files_REZ <- files[grep("REZ_", files)]

## Separate out site scores from species scores data
files_site <- files_REZ[grep("site", files_REZ)]
files_spec <- files_REZ[grep("spec", files_REZ)]

## Prep site scores for plotting
sitescores <- pbapply::pblapply(1:length(files_site), function(i){
  name1 <- sub("site_REZ_", "", sub("_sq", "", basename(fs::path_ext_remove(files_site[i]))))
  x <- openxlsx::read.xlsx(files_site[i])
  x1 <- data.frame(scores = "site", data = name1, response = do.call(c, strsplit(name1, split = "_"))[1], group = do.call(c, strsplit(name1, split = "_"))[2], label = x[,1], x[,c("Location", "Distance", "month", "LocXDist", "LocXMo", "DistXMo", "LocXDistXMo")], 
                   x[,which(colnames(x) %in% paste("PCO", 1:4, sep = ""))])
  x1$Traffic <- factor(x1$Location, levels = c("ES", "HY"), labels = c("Low", "High"))
  x1$Dist <- factor(x1$Distance, levels = c("1", "2", "3", "4", "5", "6", "7"), 
                    labels = c("50 m", "250 m", "450 m", "650 m", "850 m", "1050 m", "1250 m"))
  x1$Month <- factor(x1$month, levels = c("202205", "202206", "202207", "202208", "202209", "202210", 
                                          "202211", "202212", "202301", "202302", "202303", "202304"), 
                     labels = c("May 22", "Jun 22", "Jul 22", "Aug 22", "Sep 22", "Oct 22", 
                                "Nov 22", "Dec 22", "Jan 23", "Feb 23", "Mar 23", "Apr 23"))
  
  x1$Response <- factor(x1$response, levels = c("PDmam", "PDnat+", "PDnat", "PDcarn", "PDung", "PDfelid"), 
                        labels = c("Mammal Community", "Native Mammal Plus Community", "Native Mammal Community", "Carnivore Community", "Ungulate Community", "Bobcat and Ocelot"))
  return(x1)
})
names(sitescores) <- sub("site_REZ_", "", sub("_sq", "", basename(fs::path_ext_remove(files_site))))
## Prep species scores for plotting
specscores <- pbapply::pblapply(1:length(files_spec), function(i){
  name1 <- sub("spec_REZ_", "", sub("_sq", "", basename(fs::path_ext_remove(files_spec[i]))))
  x <- openxlsx::read.xlsx(files_spec[i])
  x1 <- data.frame(t(x[,-1]))
  colnames(x1) <- c(x[,1])
  x2 <- data.frame(scores = "spec", data = name1, response = do.call(c, strsplit(name1, split = "_"))[1], group = do.call(c, strsplit(name1, split = "_"))[2], label = rownames(x1), apply(x1[,which(colnames(x1) %in% paste("PCO", 1:4, sep = ""))], 2, as.numeric))
  x2$Species <- factor(x2$label, 
                       levels = c("badger", "bobcat", "coyote", "domestic_cat", "domestic_dog", "grey_fox", "ocelot", "raccoon", 
                                  "striped_skunk", "weasel", "addax", "armadillo", "cottontail", "domestic_horse", "feral_hog", 
                                  "javelina", "mexican_ground_squirrel", "nilgai", "opossum", "oryx", "unk_mammal", "unk_ungulate", 
                                  "waterbuck", "white.tailed_deer"), 
                       labels = c("badger", "bobcat", "coyote", "cat", "dog", "fox", "ocelot", "raccoon", 
                                  "skunk", "weasel", "addax", "armadillo", "rabbit", "horse", "hog", 
                                  "javelina", "squirrel", "nilgai", "opossum", "oryx", "mammal", "ungulate", 
                                  "waterbuck", "deer"))
  x2$Response <- factor(x2$response, levels = c("PDmam", "PDnat+", "PDnat", "PDcarn", "PDung", "PDfelid"), 
                        labels = c("Mammal Community", "Native Mammal Plus Community", "Native Mammal Community", "Carnivore Community", "Ungulate Community", "Bobcat and Ocelot"))
  return(x2)
})
names(specscores) <- sub("spec_REZ_", "", sub("_sq", "", basename(fs::path_ext_remove(files_spec))))

## Prepare data for each graph
### Combine list objects into a single object
sitescores2 <- do.call(rbind, sitescores)
specscores2 <- do.call(rbind, specscores)
### Create plotting objects that contain the site scores and species scores
LocDist <- c("PDmam", "PDnat+")
LocMo <- c("PDmam", "PDnat+", "PDcarn", "PDung")
### Location X Distance plots
multLocDist <- list(sitescores = sitescores2[sitescores2$response %in% LocDist & sitescores2$group == "LocDist",], 
                    specscores = specscores2[specscores2$response %in% LocDist & specscores2$group == "LocDist",])
multLocMo <- list(sitescores = sitescores2[sitescores2$response %in% LocMo & sitescores2$group == "LocMo",], 
                  specscores = specscores2[specscores2$response %in% LocMo & specscores2$group == "LocMo",])

## A function for creating the plots
plotPCO <- function(ds, plot.type, X, Y, LABEL, facet.plot = FALSE, size.text, size.lab, size.point){
  #ds <- multLocDist
  #ds <- multLocMo
  #plot.type <- "LocDist"
  #plot.type <- "LocMo"
  #X <- "PCO1"
  #Y <- "PCO2"
  #LABEL <- "Species"
  #facet.plot <- FALSE
  #facet.plot <- TRUE
  #size.text <- 3
  #size.lab <- 9
  #size.point <- 4
  
  # Extract the site and species scores
  sites <- ds$sitescores
  specs <- ds$specscores
  
  spec_X <- sub("av_", "", X)
  spec_Y <- sub("av_", "", Y)
  
  # Load the ggplot2 package
  require(ggplot2)
  
  # The theme to make the graph look how we want it to
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.50, margin = margin(b = 0.5, unit = "inch")), 
          plot.background = element_blank(), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab*1.15, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = -1.0, hjust = 0.50), 
          axis.title.y = element_text(angle = 90, hjust = 0.50, vjust = 2.0), 
          axis.text = element_text(size = size.lab, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.0), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = NA), 
          panel.grid.major.x = element_line(color = NA), 
          panel.grid.major.y = element_line(color = NA), 
          panel.grid.minor = element_line(color = NA), 
          panel.grid.minor.x = element_line(color = NA), 
          panel.grid.minor.y = element_line(color = NA)) +  
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab), 
          legend.title=element_text(size = size.lab), 
          legend.position = "right", 
          #legend.key = element_rect(color = "black", fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray85", color = "black"), 
          strip.text = element_text(size = size.lab), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  # Adjust fill colors and shapes
  if(plot.type == "month"){
    COL <- "Month"
    SHAPE <- NULL
    scale_color <- scale_color_manual("Month", values = RColorBrewer::brewer.pal(12, "Paired"))
    scale_shape <- NULL
  }else if(plot.type == "distance"){
    COL <- "Dist"
    SHAPE <- NULL
    scale_color <- scale_color_manual("Distance", values = RColorBrewer::brewer.pal(7, "Blues"))
    scale_shape <- NULL
  }else if(plot.type == "LocMo"){
    COL <- "Month"
    SHAPE <- "Traffic"
    scale_color <- scale_color_manual("Month", values = RColorBrewer::brewer.pal(12, "Paired"))
    scale_shape <- scale_shape_manual("Traffic Volume", values = c(15,16))
  }else if(plot.type == "MoLoc"){
    COL <- "Traffic"
    SHAPE <- "Month"
    scale_color <- scale_color_manual("Traffic Volume", values = c("#FFFF73", "#A80000"))
    scale_shape <- scale_shape_manual("Month", values = c(15, 16, 17, 18, 0, 1, 2, 5, 6, 3, 4, 8))
  }else if(plot.type == "LocDist"){
    COL <- "Dist"
    SHAPE <- "Traffic"
    scale_color <- scale_color_manual("Distance from Hwy", values = RColorBrewer::brewer.pal(7, "Blues"))
    scale_shape <- scale_shape_manual("Traffic Volume", values = c(15,16))
  }else if(plot.type == "DistLoc"){
    COL <- "Traffic"
    SHAPE <- "Dist"
    scale_color <- scale_color_manual("Traffic Volume", values = c("#FFFF73", "#A80000"))
    scale_shape <- scale_shape_manual("Distance from Hwy", values = c(15, 16, 17, 18, 3, 4, 8))
  }else{
    stop("You must choose a valid plotting type")
  }
  
  # Set the limits of the plots based on the data
  limits <- list(xlims = c(min = plyr::round_any(min(sites[,X]), accuracy = 5, f = floor), 
                           max = plyr::round_any(max(sites[,X]), accuracy = 5, f = ceiling)), 
                 ylims = c(min = plyr::round_any(min(sites[,Y]), accuracy = 5, f = floor), 
                           max = plyr::round_any(max(sites[,Y]), accuracy = 5, f = ceiling)))
  
  # Standardize the correlations to the maximum and minimum values 
  site_split <- split(sites, f = sites$Response)
  site_split2 <- site_split[which(sapply(site_split, nrow) > 0)]
  spec_split <- split(specs, f = specs$Response)[names(site_split)]
  spec_split1 <- spec_split[which(sapply(spec_split, nrow) > 0)]
  spec_split2 <- lapply(1:length(spec_split1), function(i){
    out <- data.frame(spec_split1[[i]][,c(LABEL, "Response", spec_X, spec_Y)], 
                      ifelse(spec_split1[[i]][,spec_X] >= 0, spec_split1[[i]][,spec_X] * limits$xlims["max"], spec_split1[[i]][,spec_X] * -limits$xlims["min"]), 
                      ifelse(spec_split1[[i]][,spec_Y] >= 0, spec_split1[[i]][,spec_Y] * limits$ylims["max"], spec_split1[[i]][,spec_Y] * -limits$ylims["min"]))
    colnames(out) <- c("Species", "Response", paste(c(spec_X,spec_Y), "orig", sep = "_"), spec_X, spec_Y)
    return(out)
  })
  names(spec_split2) <- names(spec_split1)
  specs2 <- do.call(rbind, spec_split2)
  
  
  # Set data to use in the plot depending on if you want facets or not
  if(isTRUE(facet.plot)){
    ds_site <- list(sites)
    ds_spec <- list(specs2)
    facets <- facet_wrap(facets = ~ Response, ncol = 2)
  }else if(isFALSE(facet.plot)){
    ds_site <- site_split2
    ds_spec <- spec_split2
    facets <- NULL
  }
  
  # Create the plot/s depending on if you want facets or not
  plot.out <- lapply(1:length(ds_site), function(i){
    gtext <- ggrepel::geom_text_repel(data = ds_spec[[i]], mapping = aes(x = !!sym(spec_X), y = !!sym(spec_Y), label = Species), 
                                      size = size.text, max.overlaps = 20)
    
    #gtext <- geom_text(data = ds_spec[[i]], mapping = aes(x = !!sym(X), y = !!sym(Y), label = Species), 
    #                   position = position_jitter(width = 0.1, height = 0.75, seed = 4), size = size.text)
    if(is.null(SHAPE)){
      gpoint <- geom_point(data = ds_site[[i]], mapping = aes(x = !!sym(X), y = !!sym(Y), col = !!sym(COL)), size = size.point)
    }else{
      gpoint <- geom_point(data = ds_site[[i]], mapping = aes(x = !!sym(X), y = !!sym(Y), col = !!sym(COL), shape = !!sym(SHAPE)), size = size.point)
    }
    ggplot(data = ds_site[[i]]) + 
      geom_segment(mapping = aes(x = 0, y = Inf, xend = 0, yend = -Inf), col = "grey80", linewidth = 1.5) + 
      geom_segment(mapping = aes(x = Inf, y = 0, xend = -Inf, yend = 0), col = "grey80", linewidth = 1.5) + 
      gpoint + gtext + 
      scale_x_continuous(name = spec_X, breaks = seq(-50,50,5), limits = c(-50,50)) + 
      scale_y_continuous(name = spec_Y, breaks = seq(-50,50,5), limits = c(-50,50)) + 
      scale_color + scale_shape + 
      coord_cartesian(xlim = limits$xlims, ylim = limits$ylim) + 
      facets + 
      theme.plot
  })
  plot.out[[1]]
  
  if(isFALSE(facet.plot)){
    names(plot.out) <- names(site_split2)
  }
  
  if(length(plot.out) == 1){
    plot.out <- plot.out[[1]]
  }else{
    message("The output will be a list because multiple plots are created")
  }
  
  return(plot.out)
  rm(sites, specs, specs2, limits, site_split, spec_split, spec_split2, ds_site, ds_spec, SHAPE, COL, facets, gpoint, gtext, scale_color, scale_shape, theme.plot, plot.out)
  #rm(ds, plot.type, X, Y, facets, size.text, size.lab)
}

## Create the plots
sizeText <- 7
sizeLab <- 20
sizePoint <- 7
p.m.LocDist.f <- plotPCO(ds = multLocDist, plot.type = "LocDist", X = "PCO1", Y = "PCO2", LABEL = "Species", facet.plot = TRUE, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
p.m.LocDist.l <- plotPCO(ds = multLocDist, plot.type = "LocDist", X = "PCO1", Y = "PCO2", LABEL = "Species", facet.plot = FALSE, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
p.m.LocDist.c <- ggpubr::ggarrange(plotlist = p.m.LocDist.l, nrow = 1, ncol = 2, labels = LETTERS[1:length(p.m.LocDist.l)], common.legend = TRUE, legend = "top") + ggpubr::bgcolor("white") + ggpubr::border(color = NA)
p.m.LocMo.f <- plotPCO(ds = multLocMo, plot.type = "LocMo", X = "PCO1", Y = "PCO2", LABEL = "Species", facet.plot = TRUE, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
p.m.LocMo.l <- plotPCO(ds = multLocMo, plot.type = "LocMo", X = "PCO1", Y = "PCO2", LABEL = "Species", facet.plot = FALSE, size.text = sizeText, size.lab = sizeLab, size.point = sizePoint)
p.m.LocMo.c <- ggpubr::ggarrange(plotlist = p.m.LocMo.l, nrow = 2, ncol = 2, labels = LETTERS[1:length(p.m.LocMo.l)], common.legend = TRUE, legend = "right") + ggpubr::bgcolor("white") + ggpubr::border(color = NA)

## Display the plots
p.m.LocDist.f
p.m.LocDist.c
p.m.LocMo.f
p.m.LocMo.c

## Save the plots
w.pres <- 15
h.pres <- 12
### Use these for publication then shrink them with powerpoint
ggsave(filename = file.path("Results", paste("PCO_LocDist_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = p.m.LocDist.c, device = "tiff", compression = "lzw", width = w.pres, height = h.pres*0.75, dpi = 600)
ggsave(filename = file.path("Results", paste("PCO_LocMo_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = p.m.LocMo.c, device = "tiff", compression = "lzw", width = w.pres, height = h.pres, dpi = 900)


# Step 2: Road Effect Zone Univariate plots from PERMANOVA results (bar graphs version) ####
## Load the data
lsmeans <- read.csv(file.path("Data_from_Others", "LSM_REZ_20240104.csv"))
lsmeans$Location[lsmeans$Location == ""] <- NA
lsmeans$Distance[lsmeans$Distance == "_"] <- NA
lsmeans[1:5,]
lsmeans$Estimate <- ifelse(lsmeans$Estimate < 0, 0, lsmeans$Estimate)

## Create factors and factor labels for relevant data
lsmeans$traffic <- factor(lsmeans$Location, levels = c("ES", "HY"), 
                          labels = c("Low", "High"))
lsmeans$Dist <- factor(lsmeans$Distance, levels = c("1", "2", "3", "4", "5", "6", "7"), 
                       labels = c("50 m", "250 m", "450 m", "650 m", "850 m", "1050 m", "1250 m"))
lsmeans$Month <- factor(lsmeans$month, c("202205", "202206", "202207", "202208", "202209", "202210", 
                                         "202211", "202212", "202301", "202302", "202303", "202304"), 
                        labels = c("May 22", "Jun 22", "Jul 22", "Aug 22", "Sep 22", "Oct 22", 
                                   "Nov 22", "Dec 22", "Jan 23", "Feb 23", "Mar 23", "Apr 23"))
lsmeans$Species <- factor(lsmeans$species, levels = c("armadillo", "badger", "bobcat", "cotton", 
                                                      "coyote", "feral_hog", "javelina", "nilgai", 
                                                      "opossum", "raccoon", "striped_skunk", "white-tailed_deer"), 
                          labels = c("Armadillo", "Badger", "Bobcat", "Cottontail", 
                                     "Coyote", "Feral hog", "Javelina", "Nilgai", 
                                     "Opossum", "Raccoon", "Striped skunk", "White-tailed deer"))

## Create error bars and confidence intervals
lsmeans$LSE <- with(lsmeans, ifelse(Estimate - StdErr < 0, 0, Estimate - StdErr))
lsmeans$USE <- with(lsmeans, Estimate + StdErr)
lsmeans$LCL <- with(lsmeans, ifelse(Estimate - qt(0.95, DF)*StdErr < 0, 0, Estimate - qt(0.95, DF)*StdErr))
lsmeans$UCL <- with(lsmeans, Estimate + qt(0.95, DF)*StdErr)

## Determine significance levels
lsmeans$sig1 <- ""
lsmeans$sig2 <- ""
lsmeans$sig3 <- ""

## Figure out what graphs need to be made 
## Because of the extreme complexity of the interactions, we need to produce multiple graphs for some species
## For Distance effects: 
### Hog
uniHogLocDist <- lsmeans[lsmeans$Effect == "Location*Distance" & lsmeans$species == "feral_hog",]
uniHogDistMo <- lsmeans[lsmeans$Effect == "Distance*month" & lsmeans$species == "feral_hog",]
ggplot(data = uniHogLocDist, mapping = aes(x = traffic, y = Estimate, fill = Dist)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw()
ggplot(data = uniHogDistMo, mapping = aes(x = Month, y = Estimate, fill = Dist)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw()
##### CHOICE: Location X Distance Interaction
### Cottontail
uniCotton <- lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species %in% c("cotton"),]
#### Location X Distance: 202301, 202302
ggplot(data = uniCotton, mapping = aes(x = traffic, y = Estimate, fill = Dist)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
#### Distance only: 202207, 202210, 202303, 202304
ggplot(data = uniCotton, mapping = aes(x = Dist, y = Estimate, fill = Dist)) + 
  stat_summary(geom = "bar", fun = "mean") + 
  #geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
##### CHOICE: February 2023 (Location X Distance Interaction)
### Javelina
uniJavelina <- lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species %in% c("javelina"),]
#### Location X Distance: 202302, 202303
ggplot(data = uniJavelina, mapping = aes(x = traffic, y = Estimate, fill = Dist)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
#### Distance only: 202206, 202208
ggplot(data = uniJavelina, mapping = aes(x = Dist, y = Estimate, fill = Dist)) + 
  stat_summary(geom = "bar", fun = "mean") + 
  #geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
##### CHOICE: February 2023 (Location X Distance Interaction)
### Skunk
uniSkunk <- lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species %in% c("striped_skunk"),]
#### Location X Distance: 202206
ggplot(data = uniSkunk, mapping = aes(x = traffic, y = Estimate, fill = Dist)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
##### CHOICE: June 2022 (Location X Distance Interaction)
## For Traffic Volume effects: 
### Badger, Coyote, Hog, Nilgai, Opossum, Raccoon
##### CHOICE: Month X Location Interaction
### Cottontail
#### Location X Distance: 202301, 202302
ggplot(data = uniCotton, mapping = aes(x = Dist, y = Estimate, fill = traffic)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
##### CHOICE: January 2023 (Distance X Location Interaction)
### Javelina
#### Location X Distance: 202302, 202303
ggplot(data = uniJavelina, mapping = aes(x = Dist, y = Estimate, fill = traffic)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
#### Location only: 202205, 202206, 202208, 202209, 202210, 202211, 202212, 202301, 202304
ggplot(data = uniJavelina, mapping = aes(x = traffic, y = Estimate, fill = traffic)) + 
  stat_summary(geom = "bar", fun = "mean") + 
  #geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
##### CHOICE: February 2023 (Distance X Location Interaction)
### Skunk
#### Location X Distance: 202206
ggplot(data = uniSkunk, mapping = aes(x = Dist, y = Estimate, fill = traffic)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
#### Location only: 202207, 202210
ggplot(data = uniSkunk, mapping = aes(x = traffic, y = Estimate, fill = traffic)) + 
  stat_summary(geom = "bar", fun = "mean") + 
  #geom_bar(stat = "identity", position = position_dodge()) + 
  facet_wrap(facets = vars(Month), ncol = 3, nrow = 4) + 
  theme_bw()
##### CHOICE: June 2022 (Distance X Location Interaction)
rm(uniHogDistMo, uniHogLocDist, uniCotton, uniJavelina, uniSkunk)

## Create datasets for each plot
### Month Effect (Armadillo, Bobcat, Deer) - Only for supplemental material or not at all
uniMonth <- list(Armadillo = lsmeans[lsmeans$Effect == "month" & lsmeans$species == "armadillo",], 
                 Bobcat = lsmeans[lsmeans$Effect == "month" & lsmeans$species == "bobcat",], 
                 "White-tailed deer" = lsmeans[lsmeans$Effect == "month" & lsmeans$species == "white-tailed_deer",])

### Traffic Volume effect (Badger, Coyote, Hog, Nilgai, Opossum, Raccoon) + (Cottontail 202302, Javelina 202302, Skunk 202206)
uniTraffic <- list(Badger = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "badger",], 
                   Coyote = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "coyote",], 
                   "Feral hog" = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "feral_hog",], 
                   Nilgai = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "nilgai",], 
                   Opossum = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "opossum",], 
                   Raccoon = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "raccoon",], 
                   "Cottontail (Feb 23)" = lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species == "cotton" & lsmeans$month == "202302",], 
                   "Javelina (Feb 23)" = lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species == "javelina" & lsmeans$month == "202302",], 
                   "Striped skunk (Jun 22)" = lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species == "striped_skunk" & lsmeans$month == "202206",])

### Distance effect (Hog) + (Cottontail 202302, Javelina 202302, Skunk 202206)
uniDistance <- list("Cottontail (Feb 23)" = lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species == "cotton" & lsmeans$month == "202302",], 
                    "Feral hog" = lsmeans[lsmeans$Effect == "Location*Distance" & lsmeans$species == "feral_hog",], 
                    "Javelina (Feb 23)" = lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species == "javelina" & lsmeans$month == "202302",], 
                    "Striped skunk (Jun 22)" = lsmeans[lsmeans$Effect == "Locati*Distanc*month" & lsmeans$species == "striped_skunk" & lsmeans$month == "202206",])

### Conference figures (badger, bobcat, coyote, raccoon)
uniTraffic.pres <- list(Badger = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "badger",],
                        Bobcat = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "bobcat",], 
                        Coyote = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "coyote",],
                        Raccoon = lsmeans[lsmeans$Effect == "Location*month" & lsmeans$species == "raccoon",])

## A function for creating plots
plotMeans <- function(ds, plot.type, errorbar, facet, size.text, size.lab){
  #ds <- uniTraffic
  #ds <- uniMonth
  #ds <- uniDistance
  
  #plot.type <- "Traffic"
  #plot.type <- "Month"
  #plot.type <- "Distance"
  
  #errorbar <- "SE"
  #size.text <- 3
  #size.lab <- 12
  
  # Load the ggplot2 package
  require(ggplot2)
  
  # The theme to make the graph look how we want it to
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.50, margin = margin(b = 0.5, unit = "inch")), 
          plot.background = element_blank(), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab*1.15, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = -1.0, hjust = 0.50), 
          axis.title.y = element_text(angle = 90, hjust = 0.50, vjust = 2.0), 
          axis.text = element_text(size = size.lab, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.0), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = NA), 
          panel.grid.major.x = element_line(color = NA), 
          panel.grid.major.y = element_line(color = NA)) +  
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab), 
          legend.title=element_text(size = size.lab), 
          legend.position = "right", 
          #legend.key = element_rect(color = "black", fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray85", color = "black"), 
          strip.text = element_text(size = size.lab), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  # We can create error bars by calculating a Lower Confidence Limit and an Upper Confidence Limit.
  ## Note: It is important to set the minimum LCL to 0 when your data are counts/density/etc.
  if(errorbar == "SE"){
    error <- geom_errorbar(mapping = aes(ymin = LSE, ymax = USE), color = "grey50", position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.15)
  }else if(errorbar == "CI"){
    error <- geom_errorbar(mapping = aes(ymin = LCL, ymax = UCL), color = "grey50", position = position_dodge(width = 0.9), width = 0.2, linewidth = 0.15)
  }
  
  ds2 <- lapply(1:length(ds), function(i){ds[[i]]$name <- names(ds)[i]; return(ds[[i]])})
  names(ds2) <- names(ds)
  
  # Adjust fill colors and create individual plots
  if(plot.type == "Month"){
    X <- "Month"
    FILL <- "Month"
    scale_x <- xlab("Month")
    scale_fill <- scale_fill_manual("Month", values = RColorBrewer::brewer.pal(12, "Paired"))
    theme_axis <- theme(axis.text.x = element_text(angle = 90))
    legend.pos <- "none"
    rows <- 1; cols <- 3
  }else if(plot.type == "Distance"){
    X <- "traffic"
    FILL <- "Dist"
    scale_x <- xlab("Traffic Volume")
    scale_fill <- scale_fill_manual("Distance from Hwy", values = RColorBrewer::brewer.pal(7, "Blues"))
    theme_axis <- NULL
    legend.pos <- "top"
    rows <- 2; cols <- 2
  }else if(plot.type == "Traffic"){
    X <- "pass"
    FILL <- "traffic"
    scale_x <- "pass"
    scale_fill <- scale_fill_manual("Traffic Volume", values = c("#FFFF73", "#A80000"))
    theme_axis <- theme(axis.text.x = element_text(angle = 90))
    legend.pos <- "top"
    if(length(ds) == 4){rows <- 2; cols <- 2}else{rows <- 3; cols <- 3}
  }else{
    stop("You must choose a valid plotting type")
  }
  
  if(isTRUE(facet)){
    ds3 <- do.call(rbind, ds2)
    
    if(X == "pass"){
      if(nrow(ds3) == 14){
        X <- "Dist"
        scale_x <- xlab("Distance from Hwy")
      }else{
        X <- "Month"
        scale_x <- xlab("Month")
      }
    }
    
    if(errorbar == "SE"){
      ylim <- max(ds3$USE)
    }else if(errorbar == "CI"){
      ylim <- max(ds3$UCI)
    }
    limit <- plyr::round_any(ylim, accuracy = 0.05, f = ceiling)
    br <- seq(0, limit, length.out = 5)
    
    plot2 <- ggplot(data = ds3, mapping = aes(x = !!sym(X), y = Estimate, fill = !!sym(FILL))) + 
      geom_bar(stat = "identity", position = position_dodge()) + 
      scale_fill + scale_x + 
      scale_y_continuous("Average Detections/Trap Night") + 
      error + 
      coord_cartesian(ylim = c(0.0, limit)) + 
      facet_wrap(facets = vars(Species), nrow = rows, ncol = cols) + 
      theme.plot + theme_axis
    
  }else{
    plot1 <- lapply(ds2, function(x){
      if(X == "pass"){
        if(nrow(x) == 14){
          X <- "Dist"
          scale_x <- xlab("Distance from Hwy")
        }else{
          X <- "Month"
          scale_x <- xlab("Month")
        }
      }
      
      if(errorbar == "SE"){
        ylim <- max(x$USE)
      }else if(errorbar == "CI"){
        ylim <- max(x$UCI)
      }
      limit <- plyr::round_any(ylim, accuracy = 0.1, f = ceiling)
      br <- seq(0, limit, length.out = 5)
      
      ggplot(data = x, mapping = aes(x = !!sym(X), y = Estimate, fill = !!sym(FILL))) + 
        geom_bar(stat = "identity", position = position_dodge()) + 
        scale_fill + scale_x + 
        scale_y_continuous(name = "", breaks = br, limits = c(0.0, 50.0), expand = c(0,0)) + 
        error + 
        coord_cartesian(ylim = c(0.0, limit)) + 
        theme.plot + theme_axis
    })
    
    plot2 <- ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = plot1, nrow = rows, ncol = cols, legend = legend.pos, common.legend = TRUE, labels = LETTERS[1:length(plot1)], font.label = list(size = size.lab)), 
                                     left = ggpubr::text_grob("Average Detections/Trap Night", rot = 90, size = size.lab*1.15)) + ggpubr::bgcolor("white") + ggpubr::border(color = NA)
    plot2
  }
  
  return(plot2)
  rm(theme.plot, error, ds2, X, FILL, scale_x, scale_fill, theme_axis, legend.pos, rows, cols, plot1, plot2)
  #rm(ds, plot.type, errorbar, size.text, size.lab)
}

## make the graphs
sizeText.pub <- 8
sizeLab.pub <- 11
sizeText.pres <- 15
sizeLab.pres <- 20
### For publication
p.u.Traffic <- plotMeans(ds = uniTraffic, plot.type = "Traffic", errorbar = "SE", facet = FALSE, size.text = sizeText.pub, size.lab = sizeLab.pub)
p.u.Distance <- plotMeans(ds = uniDistance, plot.type = "Distance", errorbar = "SE", facet = FALSE, size.text = sizeText.pub, size.lab = sizeLab.pub)

## Display the graphs
p.u.Traffic
p.u.Distance

## Save the graphs
### For publication
w.print <- 6.5
h.print <- 4
ggsave(filename = file.path("Results", paste("Bar_Traffic_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = p.u.Traffic, device = "tiff", compression = "lzw", width = w.print, height = h.print*1.5, dpi = 900)
ggsave(filename = file.path("Results", paste("Bar_Distance_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = p.u.Distance, device = "tiff", compression = "lzw", width = w.print, height = h.print, dpi = 900)


#####################################################################################################################
########################################### Part 5: Diel Activity Analysis ########################################## 
#####################################################################################################################
rm(list = ls()); gc()

# Step 0: Load the RData file for the calculated events data ####
load(file = file.path("Results", "Events_Data_20231206.RData"))


# Step 1: Check which species are eligible for diel activity analysis ####
sort(unique(AP2_1min$species))
## Final PERMANOVA analyses revealed differences in community composition at the distance level 1, 2, 3 then no difference at 4, 5, 6, 7
AP2_1min$DistInt2 <- ifelse(AP2_1min$Distance <= 3, "close", "far")
AP2_1min$LocDistInt2 <- paste(AP2_1min$Location, AP2_1min$DistInt2, sep = "_")
with(AP2_1min, table(species, LocDistInt2))
### Species that can be run using updated distance (2 categories) and Location: 
## Armadillo
## (Only 1) Badger
## Bobcat
## Cottontail
## Coyote
## Feral Hog
## Javelina
## Nilgai
## (Only 1) Opossum
## (Almost) Oryx
## Raccoon
## (Only 1) Striped Skunk
## White-tailed Deer
actspec_locdistint2 <- c("armadillo", "badger", "bobcat", "cottontail", "coyote", "feral_hog", "javelina", "nilgai", "opossum", "oryx", "raccoon", "striped_skunk", "white-tailed_deer")


# Step 2: Run activity analysis ####
## Run activity analysis and save the outputs for quick reloading
reps <- 9999
### For Publication
#### Updated Distance categories (close/far) x Location
#actLocDistInt2 <- cameraTrapping::actFun(AP2_1min, split = "LocDistInt2", include = actspec_locdistint2, return = "species", rep = reps, pp = TRUE, cores.left = 4)
#saveRDS(actLocDistInt2, file = file.path("Results", paste("Activity_LocDistInt2_", format(Sys.Date(), "%Y%m%d"), ".RDS", sep = "")))

## Load the completed activity files
### For Publication
actLocDistInt2 <- readRDS(file.path("Results", "Activity_LocDistInt2_20231208.RDS"))


# Step 3: Statistical analyses comparing overall diel activity for each species across distance and location categories ####
## Function for comparing activity levels for each species
compAct <- function(x){
  x1 <- lapply(1:length(x$activity), function(i){
    y <- x$activity[[i]]
    y1 <- sapply(y, class)
    y2 <- y[which(y1 == "actmod")]
    n <- names(x$activity)[[i]]
    l <- names(y2)
    j <- rep(1:(length(l)-1), (length(l) - 1):1)
    k <- unlist(sapply(2:length(l), function(a){a:length(l)}))
    if(length(y2) < 2){
      y3 <- data.frame(species = n, comparison = "No comparisons made", Difference = NA, SE = NA, W = NA, p = NA)
    }else{
      y3 <- data.frame(species = n, comparison = paste(l[j], l[k], sep = "-"), activity::compareAct(y2))
    }
    rownames(y3) <- NULL
    return(y3)
  })
  x2 <- do.call(rbind, x1)
  return(x2)
  rm(x1, x2)
}

## Compare activity levels for each species
### For publication
compActs2 <- compAct(actLocDistInt2)

## Compile and export Diel Activity comparisons
#openxlsx::write.xlsx(x = compActs2, file = file.path("Results", paste("Activity_overall_comps_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))


# Step 4: Statistical analysis comparing diel activity at particular times ####
## Function for comparing activity levels at a particular time for each species
### This function is based on the activity::compareTimes() function except it is designed to compare the activity level of a given species in different groups at the same time interval. 
### It uses the same test but on multiple activity distributions instead of on multiple times within an activity distribution
compTime <- function(x, time){
  #x <- actLocDistInt2  # Output from the cameraTrapping::actFun() function
  #time <- c("00:00", "12:00")
  
  # Convert time in hms format to time in radians
  time_split <- sapply(strsplit(time, split = ":"), as.numeric)
  time_split <- apply(time_split, 2, function(x){c(x, rep(0, 3 - length(x)))})
  time_rad <- 2 * pi * (apply(time_split, 2, function(x){(x[1]*3600 + x[2]*60 + x[3])/(86400)}))
  
  x1 <- lapply(1:length(x$activity), function(i){
    y <- x$activity[[i]]
    y1 <- sapply(y, class)
    y2 <- y[which(y1 == "actmod")]
    n <- names(x$activity)[[i]]
    l <- names(y2)
    j <- rep(1:(length(l)-1), (length(l) - 1):1)
    k <- unlist(sapply(2:length(l), function(a){a:length(l)}))
    if(length(y2) < 2){
      y4 <- data.frame(species = n, time = time, comparison = "No comparisons made", Difference = NA, SE = NA, W = NA, p = NA)
    }else{
      y3 <- lapply(1:length(time_rad), function(i){
        z1 <- do.call(rbind, lapply(y2, function(a){a@pdf[a@pdf[,1] %in% time_rad[i],]}))
        dif <- z1[j,2] - z1[k,2]
        vardif <- z1[j,3]^2 + z1[k,3]^2
        W <- dif^2/vardif
        prob <- 1 - stats::pchisq(W, 1)
        res <- cbind(Difference = dif, SE = sqrt(vardif), W = W, p = prob)
        z2 <- data.frame(species = n, time = time[i], comparison = paste(l[j], l[k], sep = "-"), res)
        return(z2)
        rm(z1, dif, vardif, W, prob, res, z2)
        #rm(i)
      })
      y4 <- do.call(rbind, y3)
    }
    rownames(y4) <- NULL
    return(y4)
    rm(y, y1, y2, y3, y4, n, l, j, k)
    #rm(i)
  })
  x2 <- do.call(rbind, x1)
  return(x2)
  rm(x1, x2, time_split, time_rad)
  #rm(x, time)
}

## Compare activity levels at a given time for each species
compTimes2 <- compTime(actLocDistInt2, time = c("00:00","06:00", "12:00", "18:00"))

## Compile and export time-specific diel activity comparisons
#openxlsx::write.xlsx(x = compTimes2, file = file.path("Results", paste("Activity_time_comps_", format(Sys.Date(), "%Y%m%d"), ".xlsx", sep = "")))


# Step 5: Plot the results ####
## First, extract the necessary information for use in ggplot2
### For Publication
actLocDistInt2.plot <- cameraTrapping::actPlot(actLocDistInt2, top = "species")

## Remove oryx from plots
actLocDistInt2.plot <- actLocDistInt2.plot[actLocDistInt2.plot$species != "oryx",]


## Now, rename some factors to make them more plotting-friendly
### Updated Location x Distance Categories
actLocDistInt2.plot <- data.frame(actLocDistInt2.plot, do.call(rbind, strsplit(actLocDistInt2.plot$group, split = "_")))
colnames(actLocDistInt2.plot)[c(ncol(actLocDistInt2.plot) - 1, ncol(actLocDistInt2.plot))] <- c("Traffic", "Distance")
actLocDistInt2.plot$species <- factor(actLocDistInt2.plot$species, levels = actspec_locdistint2, labels = c("Armadillo", "Badger", "Bobcat", "Cottontail", "Coyote", "Feral hog", "Javelina", "Nilgai", "Opossum", "Oryx", "Raccoon", "Striped skunk", "White-tailed deer"))
actLocDistInt2.plot$Distance <- factor(actLocDistInt2.plot$Distance, levels = c("close", "far"), labels = c("Close", "Far"))
actLocDistInt2.plot$Traffic <- factor(actLocDistInt2.plot$Traffic, levels = c("ES", "HY"), labels = c("Low", "High"))
actLocDistInt2.plot$DistTraffic <- factor(actLocDistInt2.plot$group, levels = c("ES_close", "ES_far", "HY_close", "HY_far"), labels = c("Low, Close", "Low, Far", "High, Close", "High, Far"))

## A function to help manage all the plotting options
plotAct <- function(act, centre, type, yaxis, ci, coords, facet.plot, size.text, size.lab){
  #act <- actLocDistInt2.plot
  #centre <- "night"
  #type <- "LocDistInt"
  #yaxis <- "density"
  #ci <- FALSE
  #coords <- c(0, 0.15)
  #size.text <- 4
  #size.lab <- 12
  
  # Plotting adjustments
  require(ggplot2)
  
  ## Should the center of the graph be 12 midnight (night) or 12 noon (day)?
  if(centre == "day"){
    xunit <- "x_day"
    scale_x <- scale_x_continuous("Time", limits = c(0, 2*pi), breaks = c(0, pi/2, pi, 3*pi/2, 2*pi), labels = c("00:00", "06:00", "12:00", "18:00", "00:00"))
  }else if(centre == "night"){
    xunit <- "x_night"
    scale_x <- scale_x_continuous("Time", limits = c(-pi, pi), breaks = c(-pi, -pi/2, 0, pi/2, pi), labels = c("12:00", "18:00", "00:00", "06:00", "12:00"))
  }
  
  ## Do you want to plot density or frequency on the y axis?
  if(yaxis == "density"){
    yunit <- "y_dens"
    lclunit <- "lcl_dens"
    uclunit <- "ucl_dens"
    scale_y <- scale_y_continuous("Density", limits = c(0,0.3))
  }else if(yaxis == "frequency"){
    yunit <- "y_freq"
    lclunit <- "lcl_freq"
    uclunit <- "ucl_freq"
    scale_y <- scale_y_continuous("Frequency")
  }
  
  ## Make adjustments to your color scale
  if(type == "Distance"){
    color <- "Distance"
    scale_color <- scale_color_manual("Distance", values = RColorBrewer::brewer.pal(n = 7, name = "Blues"))
  }else if(type == "Location"){
    color <- "Traffic"
    scale_color <- scale_color_manual("Traffic Volume", values = c("High" = "#A80000", "Low" = "#FFFF73"))
  }else if(type == "LocDist"){
    color <- "DistTraffic"
    scale_color <- scale_color_manual("Dist X Volume", values = c("Low, 50 m" = "#A6CEE3", "Low, 250 m" = "#B2DF8A", "Low, 450 m" = "#FB9A99", "Low, 650 m" = "#FDBF6F", "Low, 850 m" = "#CAB2D6", "Low, 1050 m" = "#FFFF99", "Low, 1250 m" = "grey70", "High, 50 m" = "#1F78B4", "High, 250 m" = "#33A02C", "High, 450 m" = "#E31A1C", "High, 650 m" = "#FF7F00", "High, 850 m" = "#6A3D9A", "High, 1050 m" = "#B15928", "High, 1250 m" = "black"))
  }else if(type == "LocDistInt"){
    color <- "DistTraffic"
    scale_color <- scale_color_manual("Dist X Volume", values = c("#FFFF99", "#B15928", "#FB9A99", "#E31A1C"))
  }
  
  ## Modify the theme for the plot
  theme.plot <- theme(text = element_text(family = "serif")) + 
    theme(plot.title = element_text(hjust = 0.5, size = size.lab*1.25, margin = margin(b = 0.5, unit = "inch")), 
          #plot.background = element_blank(), 
          plot.margin = unit(c(.1,.1,.1,.1), "inch")) +
    theme(axis.ticks = element_line(color = NA, linewidth = 1, linetype = "solid"), 
          axis.line = element_line(color = NA, linewidth = .1, linetype = "solid"), 
          axis.title=element_text(size=size.lab, margin = margin(t = 0.25, unit="inch")),  
          axis.title.x = element_text(vjust = 0), 
          axis.title.y = element_text(angle = 90, vjust = 1.5), 
          axis.text = element_text(size = size.lab*0.75, color = "black"), 
          axis.text.x = element_text(angle = 0, hjust = 0.5), 
          axis.text.y = element_text(angle = 0, hjust = 0)) + 
    theme(panel.border = element_rect(fill = NA, color = "black"), 
          panel.background = element_rect(fill = NA, color = NA), 
          panel.grid.major = element_line(color = "grey50"), 
          panel.grid.minor = element_line(color = NA), 
          panel.spacing.x = unit(0.25, units = "inch")) + 
    theme(legend.margin=margin(c(0.15,0.15,0.15,0.15), unit = "inch"), 
          legend.background = element_rect(fill = NA, color = NA), 
          legend.text=element_text(size = size.lab*0.75), 
          legend.title=element_text(size=size.lab*0.75), 
          #legend.position = "top", 
          legend.key = element_rect(color = NA, fill = NA), 
          legend.key.height = unit(0.25,"inch"), 
          legend.key.width = unit(0.25, "inch")) + 
    theme(strip.background = element_rect(fill = "gray90", color = "black"), 
          strip.text = element_text(size = size.lab*0.75), 
          strip.text.x = element_text(margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")), 
          strip.text.y = element_text(angle = -90, margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "inch")))
  
  ## Facet plots or not?
  act_split <- split(act, f = act$species)
  act_split <- act_split[which(sapply(act_split, nrow) > 0)]
  ### Number of facets, if using facet plot
  if(length(act_split) == 1){
    cols <- 1
  }else if(length(act_split) >= 2 & length(act_split) <= 6){
    cols <- 2
  }else if(length(act_split) >= 7 & length(act_split) <= 12){
    cols <- 3
  }else if(length(act_split) >= 13){
    cols <- 4
  }else{
    stop("There are no unique species in the input data")
  }
  ### Specify the faceting
  if(isTRUE(facet.plot)){
    ds_act <- list(act)
    facets <- facet_wrap(facets = vars(species), ncol = cols)
    theme.back <- theme(plot.background = element_blank())
  }else if(isFALSE(facet.plot)){
    ds_act <- act_split
    facets <- NULL
    theme.back <- theme(panel.background = element_rect(fill = "grey90"), legend.key =element_rect(fill = "grey90"))
  }
  
  ## Should you plot the confidence intervals or not?
  line.width <- 0.75
  
  if(isTRUE(ci)){
    lcl <- geom_line(aes(y = !!sym(lclunit)), linetype = "dashed", linewidth = line.width)
    ucl <- geom_line(aes(y = !!sym(uclunit)), linetype = "dashed", linewidth = line.width)
  }else if(isFALSE(ci)){
    lcl <- NULL
    ucl <- NULL
  }
  
  plot.out <- lapply(1:length(ds_act), function(i){
    ggplot(data = ds_act[[i]], mapping = aes(x = !!sym(xunit), col = !!sym(color))) + 
      geom_line(aes(y = !!sym(yunit)), linewidth = line.width) + 
      lcl + ucl + 
      scale_x + scale_y + scale_color + 
      coord_cartesian(ylim = coords) + 
      facets + theme.plot + theme.back
  })
  plot.out[[1]]
  
  if(length(plot.out) == 1){
    plot.out <- plot.out[[1]]
  }else{
    message("The output will be a list because multiple plots are created")
  }
  
  return(plot.out)
  #rm(x1, x2, xunit, scale_x, yunit, lclunit, uclunit, scale_y, scale_color, theme.plot, plot.act)
  #rm(act, centre, yaxis, type, size.text, size.lab)
}

sizeText <- 7
sizeLab <- 20

## Create plots 
### Updated Location X Distance categories
p.act.LD2.f <- plotAct(actLocDistInt2.plot, centre = "night", type = "LocDistInt", yaxis = "density", ci = F, coords = c(0, 0.20), facet.plot = TRUE, size.text = sizeText, size.lab = sizeLab)
p.act.LD2.l <- plotAct(actLocDistInt2.plot, centre = "night", type = "LocDistInt", yaxis = "density", ci = F, coords = c(0, 0.20), facet.plot = FALSE, size.text = sizeText, size.lab = sizeLab)
p.act.LD2.c <- ggpubr::ggarrange(plotlist = p.act.LD2.l, nrow = 4, ncol = 3, labels = LETTERS[1:length(p.act.LD2.l)], common.legend = TRUE, legend = "top") + ggpubr::bgcolor("white") + ggpubr::border(color = NA)

## View plots
p.act.LD2.f
p.act.LD2.c

## Save all the plots
w <- 15
h <- 12
### For updated Location X Distance categories
ggsave(filename = file.path("Results", paste("Act_LD2_f_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = p.act.LD2.f, device = "tiff", compression = "lzw", width = w, height = h, dpi = 600)
ggsave(filename = file.path("Results", paste("Act_LD2_c_", format(Sys.Date(), "%Y%m%d"), ".tif", sep = "")), 
       plot = p.act.LD2.c, device = "tiff", compression = "lzw", width = w, height = h, dpi = 600)

