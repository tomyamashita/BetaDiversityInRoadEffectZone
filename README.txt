This README.txt file was generated on 2023-01-10 by Thomas J. Yamashita


GENERAL INFORMATION

1. Title of Dataset: Influence of traffic volume on mammal beta diversity within the road effect zone
	This data is associated with the manuscript, titled Influence of traffic volume on mammal beta diversity within the road effect zone, available here: https://github.com/tomyamashita/BetaDiversityInRoadEffectZone

2. Author Information
	Thomas J. Yamashita
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville
		tjyamashta@gmail.com
		Corresponding Author
	David B. Wester
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville
	Zachary M. Wardle
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville
	Daniel G. Scognamillo
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville
	Landon R. Schofield
		East Foundation
	Michael E. Tewes
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville
	John H. Young Jr. 
		Environmental Affairs Division, Texas Department of Transportation
	Jason V. Lombardi
		Caesar Kleberg Wildlife Research Institute, Texas A&M University - Kingsville

3. Date of Data Collection: 
	Camera Data: May 2022 - April 2023

4. Geographic location of data collection: Private ranches in Willacy and Kenedy Counties, Texas, USA

5. Funding Sources: Texas Department of Transportation


DATA & FILE OVERVIEW

1. File List: 
	CTtable_REZ_20230601.xslx: CT Table describing active and total camera trap nights for the duration of the study period. Specific camera locations are redacted to protect private landowner privacy and locations of endangered and threatened species. 
	EnvData_REZ_20230509.xlsx: File containing information about experimental units and vegetation data for each site. 
	Species_List_20230509.xlsx: List of all species, including scientific names detected on cameras and those that were included in the analyses. 
	README.txt: This file

2. Relationship between files: 
	The camera names are the same in the above files and can be used to link the data

3. Files not included in this repository
	Camera data from this study. This data is available upon request. Please contact Thomas Yamashita


METHODOLOGICAL INFORMATION

1. Description of the methods used for collection/generation and processing of data: 
	Methodology for collection and processing of the data can be found in the manuscript

2. Quality Assurance Procedures: 
	Photographs were processed by trained personnel and went through multiple rounds of checks to ensure accurate detection of all animals 

3. People involved with data collection, processing, and analysis: 
	Thomas J. Yamashita (data collection, data processing, data analysis)
	Zachary M. Wardle (data collection)
	Daniel G. Scognamillo (data collection)
	Jason V. Lombardi (data collection, data analysis)


DATA SPECIFIC INFORMATION FOR: CTtable_REZ_20230601.xlsx

1. Data Type: Microsoft excel xlsx file

2. Number of Variables: 19

3. Number of Rows: 98

4. Variable List: 
	Camera: Camera Name
	Transect: The transect number for the camera
	Location: Ranch the camera was on
	utm_y: Y coordinate in UTM zone 14N for the camera. This information is redacted and is available upon request.
	utm_x: X coordinate in UTM zone 14N for the camera. This information is redacted and is available upon request.
	Setup_date: Date camera was set up
	Retrieval_date: Date camera was retrieved
	ProblemX_from: (multiple columns): Start date that a camera was inactive. X starts at 1 and increases by 1. Higher numbers indicate that a camera had more problems
	ProblemX_to: (multiple columns): End date that a camera was inactive. Just like above, X starts at 1 and increases by 1 as a camera has more problems


DATA SPECIFIC INFORMATION FOR: EnvData_REZ_20230509.xlsx
1. Data Type: Microsoft excel .xlsx file

2. Number of Variables: 6

3. Number of Rows: 98

4. Variable List: 
	Camera: Camera name
	Site: Secondary unique identifier for the camera
	Location: property the camera was on, hence traffic volume associated with the site
	Transect: transect number on each property
	Distance: Distance category of each camera (1 - 7, corresponding with distances in m from 50 - 1250. One camera was placed every 200 m). 
	CanHeight_Med: Median canopy height, calculated using publicly available LiDAR data as described in the manuscript.  


DATA SPECIFIC INFORMATION FOR: Species_List_20230509.xlsx
1. Data Type: Microsoft excel .xlsx file

2. Number of Variables: 13

3. Number of Rows: 62

4. Variable List: 
	Species: common name identifier (underscores are used instead of spaces)
	Timelapse: common name used in the Timelapse2 program for classifying species
	Name: full common name of the species
	Class: Taxonomic class
	Order: Taxonomic order
	Family: Taxonomic family
	ScientificName: Full scientific name for the species
	Tag: Identifier for native, exotic, or invasive species. Also included in this are names that represent multiple species (group)
	MammalAnalysis: Indicator if the species was included in the full mammal community analysis
	NativePlusAnalysis: Indicator if the species was included in the native plus community analysis
	CarnivoreAnalysis: Indicator if the species was included in the carnivore community analysis
	UngulateAnalysis: Indicator if the species was included in the ungulate community analysis
	IndividualAnalyses: Indicator if the species was one of the individual species analyses


DATA SPECIFIC INFORMATION FOR: RoadEffectZone_Manuscript_Rcode.R
1. Data Type: R script

5. Other Information: This script provides code for conducting all R-based analyses used in this manuscript. This code will not run properly without the camera data which is not provided in this repository. Camera data is available upon request from Thomas Yamashita. 


DATA SPECIFIC INFORMATION FOR: RoadEffectZone_Manuscript_SAScode.R
1. Data Type: SAS script

5. Other Information: This script provides code for conducting all SAS-based analyses used in this manuscript. This code will not run properly without the camera data which is not provided in this repository. Camera data is available upon request from Thomas Yamashita. 
