***************************************************************************************************
* Please read carefully the following information and instructions before running any simulation. *
***************************************************************************************************

RECOMMENDATION: It is strongly recommended to run this code only for small-scale/localized problems AND 
		not to extend the simulation for too many years.


****************************************************************************************************************
* This file explains what steps need to be taken in order to download and install the required software.       *
* Also, this shows you how to create and properly set up the workspace environment for the IBM simulation.     *
****************************************************************************************************************


DOWNLOAD AND INSTALLATION OF THE R SOFTWARE:
___________________________________________

(1) Go to the offical R software webpage: http://www.r-project.org/
(2) On the left menu, under 'Download, packages' click on the CRAN link
(3) Select a mirror near your country/state
(4) Click on 'Download R for Windows'
(5) Click on 'base' (base binary distribution). 
   Download the most recent version of R or click on 'previous releases' for older ones. 

********************************************************************************************************
NOTE: As specified on the main R script this code has been successfully tested using R 2.14.0 64-bit. 
      It is recommended that you do not use a version older than 2.13
********************************************************************************************************

(6) Double-click on the downloaded R version to start the setup and install the software



SETTING UP THE WORKING DIRECTORY AND RUNNING THE CODE:
______________________________________________________


(1) Create a folder under C:\ and call it "Temp" (full path: C:\Temp) 

(2) Copy into this folder both the "IBM_MonteCarlo.r" and "myfunctionsMC.r" files you have downloaded

(3) Inside the "C:\Temp" directory, place either a vector or a raster layer representing non-suitable habitat 
    for the termite species

**********************************************************************************************************************
NOTE: at this stage of code development, the only format accepted for vector layers is shapefile (.shp). 
For raster layers you can use any of the following: ('.grd','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img')
**********************************************************************************************************************
	 
(4) Inside the "C:\Temp" directory, place either a .txt (tab-delimited) or .csv file containing the coordinates 
    of all the colonies you want to start your simulation from (geographic or projected coordinates).The file name MUST 
    NOT HAVE BLANK SPACES NOR CONTAIN the symbol '.' Otherwise, the code will not recognize the file correctly.
    Also, if you use projected coordinate systems use label 'X' for easting and label 'Y' for northing.


******************************************************************************************************************
NOTE: make sure to include a header in the coordinate file to label the coordinates. 
Associate the labels with the correct column...the code will assume that under the column labeled, for example, 
'latitude' you pasted the right numbers corresponding to 'latitude' values, not to a 'longitude' ones!!
*****************************************************************************************************************
    
(5) On console type setwd('C:/Temp') and click ENTER. Then type: source('IBM_MonteCarlo.r') and press ENTER  
 

(6) The simulation starts: enter all parameters you are asked for and read carefully all the options you see

*************************************************************************************************************
NOTE: Should a window pop up asking to select a mirror to download packages from, pick one that is close to
      your state/country and click OK
*************************************************************************************************************

(7) While the simulation is running you should see a progress bar advancing up to the last year of the simulation. 
    At the same time, on console you can always keep track of what is happening and read possible warnings and error messages
   

(8) After the simulation is over, a message will pop up stating that 'The Simulation Is Over'. 
    Final statistics will be computed and saved to a .csv file into the main folder.

(9) Finally, go back to the console and answer the last question. In this step, you will be able to compute, if you wish,
    some occupancy threshold raster. That is, a separate raster file for each year of the simulation, whose cells will
    have a value of 1 if they have been occupied AT LEAST X% of all Monte Carlo simulation runs you chose. X% represents the threshold
    that you have to pick in order to compute the occupancy raster. For example, let's assume you picked 100 Monte Carlo runs.
    By choosing a threshold X = 50% you will get a raster for each year of the simulation, whose occupied cells (value = 1) 
    will represent those that have been occupied AT LEAST 50 out of 100 Monte Carlo runs or more.
   

    You will find all output raster files inside the 'Raster' folder. 


*******************************************************************************************************************
NOTE: Ignore all warnings concerning 'duplicate points', since with thousands/millions of individuals it is likely 
      to have coincident coordinates, especially where colony density is higher.
*******************************************************************************************************************




**************************************************************
Contact me for any further information or error messages:
f_tonini@hotmail.com
**************************************************************









