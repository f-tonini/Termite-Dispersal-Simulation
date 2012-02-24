***************************************************************************************************
* Please read carefully the following information and instructions before running any simulation. *
***************************************************************************************************

RECOMMENDATION: It is strongly recommended to run this code only for small-scale/localized problems AND 
		not to extend the simulation for too many years.


*************************************************************************************************************************
* This file explains what steps need to be taken in order to download and install the required software. *
* Also, this shows you how to create and properly set up the workspace environment for the IBM simulation.              *
*************************************************************************************************************************


DOWNLOAD AND INSTALLATION OF THE R SOFTWARE:
___________________________________________

(1) Go to the offical R software webpage: http://www.r-project.org/
(2) On the left menu, under 'Download, packages' click on the CRAN link
(3) Select a mirror near your country/state
(4) Click on 'Download R for Mac OS X'
(5) If you selected the MacOS X environment, click on the .pkg R version you want to download. 
    Download the most recent version of R or click on 'old' for older ones.  

********************************************************************************************************
NOTE: As specified on the main R script this code has been successfully tested using R 2.14.0 64-bit. 
      It is recommended that you do not use a version older than 2.13
********************************************************************************************************

(6) Double-click on the downloaded R version to start the setup and install the software



SETTING UP THE WORKING DIRECTORY AND RUNNING THE CODE:
______________________________________________________


(1) Create a folder on your Desktop and name it "Temp"

(2) Copy into this folder both the 'IBM_MonteCarlo.r' and 'myfunctionsMC.r' files you have downloaded

(3) Create a folder within 'Temp' called "Raster" 

(4) Inside the 'Temp' directory, place either a vector or a raster layer representing non-suitable habitat 
   for the termite species

**********************************************************************************************************************
NOTE: at this stage of code development, the only format accepted for vector layers is shapefile (.shp). 
For raster layers you can use any of the following: ('.grd','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img')
**********************************************************************************************************************
	 
(5) Inside the 'Temp' directory, place either a .txt (tab-delimited) or .csv file containing the coordinates 
    of all the colonies you want to start your simulation from (geographic or projected coordinates).The file name MUST 
    NOT HAVE BLANK SPACES and NOT CONTAIN the symbol '.' Otherwise, the code will not recognize the file correctly.


******************************************************************************************************************
NOTE: make sure to include a header in the coordinate file to label the coordinates. 
Associate the labels with the correct column...the code will assume that under the column labeled, for example, 
'latitude' you pasted the right numbers corresponding to 'latitude' values, not to a 'longitude' ones!!
*****************************************************************************************************************

(6) On console type setwd ('~/Desktop/Temp') and click ENTER. Then type source('IBM_MonteCarlo.r') and press ENTER 


(7) The simulation starts: enter all parameters you are asked for and read carefully all the options you see

*************************************************************************************************************
NOTE: Should a window pop up asking to select a mirror to download packages from, pick one that is close to
      your state/country and click OK
*************************************************************************************************************


(8) While the simulation is running you should see a progress bar advancing up to the last year of the simulation. 
    At the same time, on console you can always keep track of what is happening and read possible warnings and error messages
   
(9) After the simulation is over, you should read 'Simulation Is Over' on the console window 
    and find all the results in the 'Raster' folder



*******************************************************************************************************************
NOTE: Ignore all warnings concerning 'duplicate points', since with thousands/millions of individuals it is likely 
      to have coincident coordinates, especially where colony density is higher.
*******************************************************************************************************************




**************************************************************
Contact me for any further information or error messages:
f_tonini@hotmail.com
**************************************************************






