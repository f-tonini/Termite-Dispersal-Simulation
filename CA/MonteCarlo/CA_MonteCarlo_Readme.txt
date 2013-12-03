***************************************************************************************************
* Please read carefully the following information and instructions before running any simulation.    *
***************************************************************************************************

*************************************************************************************************************************
* This file explains what steps need to be taken in order to download and install the required software.                                *
* Also, this shows you how to create and properly set up the workspace environment for the CA simulation.                        *
*************************************************************************************************************************


DOWNLOAD AND INSTALLATION OF THE R SOFTWARE:
___________________________________________

(1) Go to the offical R software webpage: http://www.r-project.org/
(2) On the left menu, under 'Download, packages' click on the CRAN link
(3) Select a mirror near your country/state
(4) Click on 'Download R for Windows'
(5) Click on 'base' (base binary distribution). 
   Download the most recent version of R or click on 'previous releases' for older ones. 

********************************************************************************************************
NOTE: As specified on the main R script this code has been successfully tested using R 2.15.0 64-bit. 
      It is recommended that you do not use a version < 2.13.0
********************************************************************************************************

(6) Double-click on the downloaded R version to start the setup and install the software



SETTING UP THE WORKING DIRECTORY AND RUNNING THE CODE:
______________________________________________________


(1) Create a folder on your computer and paste the folder path into the main R script

(2) Copy into this folder both the "CA_MonteCarlo.r" and "myfunctionsMC_CA.r" files you have downloaded

(3) If you desire to run the model considering suitable/unsuitable habitat, place in the same folder either a vector or a raster layer representing non-suitable habitat 
    for the termite species. Call that layer "habitat_block" (if you choose a different name, remember to change it also in the main R script)

(3b) If you DO NOT want to use/have an unsuitable habitat layer, then remember to out-comment in the main script the lines that read
        the habitat layer

**********************************************************************************************************************
NOTE: at this stage of code development, the only format accepted for vector layers is shapefile (.shp). 
For raster layers you can use any of the following: ('.grd','.asc','.sdat','.rst','.nc','.tif','.envi','.bil','.img')
**********************************************************************************************************************
	 
(4) Inside the folder place either a .txt (tab-delimited) or .csv file containing the coordinates 
    of all the colonies you want to start your simulation from (geographic or projected coordinates).The file name MUST 
    NOT HAVE BLANK SPACES NOR CONTAIN the symbol '.' Otherwise, the code will not recognize the file correctly.
    Also, if you use projected coordinate systems use label 'X' for easting and label 'Y' for northing.


******************************************************************************************************************
NOTE: make sure to include a header in the coordinate file to label the coordinates. 
Associate the labels with the correct column...the code will assume that under the column labeled, for example, 
'latitude' you pasted the right numbers corresponding to 'latitude' values, not to a 'longitude' ones!!
*****************************************************************************************************************
    
(5) Open the R sofware and in the top menu select File >> Source R code... and point to the R script you copied in the folder you created
 

(6) The simulation starts: enter all parameters you are asked for and read carefully all the options you see

*************************************************************************************************************
NOTE: Should a window pop up asking to select a mirror to download packages from, pick one that is close to
      your state/country and click OK
*************************************************************************************************************

(7) While the simulation is running you should see a progress bar advancing up to the last year of the simulation. 
    At the same time, on console you can always keep track of what is happening and read possible warnings and error messages
   
(8) After the simulation is over, final statistics will be computed and saved to a .csv file into the main folder.
    You will also find all output simulation files inside the 'Output' folders 


**************************************************************
Contact me for any further information or error messages:
f_tonini@hotmail.com
**************************************************************



