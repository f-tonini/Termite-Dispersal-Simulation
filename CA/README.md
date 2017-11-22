[Simple](Simple): this subfolder cointains all the .r script files and readme.txt file that should be used in case you are 
            interested in running a SINGLE simulation over a given time span (e.g. from 1990-2000, etc.).

[MonteCarlo](MonteCarlo): this subfolder contains all the .r script files and readme.txt file that should be used in case you are
               interested in running MULTIPLE (MonteCarlo) simulations over the same time span (e.g. from 1990-2000, etc.).

The MonteCarlo simulation is recommended if you are interested in having a robust and more precise estimation
of the areas that could potentially be infested over time. This because a single simulation run only gives one
of many possible outcomes. Instead, a MonteCarlo simulation run for multiple times (the more the better) averages
several different outcomes for each year and return a raster with all areas occupied multiple times.

I would recommend to use [Simple](Simple) only if you want to have a first hands-on experience on the model 
input/output and run a single simulation to test the code.

On the other hand, if you are already familiar with the aforementioned concepts and want to have a more robust 
result, choose the [MonteCarlo](MonteCarlo).


**NOTE**: whether you use [Simple](Simple) or the [MonteCarlo](MonteCarlo) case, try to limit the size of the time span. The longer
         the time span, the longer the simulation time will be. The same recommendation holds for chosen number of
	 MonteCarlo runs, even though a minimum of 50 or 100 runs is suggested to have robust results.


**NOTE**: The total simulation time strongly depends on your computer its technical specs. Keep in mind that even 
	 with powerful computers (high processor speed, and RAM) the total time may vary based on the amount of
	 points of invasion you decide to start the simulation with.
