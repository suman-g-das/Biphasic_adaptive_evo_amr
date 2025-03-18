Codes accompanying paper “Epistasis-mediated compensatory evolution in a fitness landscape with adaptational tradeoffs” by Suman G Das, Muhittin Mungan, Joachim Krug.
The code adaptive-walk.f90 simulates averages over random walks on TIL landscape.

•	The code is written in Fortran and can be compiled with the gfortran compiler from terminal. 
•	The model and simulation parameters are in Line 238-243. Parameters not defined here were held constant throughout the paper.
•	The output is printed into the file “output.d”. 
•	The output has 8 columns. In order: time step, mutation number, log fitness, log null-fitness, log fitness, log resistance, fixation probability of fitter neighbors, relative cost.
•	Un-comment line 529 to print landscape averaged -u and v at end of file. 
•	The code is set for Kimura walk. To switch to unform walk, comment out Line 421 and un-comment Line 422.

The code wright-fisher-til.py simulates Wright-Fisher dynamics on the TIL fitness landscape.


•	The code is written in python Version 3.
•	The parameter “conc” is the mean value of the stress variable.
•	The parameter “fluc” is the standard deviation of stress variable.
•	The output is printed into the file “file.d”. 
•	The output has 5 columns: In order: generation number, population mean mutation number, population mean null-fitness, population mean resistance, population mean fitness. The results are averaged over landscapes. 





