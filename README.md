# Network of Coupled Crystal Oscillators
C/C++ codes for the Coupled Crystal Oscillator Project

CCOST_colored_noise.cpp simulates the current of a unidirectional
array of cyrstal oscillators of a user specified size. 
This code produces data for plotting in CCOSTts.dat. 
One can plot this using the python code ccost_ts_plot.py

CCOST_phasedrift_parallel.cpp calculates the simulated phase error for increasing
sized arrays of crystal oscillators with fixed coupling strength. 
This program creates the data file ccost_phase_error.dat and is plotted by 
phase_error_plt.py. 
This program can be compiled using the makefile. 
Currently this program is set up to take 100 samples (NIT), using 25 threads (workers). Go into 
CCOST_phasedrift_parallel.cpp to edit this. 

CCOST_phasedrift_totalvariation.cpp calculates the phase error as a function of 
coupling strength and size of oscillator array. *Note* this code is underconstruction and
is not ready to be implemented. 


