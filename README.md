Source models of 29 embedded low-mass young protostellar objects, all observed as part of the "Water in Star-forming regions Herschel" Guaranteed Time Key Program on Herschel (van Dishoeck et al. 2011, Kristensen et al. 2012). The source models are all 1D, spherically symmetric and calculated using the DUSTY code (Ivezic & Elitzur 1997). Each model consists of a density and temperature profile. Please see Kristensen et al. (2012) for a more thorough discussion (freely available here: http://arxiv.org/abs/1204.0009). If the model results are used, please cite Kristensen et al. (2012). 

The document dusty_cat.pdf contains a description of how the models were run as well as the underlying assumptions, a summary of the relevant SCUBA data, and the final results. This document should always be the starting point for anyone wanting to use the source models. 

This repository contains three directories with their respective subdirectories:

- sources: The model input and output files are stored here for each source. The most relevant files are res_source.dat and mod_input.rtb. The former contains the physical properties of each source and the latter the temperature profile. The emperature profile is stored in the column labeled Td as a function of the dimensionless parameter y, where y = r / r_in and r_in is available in the res_source.dat file. The density profile can be directly obtained either from Kristensen et al. (2012) or dusty_cat.pdf and is of the form n = n0 * (r / r0)^-p; n0, r0 and p are all available from one of the two sources. 

- input: All the input data are stored in this directory, including the SCUBA 450 and 850 micron images. A summary of the input data is also provided in each source directory. 

- scripts: The relevant IDL scripts for reproducing the above analysis. 
