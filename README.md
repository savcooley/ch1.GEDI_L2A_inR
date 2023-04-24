# ch1.GEDI_L2A_inR

Author: Savannah Cooley
Goal: Process and analyze GEDI L2A data
Description: This script runs a GEDI data processing function that filters the data then spatially subsets the data over a user-defined area of interest and saves the data as .csv files.

We applied the following noise filtering steps based on the methodology outlined by Potapov et al. (2020): (a) power beam mode only, (b) beam sensitivity ≥0.9 and (c) ≤2m range of predicted ground elevations among the five algorithms after removing the algorithm that yielded the largest outlier elevation. For (c), we deviated from Potapov et al. (2020) and chose this approach to retain some observations that otherwise would have been removed if applying the ≤2m range threshold to all six algorithms. 

