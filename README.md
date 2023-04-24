# ch1.GEDI_L2A_inR

Author: Savannah Cooley

Goal: Process and analyze GEDI L2A data

Description: This R script runs a GEDI data processing function that filters the data then spatially subsets the data over a user-defined area of interest and saves the data as .csv files.

We applied the following noise filtering steps based on the methodology outlined by Potapov et al. (2020): (a) power beam mode only, (b) beam sensitivity ≥0.9 and (c) ≤2m range of predicted ground elevations among the five algorithms after removing the algorithm that yielded the largest outlier elevation. For (c), we deviated from Potapov et al. (2020) and chose this approach to retain some observations that otherwise would have been removed if applying the ≤2m range threshold to all six algorithms. 

 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
 
Autora: Savannah Cooley

Objetivo: Procesar y analizar datos de GEDI L2A

Descripción: este script en R ejecuta una función de procesamiento de datos GEDI que filtra los datos y luego los subdivide espacialmente en un área de interés definida por el usuario y guarda los datos como archivos .csv.

Aplicamos los siguientes pasos de filtrado de ruido basados en la metodología descrita por Potapov et al. (2020): (a) modo de haz de potencia únicamente, (b) sensibilidad del haz ≥0,9 y (c) rango de elevaciones del suelo previstas de ≤2 m entre los cinco algoritmos después de eliminar el algoritmo que produjo la mayor elevación atípica. Para (c), nos desviamos de Potapov et al. (2020) y eligió este enfoque para conservar algunas observaciones que, de otro modo, se habrían eliminado si se aplicara el umbral de rango de ≤2 m a los seis algoritmos.
 
