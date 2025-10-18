# PhloemSinksModel
This repository contains the MATLAB implementation of the Phloem model with different sink profiles along the phloem pathway, used to simulate sugar transport and pressure dynamics along the phloem with distributed sink profiles. The model couples osmotic flow, sucrose transport, and xylem–phloem water exchange to explore how spatial variations in sink demand influence phloem function.


main.m: main driver script that runs the complete phloem simulation.
plotResults.m: generates the figures and post-processing visualizations.
BCI_hydro_Cohort_hourly.csv: input data from BiomeE providing hourly hydraulic and carbon flux forcing.
/fcn/: supporting routines for model equations
/initial/: scripts to approximate an initial profile for the variables
/solver/: routines to perform Newton's method
