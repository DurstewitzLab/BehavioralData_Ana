# Welcome

This folder contains the MatLab code and data files from

[Koppe, G, Mallien, AS, Berger, S, Bartsch, P., Gass, B, Durstewitz, D (2017) CACNA1C gene regulates behavioral strategies in operant rule learning. PLOS Biology.](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2000936)

Copyright: © 2017 Georgia Koppe.


1. Raw behavioral data
Behavioral data recorded during this experiment can be found in the subfolder 'data' and contains one file for each experimental animal. The title of each file contains information on which experimental group the animal belongs to with group1=Cav1.2NesCre, group2=Cav1.2fl/fl, group3=NesCre, and group4=wild-type. 

2. Figure Data
Numerical values used to generate the figures can be found in the subfolder 'figure data' which contains a 'figure_data.mat' (MATLAB) file. Within this file the data used to generate the individual figures can be found in a variable which is labeled according to the figure.

3. Code for analysis
MATLAB code used to generate figures or conduct analyses can be found in the subfolder 'code'. This folder contains several MATLAB files for generating the corresponding figures in the manuscript (e.g. make_figs_2A_2B_3A_3B.m, make_figs_3C.m, etc.), perform reinforcement learning model analyses (see subfolder 'RLmodels'), or construct bootstrap distributions (see subfolder 'Bootstraps'). Analyses are started with the files starting with 'main_'.

Please contact georgia.koppe@zi-mannheim.de for further information.
