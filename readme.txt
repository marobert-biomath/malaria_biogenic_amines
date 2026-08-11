% Readme file for the Matlab code to generate Figures 3-8. 
% Manuscript: "Mathematical modeling predicts increased malaria incidence and prevalence from 
% biogenic amine-induced changes in mosquito biology and infection."

The file baseline_biogenic_amine_model_explorer_multi_type_v20260730.m
was used to generate the results shown in Figures 3-8. 

%%% Figures 3-7 %%%

Figure 3
choice_mat(2,:) = [0 0 1 0 0 0]; 

Figure 4
choice_mat(2,:) = [0 0 0 1 1 0]; 

Figure 5
choice_mat(2,:) = [0 0 0 0 0 1]; 

Figure 6
choice_mat(2,:) = [1 0 0 0 0 0]; 

Figure 7 
choice_mat(2,:) = [1 0 1 1 1 1]; 

%%% Figure 8 
This figure includes P. Falciparum values for transmission from 
host to vector. 

Figure 8 
The P. falciparum values are generated from the baseline file 
using the following set up for choice_mat.
 
choice_mat(1,:) = [2 0 0 0 0 0]; 
choice_mat(2,:) = [3 0 1 1 1 1];

The P. yoelli values are the same as those in Figure 7.

Figure 8 was assembled using model output and data.  The file figure08_combination_figure.m describes this figure. 
This code requires two sets of files from the Matlab Stack Exchange.
These files were both downloaded 26 September 2025. 

Violinplot-Matlab v 1.0.0: https://www.mathworks.com/matlabcentral/fileexchange/170126-violinplot-matlab
plot_stats v 1.0.0: https://www.mathworks.com/matlabcentral/fileexchange/92613-statistical-significance


Supplementary Material

Figure S1 was generated with the file JULY2026_bioamines_model_for_fitting.m.

Figure S2 and S3 were generated using the file figure_si_sensitivities_final.m.  Data files 
that this code reads in were generated through model simulations with different parameter values and output
from these simulations were synthesized in this file. 