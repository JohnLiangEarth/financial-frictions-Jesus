% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

clear all
close all
clc

diary off;
delete('_log.txt');
diary('_log.txt');

run ('../3_model_NN/b1_parameters')
load '../1_steadystate/ss_results'

tic
b2_Klm
disp('model solution: completed')
toc

b9_plot

diary off;
close all

% save 'z_FinalWorkspace_gsim.mat' g_big -mat -v7.3
clear g_big
save 'z_FinalWorkspace.mat' -mat -v7.3

