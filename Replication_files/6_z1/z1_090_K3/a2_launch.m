% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution
% This script launches the routines in charge of solving the model, plotting the results and saving them

clear all
close all
clc

diary off;
delete('_log.txt');
diary('_log.txt');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% special stuff just for this folder
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

run('../../3_model_NN/b1_parameters')
load '../../1_steadystate/ss_results'

load_PLM       = 'load ''../z1_090/z_FinalWorkspace.mat'' PLM PLM_2 PLM_finegrid PLM_finegrid_2 PLM_visits;';
load_LR_errors = 'load ''../../../2_model_LR/z_FinalWorkspace.mat'' Y Y_fit;';

wePLM=0.1;

z1=0.9;
z2            = 1 + la2/la1 * (1-z1);

maxitHJB = 1000;                             % max number of iterations for the HJB
weHJB=0.03;

Bmin=0.04;
Bmax=1;
Nmin=3;
Nmax=4;

nval_a        = 501;                        % number of points in amin-to-amax range (individual savings)
nval_z        = 2;                          % number of options for z (the idiosincratic shock)
nval_B        = 3;                          % number of points in Bmin-to-Bmax range (aggregate savings), on the coarse grid used for the HJB
nval_N        = 51;                         % number of points in Nmin-to-Nmax range (aggregate equity) , on the coarse grid used for the HJB

b1_parameters_redo


% copyfile '../../3_model_NN/b2_Klm.m';     % this one is different because here we switch from nval_B=4 to nval_B=3
  copyfile '../../3_model_NN/b3_HJB.m';
  copyfile '../../3_model_NN/b5_KFE.m';
  copyfile '../../3_model_NN/b5_use_parallel.m';
  copyfile '../../3_model_NN/b7_PLM.m';
  copyfile '../../3_model_NN/b7_PLM_iter.m';
  copyfile '../../3_model_NN/b9_plot.m';
  copyfile '../../3_model_NN/f1_NN_loss.m';
  copyfile '../../3_model_NN/f2_NN_eval.m';
  copyfile '../../3_model_NN/f5_NN_gradient.m';
  copyfile '../../3_model_NN/f8_KFE_sim.m';  % this file controls if the parallel toolbox is used or not
  
  cd graphs_v3
  copyfile '../../../3_model_NN/graphs_v3/a2_launch.m';
  copyfile '../../../3_model_NN/graphs_v3/c2_PLM_NPLM_phase.m';
  copyfile '../../../3_model_NN/graphs_v3/c3_sss.m';
  copyfile '../../../3_model_NN/graphs_v3/c4_errors_zone.m';
  copyfile '../../../3_model_NN/graphs_v3/c5_IRF.m';
  copyfile '../../../3_model_NN/graphs_v3/c5_IRF_sim.m';
  cd ..



%% %%%%%%% %%
% solve model
%% %%%%%%% %%

tic
b2_Klm
disp('model solution: completed')
toc

% plot and save results

warning('off','all')

b9_plot

save 'z_FinalWorkspace.mat' Y Y_fit Bsim Nsim -mat -v7.3

cd graphs_v3
a2_launch
close all
cd ..

warning('on','all')

close all

% save 'z_FinalWorkspace_gsim.mat' g_big -mat -v7.3
clear g_big
save 'z_FinalWorkspace.mat' -mat -v7.3

diary off;


