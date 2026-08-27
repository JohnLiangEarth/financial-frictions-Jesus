% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all
clear

copyfile '../3_model_NN/f2_NN_eval.m';

%% run likelihood for every solution of the model, and save results on .mat

for it_estim=130:2:150
    if it_estim ~= 140
        eval(['load ../4_sigma/sigma_' num2str(it_estim) '/z_FinalWorkspace.mat'])
    else
        load ../3_model_NN/z_FinalWorkspace.mat
    end
    load '../z_obs_data_bandpass.mat'
    Y_ss = Zeta*(B_ss+N_ss).^alpha;    % we renormalize data around the deterministic steady state (becuase it's the one that's constant across all model solutions)
    Y_obs = Y_obs*Y_ss;                % in the data file, Y_obs is defined as Y_t/Ytrend_t, to make this renormalization easy
    d5_likelihood; % calculate loglik
    close all
    eval(['save ''z_loglik_' num2str(it_estim) '.mat'' estim_likelihood_logsum'])
end


%% now open all those estimation results and plot loglik graph

clear

sigma_grid =nan(11,1);
loglik_grid=nan(11,1);

for it_estim=130:2:150
    eval(['load ''z_loglik_' num2str(it_estim) '.mat'''])
    sigma_grid(it_estim/2-64)=it_estim/10000;
    loglik_grid(it_estim/2-64)=estim_likelihood_logsum;
end

lik_grid = exp(loglik_grid);


myfig=figure(5);
set(myfig, 'Position', [0 0 600 600])
plot(sigma_grid,loglik_grid,'-*b','linewidth',2);
hold on
plot([0.0140,0.0140],[-10000,+1000],'r--','linewidth',2);
title('Log-likelihood for different values of $\sigma$ (data: 1984-2017 GDP)', 'interpreter','latex','FontSize',14);
xlabel('$\sigma$','interpreter','latex','FontSize',12)
ylabel('log-likelihood','interpreter','latex','FontSize',12)
ylim([floor(min(loglik_grid*10))/10,ceil(max(loglik_grid*10))/10])
ax = gca;
ax.YAxis.Exponent = 0;

print -dpdf h05_loglik_sigma
savefig(myfig,'h05_loglik_sigma.fig');

title(' ', 'interpreter','latex','FontSize',14);
print -dpdf g05_loglik_sigma


load '../z_obs_data_bandpass.mat'

[maxval,maxpos]=max(loglik_grid);
loglik_deriv2 = (loglik_grid(maxpos-1)-2*loglik_grid(maxpos)+loglik_grid(maxpos+1))./((sigma_grid(maxpos)-sigma_grid(maxpos-1))^2);
lik_deriv2    = (   lik_grid(maxpos-1)-2*   lik_grid(maxpos)+   lik_grid(maxpos+1))./((sigma_grid(maxpos)-sigma_grid(maxpos-1))^2);

maxlik_stdev = 1/(sqrt(-loglik_deriv2*size(Y_obs,1)));

disp(' ')
disp(['Standard deviation of the maximum likelihood estimator: ' num2str(maxlik_stdev)])