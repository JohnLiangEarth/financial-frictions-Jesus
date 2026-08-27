% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

nval_IRF=100/dt;

shocksize = -2/sqrt(dt);

g_sss=g_hlsss;
N_sss=N_hlsss;

Birf     = zeros(nval_IRF,1);
BirfposU = zeros(nval_IRF,1);
BirfposD = zeros(nval_IRF,1);

Nirf     = zeros(nval_IRF,1);
NirfposU = zeros(nval_IRF,1);
NirfposD = zeros(nval_IRF,1);

cirf     = zeros(nval_a,nval_z,nval_IRF);
Cirf     = zeros(nval_IRF,1);
rirf     = zeros(nval_IRF,1);

girf     = zeros(nval_a,nval_z,nval_IRF);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% baseline scenario, at sss

e_irf = zeros(nval_IRF,1);
g1=g_sss;
Nirf(1)=N_sss;

c5_IRF_sim

irf_sssA_B = Birf;
irf_sssA_N = Nirf;
irf_sssA_r = rirf;
irf_sssA_g = girf;
irf_sssA_c = cirf;
irf_sssA_C = Cirf;

% alternative scenario with shock, at sss

e_irf(1) = shocksize;
g1=g_sss;
Nirf(1)=N_sss;

c5_IRF_sim

irf_sssB_B = Birf;
irf_sssB_N = Nirf;
irf_sssB_r = rirf;
irf_sssB_g = girf;
irf_sssB_c = cirf;
irf_sssB_C = Cirf;


% baseline scenario, at ss

e_irf = zeros(nval_IRF,1);
g1=g_ss;
Nirf(1)=N_ss;

c5_IRF_sim

irf_ssA_B = Birf;
irf_ssA_N = Nirf;
irf_ssA_r = rirf;
irf_ssA_g = girf;
irf_ssA_c = cirf;
irf_ssA_C = Cirf;

% alternative scenario with shock, at ss

e_irf(1) = shocksize;
g1=g_ss;
Nirf(1)=N_ss;

c5_IRF_sim

irf_ssB_B = Birf;
irf_ssB_N = Nirf;
irf_ssB_r = rirf;
irf_ssB_g = girf;
irf_ssB_c = cirf;
irf_ssB_C = Cirf;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% CALCULATE SECONDARY VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

irf_sssA_K     = irf_sssA_B+irf_sssA_N   ;
irf_sssA_Y     = Zeta * irf_sssA_K.^alpha   ;
irf_sssA_w     = (1-alpha)*irf_sssA_Y   ;
irf_sssA_rk    = alpha*irf_sssA_Y./irf_sssA_K - delta   ;
irf_sssA_Chat  = rhohat*irf_sssA_N   ;
irf_sssA_KF    = [irf_sssA_K(2:end) ; irf_sssA_K(end)];

irf_sssB_K     = irf_sssB_B+irf_sssB_N   ;
irf_sssB_Y     = Zeta * irf_sssB_K.^alpha   ;
irf_sssB_w     = (1-alpha)*irf_sssB_Y   ;
irf_sssB_rk    = alpha*irf_sssB_Y./irf_sssB_K - delta   ;
irf_sssB_Chat  = rhohat*irf_sssB_N   ;
irf_sssB_KF    = [irf_sssB_K(2:end) ; irf_sssB_K(end)];

irf_ssA_K     = irf_ssA_B+irf_ssA_N   ;
irf_ssA_Y     = Zeta * irf_ssA_K.^alpha   ;
irf_ssA_w     = (1-alpha)*irf_ssA_Y   ;
irf_ssA_rk    = alpha*irf_ssA_Y./irf_ssA_K - delta   ;
irf_ssA_Chat  = rhohat*irf_ssA_N   ;
irf_ssA_KF    = [irf_ssA_K(2:end) ; irf_ssA_K(end)];

irf_1lsssB_K     = irf_ssB_B+irf_ssB_N   ;
irf_ssB_Y     = Zeta * irf_1lsssB_K.^alpha   ;
irf_ssB_w     = (1-alpha)*irf_ssB_Y   ;
irf_ssB_rk    = alpha*irf_ssB_Y./irf_1lsssB_K - delta   ;
irf_ssB_Chat  = rhohat*irf_ssB_N   ;
irf_ssB_KF    = [irf_1lsssB_K(2:end) ; irf_1lsssB_K(end)];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE DIFFERENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IRF_sss_Y    = irf_sssB_Y    ./ irf_sssA_Y    *100-100;
IRF_sss_c    = irf_sssB_c    ./ irf_sssA_c    *100-100;
IRF_sss_C    = irf_sssB_C    ./ irf_sssA_C    *100-100;
IRF_sss_Chat = irf_sssB_Chat ./ irf_sssA_Chat *100-100;
IRF_sss_K    = irf_sssB_K    ./ irf_sssA_K    *100-100;
IRF_sss_B    = irf_sssB_B    ./ irf_sssA_B    *100-100;
IRF_sss_N    = irf_sssB_N    ./ irf_sssA_N    *100-100;
IRF_sss_w    = irf_sssB_w    ./ irf_sssA_w    *100-100;
IRF_sss_r    = irf_sssB_r    -  irf_sssA_r    ;
IRF_sss_rk   = irf_sssB_rk   -  irf_sssA_rk   ;
IRF_sss_g    = irf_sssB_g    -  irf_sssA_g    ;


IRF_ss_Y    = irf_ssB_Y    ./ irf_ssA_Y    *100-100;
IRF_ss_c    = irf_ssB_c    ./ irf_ssA_c    *100-100;
IRF_ss_C    = irf_ssB_C    ./ irf_ssA_C    *100-100;
IRF_ss_Chat = irf_ssB_Chat ./ irf_ssA_Chat *100-100;
IRF_ss_K    = irf_1lsssB_K    ./ irf_ssA_K    *100-100;
IRF_ss_B    = irf_ssB_B    ./ irf_ssA_B    *100-100;
IRF_ss_N    = irf_ssB_N    ./ irf_ssA_N    *100-100;
IRF_ss_w    = irf_ssB_w    ./ irf_ssA_w    *100-100;
IRF_ss_r    = irf_ssB_r    -  irf_ssA_r    ;
IRF_ss_rk   = irf_ssB_rk   -  irf_ssA_rk   ;
IRF_ss_g    = irf_ssB_g    -  irf_ssA_g    ;





%%
%%%%%%%%%%%%%%%
% plot graphs %
%%%%%%%%%%%%%%%

%%
% basic IRF


myfig=figure(50);
set(myfig, 'Position', [0 0 800 800])

subplot(2,2,1);
plot([1:nval_IRF]*dt, IRF_ss_C,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_sss_C,'-.','Color',[1,0.1,0.1],'linewidth',1);
title('(a) HH consumption, $c$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(2,2,2);
plot([1:nval_IRF]*dt, IRF_ss_Chat,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_sss_Chat,'-.','Color',[1,0.1,0.1],'linewidth',1);
title('(b) Expert consumption, $\hat{c}$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(2,2,3);
plot([1:nval_IRF]*dt, IRF_ss_B,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_sss_B,'-.','Color',[1,0.1,0.1],'linewidth',1);
title('(c) Debt, $B$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(2,2,4);
plot([1:nval_IRF]*dt, IRF_ss_N,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:nval_IRF]*dt, IRF_sss_N,'-.','Color',[1,0.1,0.1],'linewidth',1);
title('(d) Equity, $N$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on


legend({'At deterministic SS','At stochastic SS'},'Location','best', 'interpreter','latex','FontSize',10)
% legend boxoff


print -dpdf g50_IRF
savefig(myfig,'g50_IRF.fig');

print -dpdf g94_IRF
savefig(myfig,'g94_IRF.fig');


