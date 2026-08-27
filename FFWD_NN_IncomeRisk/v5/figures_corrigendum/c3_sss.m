% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

sss_nval_sim=5000/dt;

tic

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first calculate some more variables for the deterministic SS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_ss = B_ss + N_ss;
r_ss = alpha * Zeta * (K_ss.^(alpha-1)) - delta - sigma2*(K_ss./N_ss);

B_ss_posD =floor((B_ss-Bmin)/dB)+1;
B_ss_posD = max(B_ss_posD,1);
B_ss_posD = min(B_ss_posD,nval_B-1);
B_ss_posU = B_ss_posD+1;
wB=(B_grid(B_ss_posU)-B_ss)/dB;

N_ss_posD =floor((N_ss-Nmin)/dN)+1;
N_ss_posD = max(N_ss_posD,1);
N_ss_posD = min(N_ss_posD,nval_N-1);
N_ss_posU = N_ss_posD+1;
wN=(N_grid(N_ss_posU)-N_ss)/dN;

c_DD=squeeze(c(:,:,B_ss_posD,N_ss_posD));
c_DU=squeeze(c(:,:,B_ss_posD,N_ss_posU));
c_UD=squeeze(c(:,:,B_ss_posU,N_ss_posD));
c_UU=squeeze(c(:,:,B_ss_posU,N_ss_posU));
c_ss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
w_DD=squeeze(w(:,:,B_ss_posD,N_ss_posD));
w_DU=squeeze(w(:,:,B_ss_posD,N_ss_posU));
w_UD=squeeze(w(:,:,B_ss_posU,N_ss_posD));
w_UU=squeeze(w(:,:,B_ss_posU,N_ss_posU));
w_ss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
s_ss = w_ss.*squeeze(z(:,:,1,1)) + r_ss*squeeze(a(:,:,1,1)) - c_ss;


%%
%%%%%%%%%%%%%%%%%%%%%%%
% refine HL-SSS point %
%%%%%%%%%%%%%%%%%%%%%%%

% find the point in the simulation that is closest to the HL-SSS, to start the simulation from there

mydistance = sqrt((Bsim-B_hlsss).^2 + (Nsim-N_hlsss).^2);
[mydistance_,myindex]=min(mydistance);
clear mydistance mydistance_

% now simulate without any shocks

for t=1:sss_nval_sim
    
    if t==1
        g1=squeeze(g_big(:,:,myindex));
        sss_Nsim=Nsim(myindex);
    else
        
        g0=g1;

        myA=A1{sss_BposD,sss_NposD};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposD,sss_NposU};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU,sss_NposD};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU,sss_NposU};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        g1 = wB*wN*g1DD + wB*(1-wN)*g1DU + (1-wB)*wN*g1UD + (1-wB)*(1-wN)*g1UU;

        sss_Nsim = sss_Nsim + ( alpha*Zeta*((sss_Bsim+sss_Nsim).^alpha) - delta*(sss_Bsim+sss_Nsim) - sss_rsim*sss_Bsim - rhohat*sss_Nsim )*dt;

    end

    sss_Bsim=sum(sum(g1.*a2*da));
    sss_Bsim=max([sss_Bsim Bmin]);
    sss_Bsim=min([sss_Bsim Bmax]);
    sss_BposD =floor((sss_Bsim-Bmin)/dB)+1;
    sss_BposU = ceil((sss_Bsim-Bmin)/dB)+1;
    wB=(B_grid(sss_BposU)-sss_Bsim)/dB;

    sss_Nsim=max([sss_Nsim Nmin]);
    sss_Nsim=min([sss_Nsim Nmax]);
    sss_NposD =floor((sss_Nsim-Nmin)/dN)+1;
    sss_NposU = ceil((sss_Nsim-Nmin)/dN)+1;
    wN=(N_grid(sss_NposU)-sss_Nsim)/dN;
    
    sss_rsim = alpha * Zeta * ((sss_Bsim+sss_Nsim).^(alpha-1)) - delta - sigma2*((sss_Bsim+sss_Nsim)./sss_Nsim);
    
end


disp('HL-SSS refinement')
disp(['B_hlsss went from ' num2str(B_hlsss) ' to ' num2str(sss_Bsim)])
disp(['N_hlsss went from ' num2str(N_hlsss) ' to ' num2str(sss_Nsim)])
disp(' ')

B_hlsss=sss_Bsim;
N_hlsss=sss_Nsim;
K_hlsss=B_hlsss+N_hlsss;
g_hlsss=g1;

r_hlsss = sss_rsim;
c_DD=squeeze(c(:,:,sss_BposD,sss_NposD));
c_DU=squeeze(c(:,:,sss_BposD,sss_NposU));
c_UD=squeeze(c(:,:,sss_BposU,sss_NposD));
c_UU=squeeze(c(:,:,sss_BposU,sss_NposU));
c_hlsss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
w_DD=squeeze(w(:,:,sss_BposD,sss_NposD));
w_DU=squeeze(w(:,:,sss_BposD,sss_NposU));
w_UD=squeeze(w(:,:,sss_BposU,sss_NposD));
w_UU=squeeze(w(:,:,sss_BposU,sss_NposU));
w_hlsss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
s_hlsss = w_hlsss.*squeeze(z(:,:,1,1)) + r_hlsss*squeeze(a(:,:,1,1)) - c_hlsss;
V_DD=squeeze(V1(:,:,sss_BposD,sss_NposD));
V_DU=squeeze(V1(:,:,sss_BposD,sss_NposU));
V_UD=squeeze(V1(:,:,sss_BposU,sss_NposD));
V_UU=squeeze(V1(:,:,sss_BposU,sss_NposU));
V1_hlsss = wB*wN*V_DD + wB*(1-wN)*V_DU + (1-wB)*wN*V_UD + (1-wB)*(1-wN)*V_UU;


%%
%%%%%%%%%%%%%%%%%%%%%%%
% refine LL-SSS point %
%%%%%%%%%%%%%%%%%%%%%%%

% find the point in the simulation that is closest to the LL-SSS, to start the simulation from there

mydistance = sqrt((Bsim-B_llsss).^2 + (Nsim-N_llsss).^2);
[mydistance_,myindex]=min(mydistance);
clear mydistance mydistance_

% now simulate without any shocks

for t=1:sss_nval_sim
    
    if t==1
        g1=squeeze(g_big(:,:,myindex));
        sss_Nsim=Nsim(myindex);
    else
        
        g0=g1;

        myA=A1{sss_BposD,sss_NposD};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposD,sss_NposU};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1DU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU,sss_NposD};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UD=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        myA=A1{sss_BposU,sss_NposU};
        g1_stacked=(speye(nval_a*nval_z)-myA'*dt) \g0(:);
        g1_stacked=g1_stacked/sum(g1_stacked*da);
        g1UU=[g1_stacked(1:nval_a) g1_stacked(nval_a+1:nval_a*2)];

        g1 = wB*wN*g1DD + wB*(1-wN)*g1DU + (1-wB)*wN*g1UD + (1-wB)*(1-wN)*g1UU;

        sss_Nsim = sss_Nsim + ( alpha*Zeta*((sss_Bsim+sss_Nsim).^alpha) - delta*(sss_Bsim+sss_Nsim) - sss_rsim*sss_Bsim - rhohat*sss_Nsim )*dt;

    end

    sss_Bsim=sum(sum(g1.*a2*da));
    sss_Bsim=max([sss_Bsim Bmin]);
    sss_Bsim=min([sss_Bsim Bmax]);
    sss_BposD =floor((sss_Bsim-Bmin)/dB)+1;
    sss_BposU = ceil((sss_Bsim-Bmin)/dB)+1;
    wB=(B_grid(sss_BposU)-sss_Bsim)/dB;

    sss_Nsim=max([sss_Nsim Nmin]);
    sss_Nsim=min([sss_Nsim Nmax]);
    sss_NposD =floor((sss_Nsim-Nmin)/dN)+1;
    sss_NposU = ceil((sss_Nsim-Nmin)/dN)+1;
    wN=(N_grid(sss_NposU)-sss_Nsim)/dN;
    
    sss_rsim = alpha * Zeta * ((sss_Bsim+sss_Nsim).^(alpha-1)) - delta - sigma2*((sss_Bsim+sss_Nsim)./sss_Nsim);
    
end


disp('HL-SSS refinement')
disp(['B_llsss went from ' num2str(B_llsss) ' to ' num2str(sss_Bsim)])
disp(['N_llsss went from ' num2str(N_llsss) ' to ' num2str(sss_Nsim)])
disp(' ')

B_llsss=sss_Bsim;
N_llsss=sss_Nsim;
K_llsss=B_llsss+N_llsss;
g_llsss=g1;

r_llsss = sss_rsim;
c_DD=squeeze(c(:,:,sss_BposD,sss_NposD));
c_DU=squeeze(c(:,:,sss_BposD,sss_NposU));
c_UD=squeeze(c(:,:,sss_BposU,sss_NposD));
c_UU=squeeze(c(:,:,sss_BposU,sss_NposU));
c_llsss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
w_DD=squeeze(w(:,:,sss_BposD,sss_NposD));
w_DU=squeeze(w(:,:,sss_BposD,sss_NposU));
w_UD=squeeze(w(:,:,sss_BposU,sss_NposD));
w_UU=squeeze(w(:,:,sss_BposU,sss_NposU));
w_llsss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
s_llsss = w_llsss.*squeeze(z(:,:,1,1)) + r_llsss*squeeze(a(:,:,1,1)) - c_llsss;
V_DD=squeeze(V1(:,:,sss_BposD,sss_NposD));
V_DU=squeeze(V1(:,:,sss_BposD,sss_NposU));
V_UD=squeeze(V1(:,:,sss_BposU,sss_NposD));
V_UU=squeeze(V1(:,:,sss_BposU,sss_NposU));
V1_llsss = wB*wN*V_DD + wB*(1-wN)*V_DU + (1-wB)*wN*V_UD + (1-wB)*(1-wN)*V_UU;


%%
%%%%%%%%%%%%%%%%%%%%
% plot SSS results %
%%%%%%%%%%%%%%%%%%%%

% report welfare

ss_BposD =floor((B_ss-Bmin)/dB)+1;
ss_BposU = ceil((B_ss-Bmin)/dB)+1;
wB=(B_grid(ss_BposU)-B_ss)/dB;

ss_NposD =floor((N_ss-Nmin)/dN)+1;
ss_NposU = ceil((N_ss-Nmin)/dN)+1;
wN=(N_grid(ss_NposU)-N_ss)/dN;

V_DD=squeeze(V1(:,:,ss_BposD,ss_NposD));
V_DU=squeeze(V1(:,:,ss_BposD,ss_NposU));
V_UD=squeeze(V1(:,:,ss_BposU,ss_NposD));
V_UU=squeeze(V1(:,:,ss_BposU,ss_NposU));
V1_ss = wB*wN*V_DD + wB*(1-wN)*V_DU + (1-wB)*wN*V_UD + (1-wB)*(1-wN)*V_UU;

disp('Average welfare at each SSS')
disp([ 'Average welfare at DSS:    ' num2str(sum(sum(V1_ss.*g_ss*da))) ])
disp([ 'Average welfare at LL-SSS: ' num2str(sum(sum(V1_llsss.*g_llsss*da))) ])
disp([ 'Average welfare at HL-SSS: ' num2str(sum(sum(V1_hlsss.*g_hlsss*da))) ])
disp(' ')

disp('Consumption equivalent of shift from DSS to each SSS')
disp([ 'LL-SSS: ' num2str( (sum(sum(V1_llsss.*g_llsss*da)) / sum(sum(V1_ss.*g_ss*da)))^(1/(1-gamma)) -1 ) ])
disp([ 'HL-SSS: ' num2str( (sum(sum(V1_hlsss.*g_hlsss*da)) / sum(sum(V1_ss.*g_ss*da)))^(1/(1-gamma)) -1 ) ])
disp(' ')

disp('Same, as percentage of aggregate consumption at each SSS')
disp([ 'LL-SSS: ' num2str( ((sum(sum(V1_llsss.*g_llsss*da)) / sum(sum(V1_ss.*g_ss*da)))^(1/(1-gamma)) -1 )/sum(sum(c_llsss.*g_llsss*da))*100) '%'])
disp([ 'HL-SSS: ' num2str( ((sum(sum(V1_hlsss.*g_hlsss*da)) / sum(sum(V1_ss.*g_ss*da)))^(1/(1-gamma)) -1 )/sum(sum(c_hlsss.*g_hlsss*da))*100) '%'])
disp(' ')


% distribution of assets 

myfig=figure(40);
set(myfig, 'Position', [0 0 800 400])

subplot(1,2,1);
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_llsss(2:ceil(size(g_ss,1)/2),1)),'-','Color',[1,0.1,0.1],'Linewidth',2)
hold on
% plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_hlsss(2:ceil(size(g_ss,1)/2),1)),'-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(2:ceil(size(g_ss,1)/2)),squeeze(g_ss   (2:ceil(size(g_ss,1)/2),1)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
plot(a_grid(1),squeeze(g_llsss(1,1)),'o','Color',[1,0.1,0.1],'Linewidth',2)
% plot(a_grid(1),squeeze(g_hlsss(1,1)),'o','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1),squeeze(g_ss   (1,1)),'o','Color',[0.1,0.3,1],'Linewidth',2)
title('(a) Low-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

subplot(1,2,2);
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_llsss(1:ceil(size(g_ss,1)/2),2)),'-','Color',[1,0.1,0.1],'Linewidth',2)
hold on
% plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_hlsss(1:ceil(size(g_ss,1)/2),2)),'-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(g_ss   (1:ceil(size(g_ss,1)/2),2)),'-.','Color',[0.1,0.3,1],'Linewidth',2)
title('(b) High-$z$ households', 'interpreter','latex','FontSize',14);
xlabel('assets ($a$)', 'interpreter','latex','FontSize',14);
grid

legend({'Stochastic SS','Deterministic SS',},'Location','northeast', 'interpreter','latex','FontSize',12)

print -dpdf g40_SSS
savefig(myfig,'g40_SSS.fig');

print -dpdf g93_SSS
savefig(myfig,'g93_SSS.fig');

disp('Time refining SSS points and plotting SSS graphs:')
toc
disp(' ')

%%




disp('Gini coefficients:')

g0=ones(size(g_hlsss(:,1)))/nval_a/da; % uniform, just for testing my Gini calculations

% Gini, method 1
gini_1 = g0*da;
gini_2 = repmat(gini_1,1,nval_a);
gini_3 = repmat(a_grid,1,nval_a);
gini_4 = gini_2 .* gini_2';
gini_5 = abs(gini_3 - gini_3');
g_gini = 1 / (2*sum(gini_1.*a_grid)) * sum(sum( gini_4.*gini_5 ));     % zero means perfect equality, one means total inequality

% Gini, method 2

gini_1 = g0*da;
gini_2 = g0 .* a_grid;
gini_3 = cumsum(gini_2);
gini_4 = gini_3 / gini_3(end);
gini_4(2:end)=(gini_4(1:end-1)+gini_4(2:end))/2;
gini_5 = gini_1 .* gini_4;
g_gini_= 1-2*sum(gini_5);

disp(['test (uniform, should be 1/3):    ' num2str(g_gini) '    '  num2str(g_gini_)])

myfig = figure(41);
plot(cumsum(gini_1),gini_4)
hold on



g0 = g_hlsss(:,1)+g_hlsss(:,2);

% Gini, method 1
gini_1 = g0*da;
gini_2 = repmat(gini_1,1,nval_a);
gini_3 = repmat(a_grid,1,nval_a);
gini_4 = gini_2 .* gini_2';
gini_5 = abs(gini_3 - gini_3');
g_gini = 1 / (2*sum(gini_1.*a_grid)) * sum(sum( gini_4.*gini_5 ));     % zero means perfect equality, one means total inequality

% Gini, method 2

gini_1 = g0*da;
gini_2 = g0 .* a_grid;
gini_3 = cumsum(gini_2);
gini_4 = gini_3 / gini_3(end);
gini_4(2:end)=(gini_4(1:end-1)+gini_4(2:end))/2;
gini_5 = gini_1 .* gini_4;
g_gini_= 1-2*sum(gini_5);

disp(['HL-SSS:    ' num2str(g_gini) '    '  num2str(g_gini_)])

plot(cumsum(gini_1),gini_4)



g0 = g_llsss(:,1) + g_llsss(:,2);

% Gini, method 1
gini_1 = g0*da;
gini_2 = repmat(gini_1,1,nval_a);
gini_3 = repmat(a_grid,1,nval_a);
gini_4 = gini_2 .* gini_2';
gini_5 = abs(gini_3 - gini_3');
g_gini = 1 / (2*sum(gini_1.*a_grid)) * sum(sum( gini_4.*gini_5 ));     % zero means perfect equality, one means total inequality

% Gini, method 2

gini_1 = g0*da;
gini_2 = g0 .* a_grid;
gini_3 = cumsum(gini_2);
gini_4 = gini_3 / gini_3(end);
gini_4(2:end)=(gini_4(1:end-1)+gini_4(2:end))/2;
gini_5 = gini_1 .* gini_4;
g_gini_= 1-2*sum(gini_5);

disp(['LL-SSS:    ' num2str(g_gini) '    '  num2str(g_gini_)])

plot(cumsum(gini_1),gini_4)



g0 = g_ss(:,1) + g_ss(:,2);

% Gini, method 1
gini_1 = g0*da;
gini_2 = repmat(gini_1,1,nval_a);
gini_3 = repmat(a_grid,1,nval_a);
gini_4 = gini_2 .* gini_2';
gini_5 = abs(gini_3 - gini_3');
g_gini = 1 / (2*sum(gini_1.*a_grid)) * sum(sum( gini_4.*gini_5 ));     % zero means perfect equality, one means total inequality

% Gini, method 2

gini_1 = g0*da;
gini_2 = g0 .* a_grid;
gini_3 = cumsum(gini_2);
gini_4 = gini_3 / gini_3(end);
gini_4(2:end)=(gini_4(1:end-1)+gini_4(2:end))/2;
gini_5 = gini_1 .* gini_4;
g_gini_= 1-2*sum(gini_5);

disp(['DSS:    ' num2str(g_gini) '    '  num2str(g_gini_)])

plot(cumsum(gini_1),gini_4)
title('Lorenz curves (used for calculating Gini coefficients)', 'interpreter','latex','FontSize',14);
legend({'Uniform','HL-SSS','LL-SSS','DSS',},'Location','northwest', 'interpreter','latex','FontSize',12)


%% variance of a

myfig=figure(43);
hold on

% at the DSS

g0 = g_ss(:,1)+g_ss(:,2);
plot(a_grid,g0,'-.','Color',[0.1,0.3,1],'Linewidth',2);
g0_chck=sum(g0*da); % comprobar que suma uno
g0_mean=sum(g0.*a_grid*da);
g0_vari=sum(g0.*a_grid.*a_grid*da)-sum(g0.*a_grid*da)^2;
g0_stdv=g0_vari^0.5;

amean_dss=g0_mean;
astdv_dss=g0_stdv


% at the HLSSS

g0 = g_hlsss(:,1)+g_hlsss(:,2);
plot(a_grid,g0,'-','Color',[0,0.5,0.1],'Linewidth',2);
g0_chck=sum(g0*da); % comprobar que suma uno
g0_mean=sum(g0.*a_grid*da);
g0_vari=sum(g0.*a_grid.*a_grid*da)-sum(g0.*a_grid*da)^2;
g0_stdv=g0_vari^0.5;

amean_hlsss=g0_mean;
astdv_hlsss=g0_stdv


% at the LLSSS

g0 = g_llsss(:,1)+g_llsss(:,2);
plot(a_grid,g0,'--','Color',[1,0.1,0.1],'Linewidth',2);
g0_chck=sum(g0*da); % comprobar que suma uno
g0_mean=sum(g0.*a_grid*da);
g0_vari=sum(g0.*a_grid.*a_grid*da)-sum(g0.*a_grid*da)^2;
g0_stdv=g0_vari^0.5;

amean_llsss=g0_mean;
astdv_llsss=g0_stdv


% over the full simulation

g00=mean(g_big,3);
g0 = g00(:,1)+g00(:,2);
plot(a_grid,g0,'-','Color',[0,0,0],'Linewidth',2);
g0_chck=sum(g0*da); % comprobar que suma uno
g0_mean=sum(g0.*a_grid*da);
g0_vari=sum(g0.*a_grid.*a_grid*da)-sum(g0.*a_grid*da)^2;
g0_stdv=g0_vari^0.5;

amean_sim=g0_mean;
astdv_sim=g0_stdv



legend({'Deterministic SS','HL-SSS','LL-SSS','simulation'},'Location','northeast', 'interpreter','latex','FontSize',12)

xlim([0 8])
print -dpdf g43_adist
savefig(myfig,'g43_adist.fig');



%% variance of w*z

% first just check that the unemployment rate is constant, then save it
g_big_check = squeeze(sum(g_big,1));
g_big_check(1,:)=g_big_check(1,:)-sum(g_ss(:,1));
g_big_check(2,:)=g_big_check(2,:)-sum(g_ss(:,2));
max(max(g_big_check));
min(min(g_big_check));
clear g_big_check;

UR=sum(g_ss(:,1))*da;

% at the DSS
w_ss_=w_ss(1,1);
wzmean_dss=(w_ss_*z1)*UR+(w_ss_*z2)*(1-UR);
wzvari_dss=((w_ss_*z1)^2)*UR+((w_ss_*z2)^2)*(1-UR) - (wzmean_dss^2);
wzstdv_dss=wzvari_dss^0.5

% over the full simulation
wsim = (1-alpha)*Ysim;
nsim = size(wsim,1);
wzmean_sim=sum( (wsim*z1)    *UR+ (wsim*z2)    *(1-UR))/nsim;
wzvari_sim=sum(((wsim*z1).^2)*UR+((wsim*z2).^2)*(1-UR))/nsim - (wzmean_sim^2);
wzstdv_sim=wzvari_sim^0.5


%% variance of r*a

% at the DSS
g0 = g_ss(:,1)+g_ss(:,2);
ramean_dss=sum(g0*r_ss     .*a_grid        *da);
ravari_dss=sum(g0*r_ss*r_ss.*a_grid.*a_grid*da)-ramean_dss^2;
rastdv_dss=ravari_dss^0.5

% over the full simulation

g0=squeeze(sum(g_big,2));
ramean_sim=sum(sum(g0.*repmat(rsim',nval_a,1)                        .*repmat(a_grid,1,nsim)                       *da))/nsim;
ravari_sim=sum(sum(g0.*repmat(rsim',nval_a,1).*repmat(rsim',nval_a,1).*repmat(a_grid,1,nsim).*repmat(a_grid,1,nsim)*da))/nsim-ramean_sim^2;
rastdv_sim=ravari_sim^0.5

wstdv_sim = (sum(wsim.^2)/nsim - (sum(wsim)/nsim)^2)^0.5
rstdv_sim = (sum(rsim.^2)/nsim - (sum(rsim)/nsim)^2)^0.5


%% variance of total income, wzra=w*z+r*a

% at the DSS
wzra_dss = repmat(w_ss_,nval_a,2).*repmat([z1 z2],nval_a,1) + r_ss.*repmat(a_grid,1,2);
wzramean_dss=sum(sum( wzra_dss    .*g_ss*da));
wzravari_dss=sum(sum((wzra_dss.^2).*g_ss*da)) - (wzramean_dss^2);
wzrastdv_dss=wzravari_dss^0.5

% over the full simulation (first transform g_big from 501*2*nsim to 1002*nsim)
g_big_=[squeeze(g_big(:,1,:)) ; squeeze(g_big(:,2,:))];
wzra_sim = repmat(wsim',nval_a*2,1).*repmat([repmat(z1,nval_a,1);repmat(z2,nval_a,1)],1,nsim) + repmat(rsim',nval_a*2,1).*repmat([a_grid;a_grid],1,nsim);

wzra_sim1 = repmat(rsim',nval_a*2,1).*repmat([a_grid;a_grid],1,nsim);
wzra_sim2 = repmat(wsim',nval_a*2,1).*repmat([repmat(z1,nval_a,1);repmat(z2,nval_a,1)],1,nsim) ;

wzramean_sim=sum(sum( wzra_sim    .*g_big_))*da/nsim;
wzravari_sim=sum(sum((wzra_sim.^2).*g_big_))*da/nsim - (wzramean_sim^2);
wzrastdv_sim=wzravari_sim^0.5

wzracov_sim = sum(sum((wzra_sim1.*wzra_sim2).*g_big_))*da/nsim - (wzmean_sim*ramean_sim);
wzravari_check = ravari_sim + wzvari_sim + 2* wzracov_sim;


%% variance of individual consumption c

% at the DSS
cmean_dss=sum(sum( c_ss    .*g_ss*da));
cvari_dss=sum(sum((c_ss.^2).*g_ss*da)) - (cmean_dss^2);
cstdv_dss=cvari_dss^0.5

% over the full simulation (first I have to calculate c(:,:,Bt,Nt) for each period)

csim=nan(2*nval_a,nsim);
for t=1:nsim
    BposD = floor((Bsim(t)-Bmin)/dB)+1;
    BposU =  ceil((Bsim(t)-Bmin)/dB)+1;
    wB=(B_grid(BposU)-Bsim(t))/dB;

    NposD = floor((Nsim(t)-Nmin)/dN)+1;
    NposU =  ceil((Nsim(t)-Nmin)/dN)+1;
    wN=(N_grid(NposU)-Nsim(t))/dN;

    c1 = wB*wN*c(:,1,BposD,NposD) + wB*(1-wN)*c(:,1,BposD,NposU) + (1-wB)*wN*c(:,1,BposU,NposD) + (1-wB)*(1-wN)*c(:,1,BposU,NposU);
    c2 = wB*wN*c(:,2,BposD,NposD) + wB*(1-wN)*c(:,2,BposD,NposU) + (1-wB)*wN*c(:,2,BposU,NposD) + (1-wB)*(1-wN)*c(:,2,BposU,NposU);
    
    csim(:,t) = [c1 ; c2];
end

cmean_sim=sum(sum( csim    .*g_big_))*da/nsim;
cvari_sim=sum(sum((csim.^2).*g_big_))*da/nsim - (cmean_sim^2);
cstdv_sim=cvari_sim^0.5

clear g_big_ wzra_sim wzra_dss csim c1 c2 BposD BposU NposD NposU wB wN;


%% display

disp(' ')
disp(' ')
disp(' ')
disp('standard deviation of individual variables')
disp(' first row is "at DSS"')
disp('second row is "over the whole simulation"')
disp(' ')
disp('       a        w*z       r*a    w*z+r*a       c')
disp([ astdv_dss wzstdv_dss rastdv_dss wzrastdv_dss cstdv_dss])
disp([ astdv_sim wzstdv_sim rastdv_sim wzrastdv_sim cstdv_sim])
disp(' ')
disp(' ')
disp('ratio of standard deviations of individual variables')
disp('over the whole simulation / at the DSS')
disp('       a        w*z       r*a    w*z+r*a       c')
disp([ astdv_sim/astdv_dss wzstdv_sim/wzstdv_dss rastdv_sim/rastdv_dss wzrastdv_sim/wzrastdv_dss cstdv_sim/cstdv_dss])



