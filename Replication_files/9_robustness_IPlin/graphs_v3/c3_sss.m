% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

tic
disp('Refining SSS points and plotting SSS graphs')
disp(' ')


sss_nval_sim=5000/dt;

sss_Bsim    = zeros(sss_nval_sim,1);
sss_BposD   = zeros(sss_nval_sim,1);
sss_BposU   = zeros(sss_nval_sim,1);

sss_Nsim    = zeros(sss_nval_sim,1);
sss_NposU   = zeros(sss_nval_sim,1);
sss_NposD   = zeros(sss_nval_sim,1);

sss_rsim    = zeros(sss_nval_sim,1);

sss_g_big   = zeros(nval_a,nval_z,sss_nval_sim+1);
sss_g_dif   = zeros(sss_nval_sim,1);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first calculate some more variables for the deterministic SS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_ss = B_ss + N_ss;
r_ss = alpha * Zeta * (K_ss.^(alpha-1)) - delta - sigma2*(K_ss./N_ss);

B_ss_posD =floor((B_ss-Bmin)/dB)+1;
B_ss_posU = ceil((B_ss-Bmin)/dB)+1;
wB=(B_grid(B_ss_posU)-B_ss)/dB;

N_ss_posD =floor((N_ss-Nmin)/dN)+1;
N_ss_posU = ceil((N_ss-Nmin)/dN)+1;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refine Baseline-SSS point %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mydistance = sqrt((Bsim-B_bsss).^2 + (Nsim-N_bsss).^2);
[mydistance_,myindex]=min(mydistance);

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


sss_Ksim=sss_Bsim+sss_Nsim;

sss_g_big_average = squeeze(sss_g_big(:,1,:) + sss_g_big(:,2,:));

disp('baseline SSS refinement')
disp(['B_bsss went from ' num2str(B_bsss) ' to ' num2str(sss_Bsim)])
disp(['N_bsss went from ' num2str(N_bsss) ' to ' num2str(sss_Nsim)])
disp(' ')

B_bsss=sss_Bsim;
N_bsss=sss_Nsim;
K_bsss=B_bsss+N_bsss;
g_bsss=g1;

r_bsss = sss_rsim;
c_DD=squeeze(c(:,:,sss_BposD,sss_NposD));
c_DU=squeeze(c(:,:,sss_BposD,sss_NposU));
c_UD=squeeze(c(:,:,sss_BposU,sss_NposD));
c_UU=squeeze(c(:,:,sss_BposU,sss_NposU));
c_bsss = wB*wN*c_DD + wB*(1-wN)*c_DU + (1-wB)*wN*c_UD + (1-wB)*(1-wN)*c_UU;
w_DD=squeeze(w(:,:,sss_BposD,sss_NposD));
w_DU=squeeze(w(:,:,sss_BposD,sss_NposU));
w_UD=squeeze(w(:,:,sss_BposU,sss_NposD));
w_UU=squeeze(w(:,:,sss_BposU,sss_NposU));
w_bsss = wB*wN*w_DD + wB*(1-wN)*w_DU + (1-wB)*wN*w_UD + (1-wB)*(1-wN)*w_UU;
s_bsss = w_bsss.*squeeze(z(:,:,1,1)) + r_bsss*squeeze(a(:,:,1,1)) - c_bsss;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refine LowLeverage-SSS point %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mydistance = sqrt((Bsim-B_llsss).^2 + (Nsim-N_llsss).^2);
[mydistance_,myindex]=min(mydistance);

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


sss_Ksim=sss_Bsim+sss_Nsim;

sss_g_big_average = squeeze(sss_g_big(:,1,:) + sss_g_big(:,2,:));

disp('low-leverage SSS refinement')
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

%%
%%%%%%%%%%%%%%%%%%%%
% plot SSS results %
%%%%%%%%%%%%%%%%%%%%


% distribution of assets 

myfig=figure(40);
set(myfig, 'Position', [0 0 600 400])

gsum_llsss = sum(g_llsss,2);
gsum_bsss  = sum(g_bsss,2);
gsum_ss    = sum(g_ss,2);

plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(gsum_llsss(1:ceil(size(g_ss,1)/2))),'--','Color',[1,0.1,0.1],'Linewidth',2)
hold on
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(gsum_bsss (1:ceil(size(g_ss,1)/2))),'-','Color',[0,0.5,0.1],'Linewidth',2)
plot(a_grid(1:ceil(size(g_ss,1)/2)),squeeze(gsum_ss   (1:ceil(size(g_ss,1)/2))),'-.','Color',[0.1,0.3,1],'Linewidth',2)
title('Distribution of assets ($g$)', 'interpreter','latex','FontSize',14);
xlabel('Assets ($a$)', 'interpreter','latex','FontSize',14);

legend({'Low-leverage stochastic steady state','Baseline stochastic steady state','Deterministic steady state',},'Location','northeast', 'interpreter','latex','FontSize',12)

grid

print -dpdf g40_SSS
print -depsc g40_SSS
savefig(myfig,'g40_SSS.fig');

toc
disp(' ')
