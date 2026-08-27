% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

copyfile '../f2_NN_eval.m';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use NN1 to calculate the PLM on the 101x101 grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X1m=zeros(nval_BB*nval_NN,1);
X2m=zeros(nval_BB*nval_NN,1);

for iN=1:nval_NN
    for iB=1:nval_BB
        X1m(iN*nval_BB-nval_BB+iB,1)=BB_grid(iB);
        X2m(iN*nval_BB-nval_BB+iB,1)=NN_grid(iN);
    end
end

Yr_mat2 = F(X1m,X2m);
PLM_surf= reshape(Yr_mat2,[nval_BB nval_NN]);

% visited part of that PLM
PLM_surf_v = PLM_surf;
for iN=1:nval_NN
    for iB=1:nval_BB
       if isnan(PLM_visits(iB,iN))
           PLM_surf_v(iB,iN)=NaN;      % put NaN on non-visited points
       end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the equation for N, without shocks, to calculate the N-PLM on the 101x101 grid %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_surf  =   alpha*Zeta*((BB_grid_2D+NN_grid_2D).^(alpha-1)) - delta - sigma2*((BB_grid_2D+NN_grid_2D)./NN_grid_2D);
NPLM_surf = alpha*Zeta*((BB_grid_2D+NN_grid_2D).^alpha)     - delta*(BB_grid_2D+NN_grid_2D) - r_surf.*BB_grid_2D - rhohat*NN_grid_2D ;

% visited part of that NPLM
NPLM_surf_v = NPLM_surf;
for iN=1:nval_NN
    for iB=1:nval_BB
       if isnan(PLM_visits(iB,iN))
           NPLM_surf_v(iB,iN)=NaN;     % put NaN on non-visited points
       end
    end
end

%%
%%%%%%%%%%%%%
% find PLM0 %
%%%%%%%%%%%%%

PLM0=zeros(nval_BB,2);

for it1=1:nval_BB
    it3=1;
    for it2=1:nval_NN
        if PLM_surf(it1,it2)>0
            it3=it2;
        end
    end
    if it3>=nval_NN
        it3=nval_NN-1;
    end
    myweight=PLM_surf(it1,it3+1)/(PLM_surf(it1,it3+1)-PLM_surf(it1,it3));
    PLM0(it1,1)=BB_grid(it1);
    PLM0(it1,2)=NN_grid(it3)*myweight+NN_grid(it3+1)*(1-myweight);
end


%%%%%%%%%%%%%%
% find NPLM0 %
%%%%%%%%%%%%%%

NPLM0=zeros(nval_BB,2);

for it1=1:nval_BB
    it3=1;
    for it2=1:nval_NN
        if NPLM_surf(it1,it2)>0
            it3=it2;
        end
    end
    if it3>=nval_NN
        it3=nval_NN-1;
    end
    myweight=NPLM_surf(it1,it3+1)/(NPLM_surf(it1,it3+1)-NPLM_surf(it1,it3));
    NPLM0(it1,1)=BB_grid(it1);
    NPLM0(it1,2)=NN_grid(it3)*myweight+NN_grid(it3+1)*(1-myweight);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find all Stochastic Steady States %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find SSS points (it's just for graphs, so: using linear approximation around zero of PLM and NPLM surfaces)
% I'll use a very long no-shocks simulation to refine these points, before using them to run IRFs etc
PLM0dif=PLM0;
PLM0dif(:,2)=PLM0(:,2)-NPLM0(:,2);
SSS_points=[];
for it1=2:size(PLM0,1)
    if (PLM0dif(it1,2)*PLM0dif(it1-1,2))<0
        myweight=PLM0dif(it1,2)/(PLM0dif(it1,2)-PLM0dif(it1-1,2)); % weight of the it1-1 point
        SSS_B = myweight*PLM0dif(it1-1,1)+(1-myweight)*PLM0dif(it1,1);
        SSS_N1 = myweight* PLM0(it1-1,2)+(1-myweight)* PLM0(it1,2);
        SSS_N2 = myweight*NPLM0(it1-1,2)+(1-myweight)*NPLM0(it1,2);
        SSS_points=[SSS_points ; [SSS_B (SSS_N1+SSS_N2)/2] ];
    end
end

% fill BSSS, LLSSS and NSSSS (baseline, high-leverage and non-sable SSS) using the values we just found

if size(SSS_points,1)==1       % if there's only one, it's a stable SSS
    B_llsss = SSS_points(1,1);
    N_llsss = SSS_points(1,2);
    B_nssss = SSS_points(1,1);
    N_nssss = SSS_points(1,2);
    B_bsss  = SSS_points(1,1);
    N_bsss  = SSS_points(1,2);
elseif size(SSS_points,1)==2   % if there are two, it's one stable and one non-stable SSS
    B_llsss = SSS_points(2,1);
    N_llsss = SSS_points(2,2);
    B_nssss = SSS_points(1,1);
    N_nssss = SSS_points(1,2);
    B_bsss  = SSS_points(2,1);
    N_bsss  = SSS_points(2,2);
else                           % if there are three, then we have ll, ns and b SSS
    B_llsss = SSS_points(1,1);
    N_llsss = SSS_points(1,2);
    B_nssss = SSS_points(2,1);
    N_nssss = SSS_points(2,2);
    B_bsss  = SSS_points(3,1);
    N_bsss  = SSS_points(3,2);
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define special points that will be used in the plots %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B1 = B_llsss;
N1 = N_llsss;

B2 = B_bsss;
N2 = N_bsss;

B3 = 2.20;         % this one is an arbitrary high-leverage point
N3 = 1.50;


mywidth=0.01;

B1_Nmin=inf;
B1_Nmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Bsim(t)-B1)<mywidth
        B1_Nmin=min([B1_Nmin Nsim(t)]);
        B1_Nmax=max([B1_Nmax Nsim(t)]);
    end
end

B2_Nmin=inf;
B2_Nmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Bsim(t)-B2)<mywidth
        B2_Nmin=min([B2_Nmin Nsim(t)]);
        B2_Nmax=max([B2_Nmax Nsim(t)]);
    end
end

B3_Nmin=inf;
B3_Nmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Bsim(t)-B3)<mywidth
        B3_Nmin=min([B3_Nmin Nsim(t)]);
        B3_Nmax=max([B3_Nmax Nsim(t)]);
    end
end

if B3_Nmin==inf
    B3 = B2;
    N3 = N2;
    B3_Nmin=B2_Nmin;
    B3_Nmax=B2_Nmax;
end

N1_Bmin=inf;
N1_Bmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Nsim(t)-N1)<mywidth
        N1_Bmin=min([N1_Bmin Bsim(t)]);
        N1_Bmax=max([N1_Bmax Bsim(t)]);
    end
end

N2_Bmin=inf;
N2_Bmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Nsim(t)-N2)<mywidth
        N2_Bmin=min([N2_Bmin Bsim(t)]);
        N2_Bmax=max([N2_Bmax Bsim(t)]);
    end
end

N3_Bmin=inf;
N3_Bmax=-inf;

for it1=1:size(nsim_used,1)
    t=nsim_used(it1);
    if abs(Nsim(t)-N3)<mywidth
        N3_Bmin=min([N3_Bmin Bsim(t)]);
        N3_Bmax=max([N3_Bmax Bsim(t)]);
    end
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use NN1 to calculate the PLM on several specific cuts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cut at B1/B2/B3, full range

X1m = B1*ones(size(BB_grid));
X2m = NN_grid;

Yr_mat2 = F(X1m,X2m);
PLM_B1full = Yr_mat2;


X1m = B2*ones(size(BB_grid));
Yr_mat2 = F(X1m,X2m);
PLM_B2full = Yr_mat2;


X1m = B3*ones(size(BB_grid));
Yr_mat2 = F(X1m,X2m);
PLM_B3full = Yr_mat2;


% cut at B1, visited range

X1m = B1*ones(101,1);
X2m = B1_Nmin:(B1_Nmax-B1_Nmin)/100:B1_Nmax;
X2m=X2m';
Yr_mat2 = F(X1m,X2m);
PLM_B1visit = Yr_mat2;


% cut at B2, visited range

X1m = B2*ones(101,1);
X2m = B2_Nmin:(B2_Nmax-B2_Nmin)/100:B2_Nmax;
X2m=X2m';
Yr_mat2 = F(X1m,X2m);
PLM_B2visit = Yr_mat2;


% cut at B3, visited range

X1m = B3*ones(101,1);
X2m = B3_Nmin:(B3_Nmax-B3_Nmin)/100:B3_Nmax;
X2m=X2m';
Yr_mat2 = F(X1m,X2m);
PLM_B3visit = Yr_mat2;


% cut at N1/N2/N3, full range

X1m = BB_grid;
X2m = N1*ones(size(NN_grid));
Yr_mat2 = F(X1m,X2m);
PLM_N1full = Yr_mat2;


X2m = N2*ones(size(NN_grid));
Yr_mat2 = F(X1m,X2m);
PLM_N2full = Yr_mat2;


X2m = N3*ones(size(NN_grid));
Yr_mat2 = F(X1m,X2m);
PLM_N3full = Yr_mat2;


% cut at N1, visited range

X1m = N1_Bmin:(N1_Bmax-N1_Bmin)/100:N1_Bmax;
X2m = N1*ones(101,1);
X1m=X1m';
Yr_mat2 = F(X1m,X2m);
PLM_N1visit = Yr_mat2;


% cut at N2, visited range

X1m = N2_Bmin:(N2_Bmax-N2_Bmin)/100:N2_Bmax;
X2m = N2*ones(101,1);
X1m=X1m';
Yr_mat2 = F(X1m,X2m);
PLM_N2visit = Yr_mat2;


% cut at N3, visited range

X1m = N3_Bmin:(N3_Bmax-N3_Bmin)/100:N3_Bmax;
X2m = N3*ones(101,1);
X1m=X1m';
Yr_mat2 = F(X1m,X2m);
PLM_N3visit = Yr_mat2;


% just the three points

X1m = [B1;B2;B3];
X2m = [N1;N2;N3];
Yr_mat2 = F(X1m,X2m);
PLM_123 = Yr_mat2;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the equation for N, without shocks, to calculate the N-PLM on several specific cuts %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cut at B1/B2/B3, full range

NN_range = NN_grid;

BB_range = B1*ones(size(BB_grid));
r_surf  =      alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_B1full = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;

BB_range = B2*ones(size(BB_grid));
r_surf  =        alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_B2full = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;

BB_range = B3*ones(size(BB_grid));
r_surf  =      alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_B3full = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;


% cut at B1, visited range

BB_range = B1*ones(101,1);
NN_range = B1_Nmin:(B1_Nmax-B1_Nmin)/100:B1_Nmax;
NN_range = NN_range';

r_surf =       alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_B1visit = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;
B1_Nvisit = NN_range;


% cut at B2, visited range

BB_range = B2*ones(101,1);
NN_range = B2_Nmin:(B2_Nmax-B2_Nmin)/100:B2_Nmax;
NN_range = NN_range';

r_surf =       alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_B2visit = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;
B2_Nvisit = NN_range;


% cut at B3, visited range

BB_range = B3*ones(101,1);
NN_range = B3_Nmin:(B3_Nmax-B3_Nmin)/100:B3_Nmax;
NN_range = NN_range';

r_surf  =       alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_B3visit = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;
B3_Nvisit = NN_range;


% cut at N1/N2/N3, full range

BB_range = BB_grid;

NN_range = N1*ones(size(NN_grid));
r_surf =      alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_N1full = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;


NN_range = N2*ones(size(NN_grid));
r_surf =      alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_N2full = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;


NN_range = N3*ones(size(NN_grid));
r_surf =      alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_N3full = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;


% cut at N1, visited range

BB_range = N1_Bmin:(N1_Bmax-N1_Bmin)/100:N1_Bmax;
NN_range = N1*ones(101,1);
BB_range = BB_range';

r_surf =       alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_N1visit = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;
N1_Bvisit = BB_range;


% cut at N2, visited range

BB_range = N2_Bmin:(N2_Bmax-N2_Bmin)/100:N2_Bmax;
NN_range = N2*ones(101,1);
BB_range = BB_range';

r_surf =       alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_N2visit = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;
N2_Bvisit = BB_range;


% cut at N3, visited range

BB_range = N3_Bmin:(N3_Bmax-N3_Bmin)/100:N3_Bmax;
NN_range = N3*ones(101,1);
BB_range = BB_range';

r_surf =       alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_N3visit = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;
N3_Bvisit = BB_range;


% just the three points

BB_range = [B1;B2;B3];
NN_range = [N1;N2;N3];

r_surf =   alpha*Zeta*((BB_range+NN_range).^(alpha-1)) - delta - sigma2*((BB_range+NN_range)./NN_range);
NPLM_123 = alpha*Zeta*((BB_range+NN_range).^alpha)     - delta*(BB_range+NN_range) - r_surf.*BB_range - rhohat*NN_range ;


%%
%%%%%%%%%%%%%%%%%%%%%%%
% Now plot everything %
%%%%%%%%%%%%%%%%%%%%%%%


% PLM: h(B,N)

myfig=figure(10);
set(myfig, 'Position', [0 0 800 1000])

subplot(3,2,1);
plot(NN_grid,PLM_B1full,'-','Color',[1,0.1,0.1]);
hold on;
plot(NN_grid,PLM_B2full,'-','Color',[0,0.5,0.1]);
plot(NN_grid,PLM_B3full,'-','Color',[0.1,0.3,1]);
plot(B1_Nvisit,PLM_B1visit,'-','Color',[1,0.1,0.1],'linewidth',3);
plot(B2_Nvisit,PLM_B2visit,'-','Color',[0,0.5,0.1],'linewidth',3);
plot(B3_Nvisit,PLM_B3visit,'-','Color',[0.1,0.3,1],'linewidth',3);
plot([N1;N2;N3],PLM_123,'.w');
title('(a) $h(B,N)$ for different values of $B$', 'interpreter','latex','FontSize',14);
xlabel('equity ($N$)', 'interpreter','latex','FontSize',14);
ylabel('$h(B,N)$', 'interpreter','latex','FontSize',14);
legend({'B at LL-SSS value','B at baseline SSS value','B at high-leverage point'},'Location','northeast', 'interpreter','latex','FontSize',10)

subplot(3,2,2);
plot(BB_grid,PLM_N1full,'-','Color',[1,0.1,0.1]);
hold on;
plot(BB_grid,PLM_N2full,'-','Color',[0,0.5,0.1]);
plot(BB_grid,PLM_N3full,'-','Color',[0.1,0.3,1]);
plot(N1_Bvisit,PLM_N1visit,'-','Color',[1,0.1,0.1],'linewidth',3);
plot(N2_Bvisit,PLM_N2visit,'-','Color',[0,0.5,0.1],'linewidth',3);
plot(N3_Bvisit,PLM_N3visit,'-','Color',[0.1,0.3,1],'linewidth',3);
plot([B1;B2;B3],PLM_123,'.w');
title('(b) $h(B,N)$ for different values of $N$', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('$h(B,N)$', 'interpreter','latex','FontSize',14);
legend({'N at LL-SSS value','N at baseline SSS value','N at high-leverage point'},'Location','northeast', 'interpreter','latex','FontSize',10)

subplot(3,2,[3:6]);
mesh(BB_grid,NN_grid,PLM_surf');
hold on
surf(BB_grid,NN_grid,PLM_surf_v');
plot(PLM0(:,1),PLM0(:,2),'-','Color',[1,0.1,0.1]);
title('(c) The perceived law of motion, $h(B,N)$', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
ylim([Nmin Nmax])
view([60, 15])

print -dpdf g10_PLM
print -depsc g10_PLM
savefig(myfig,'g10_PLM.fig');


%%
% NPLM: dN(B,N)

myfig=figure(15);
set(myfig, 'Position', [0 0 800 1000])

subplot(3,2,1);
plot(NN_grid,NPLM_B1full,'-','Color',[1,0.1,0.1]);
hold on;
plot(NN_grid,NPLM_B2full,'-','Color',[0,0.5,0.1]);
plot(NN_grid,NPLM_B3full,'-','Color',[0.1,0.3,1]);
plot(B1_Nvisit,NPLM_B1visit,'-','Color',[1,0.1,0.1],'linewidth',3);
plot(B2_Nvisit,NPLM_B2visit,'-','Color',[0,0.5,0.1],'linewidth',3);
plot(B3_Nvisit,NPLM_B3visit,'-','Color',[0.1,0.3,1],'linewidth',3);
plot([N1;N2;N3],NPLM_123,'.w');
title('(a) $dN(B,N)$ for different values of $B$', 'interpreter','latex','FontSize',14);
xlabel('equity ($N$)', 'interpreter','latex','FontSize',14);
ylabel('$dN(B,N)$', 'interpreter','latex','FontSize',14);
legend({'B at LL-SSS value','B at baseline SSS value','B at high-leverage point'},'Location','northeast', 'interpreter','latex','FontSize',10)

subplot(3,2,2);
plot(BB_grid,NPLM_N1full,'-','Color',[1,0.1,0.1]);
hold on;
plot(BB_grid,NPLM_N2full,'-','Color',[0,0.5,0.1]);
plot(BB_grid,NPLM_N3full,'-','Color',[0.1,0.3,1]);
plot(N1_Bvisit,NPLM_N1visit,'-','Color',[1,0.1,0.1],'linewidth',3);
plot(N2_Bvisit,NPLM_N2visit,'-','Color',[0,0.5,0.1],'linewidth',3);
plot(N3_Bvisit,NPLM_N3visit,'-','Color',[0.1,0.3,1],'linewidth',3);
plot([B1;B2;B3],NPLM_123,'.w');
title('(b) $dN(B,N)$ for different values of $N$', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('$dN(B,N)$', 'interpreter','latex','FontSize',14);
legend({'N at LL-SSS value','N at baseline SSS value','N at high-leverage point'},'Location','northeast', 'interpreter','latex','FontSize',10)

subplot(3,2,[3:6]);
mesh(BB_grid,NN_grid,NPLM_surf');
hold on
surf(BB_grid,NN_grid,NPLM_surf_v');
plot(NPLM0(:,1),NPLM0(:,2),'-','Color',[1,0.1,0.1]);
title('(c) Law of motion for N, $dN(B,N)$', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
ylim([Nmin Nmax])
view([60, 15])

print -dpdf g15_NPLM
print -depsc g15_NPLM
savefig(myfig,'g15_NPLM.fig');

%%
% phase diagram

myfig=figure(30);
set(myfig, 'Position', [0 0 800 800])

plot(PLM0(:,1),PLM0(:,2),'-','Color',[0.1,0.3,1],'linewidth',2);
hold on;
plot(NPLM0(:,1),NPLM0(:,2),'--r','linewidth',2);
title('Phase diagram', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
ylim([Nmin Nmax])

% upper left arrow
arrow_ = [0.16 0.82];
arrowH = [arrow_(1) arrow_(1)+0.05];
arrowV = [arrow_(2) arrow_(2)];
annotation('arrow',arrowH,arrowV);
arrowH = [arrow_(1) arrow_(1)];
arrowV = [arrow_(2) arrow_(2)-0.05];
annotation('arrow',arrowH,arrowV);

% upper left second arrow
arrow_ = [0.355 0.600];
arrowH = [arrow_(1) arrow_(1)-0.05];
arrowV = [arrow_(2) arrow_(2)];
annotation('arrow',arrowH,arrowV);
arrowH = [arrow_(1) arrow_(1)];
arrowV = [arrow_(2) arrow_(2)+0.05];
annotation('arrow',arrowH,arrowV);

% upper right arrow
arrow_ = [0.75 0.75];
arrowH = [arrow_(1) arrow_(1)-0.05];
arrowV = [arrow_(2) arrow_(2)];
annotation('arrow',arrowH,arrowV);
arrowH = [arrow_(1) arrow_(1)];
arrowV = [arrow_(2) arrow_(2)-0.05];
annotation('arrow',arrowH,arrowV);

% central arrow
arrow_ = [0.52 0.45];
arrowH = [arrow_(1) arrow_(1)+0.05];
arrowV = [arrow_(2) arrow_(2)];
annotation('arrow',arrowH,arrowV);
arrowH = [arrow_(1) arrow_(1)];
arrowV = [arrow_(2) arrow_(2)-0.05];
annotation('arrow',arrowH,arrowV);

% lower left arrow
arrow_ = [0.35 0.35];
arrowH = [arrow_(1) arrow_(1)+0.05];
arrowV = [arrow_(2) arrow_(2)];
annotation('arrow',arrowH,arrowV);
arrowH = [arrow_(1) arrow_(1)];
arrowV = [arrow_(2) arrow_(2)+0.05];
annotation('arrow',arrowH,arrowV);

% lower right arrow
arrow_ = [0.75 0.15];
arrowH = [arrow_(1) arrow_(1)-0.05];
arrowV = [arrow_(2) arrow_(2)];
annotation('arrow',arrowH,arrowV);
arrowH = [arrow_(1) arrow_(1)];
arrowV = [arrow_(2) arrow_(2)+0.05];
annotation('arrow',arrowH,arrowV);

plot(SSS_points(:,1),SSS_points(:,2),' Ok','linewidth',2)
plot(B_ss,N_ss,' s','Color',[0,0.5,0.1],'linewidth',5,'markers',10)
plot(B3,N3,' .k','linewidth',2)

text(1.05,2.70,'Low-leverage stable SSS', 'interpreter','latex','HorizontalAlignment','left','FontSize',12)
text(1.45,2.30,'Unstable SSS', 'interpreter','latex','HorizontalAlignment','left','FontSize',12)
text(2.00,1.78,'Baseline stable SSS', 'interpreter','latex','HorizontalAlignment','left','FontSize',12)
text(2.12,1.55,'Arbitrary high leverage point', 'interpreter','latex','HorizontalAlignment','left','FontSize',12)

legend({'$h(B,N)=0$','$dN(B,N)=0$','Stochastic steady state','Deterministic steady state'},'Location','southwest', 'interpreter','latex','FontSize',12)
grid

print -dpdf g30_phase
print -depsc g30_phase
savefig(myfig,'g30_phase.fig');



