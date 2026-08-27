% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

copyfile '../../../3_model_NN/f2_NN_eval.m';


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

X1mr=(X1m-X1r3)/X1r4;
X2mr=(X2m-X2r3)/X2r4;

Xmr = [X1mr X2mr];

Yr_mat2 = f2_NN_eval(Xmr,network_width,NN1_flat);
Y_mat2 = (Yr_mat2*Yr4)+Yr3;
PLM_surf= reshape(Y_mat2,[nval_BB nval_NN]);

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

PLM0=NaN(nval_BB,2);

for it1=1:nval_BB
    found=0;
    it3 = NaN;
    for it2=1:nval_NN
        if PLM_surf(it1,it2)>0
            it3=it2;
            found=1;
        end
    end
    if it3>=nval_NN
        found = 0;
    end
    if found == 1
        myweight=PLM_surf(it1,it3+1)/(PLM_surf(it1,it3+1)-PLM_surf(it1,it3));
        PLM0(it1,1)=BB_grid(it1);
        PLM0(it1,2)=NN_grid(it3)*myweight+NN_grid(it3+1)*(1-myweight);
    end
end


%%%%%%%%%%%%%%
% find NPLM0 %
%%%%%%%%%%%%%%

NPLM0=NaN(nval_BB,2);

for it1=1:nval_BB
    found=0;
    it3 = NaN;
    for it2=1:nval_NN
        if NPLM_surf(it1,it2)>0
            it3=it2;
            found=1;
        end
    end
    if it3>=nval_NN
        found = 0;
    end
    if found == 1
        myweight=NPLM_surf(it1,it3+1)/(NPLM_surf(it1,it3+1)-NPLM_surf(it1,it3));
        NPLM0(it1,1)=BB_grid(it1);
        NPLM0(it1,2)=NN_grid(it3)*myweight+NN_grid(it3+1)*(1-myweight);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find all Stochastic Steady States %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find SSS points (it's just for graphs, so: using linear approximation around zero of PLM and NPLM surfaces)
% Later I'll use a very long no-shocks simulation to refine these points, before using them to run IRFs etc
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



