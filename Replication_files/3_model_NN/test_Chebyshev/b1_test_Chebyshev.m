% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution


Y=[];
X1=[];
X2=[];

for it_sim=1:multi_sim
    Y = [Y ; (Bsim((it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim)-Bsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1))/dt];
    X1= [X1; Bsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1)];
    X2= [X2; Nsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1)];
end

order_cheb = 6; %order of the chebyshev approximation
Xn = f1_cheb2(X1,X2,Bmin,Bmax,Nmin,Nmax,order_cheb);

beta = (Xn'*Xn)^-1*Xn'*Y;
Y_fit = Xn*beta;
Y_RMSE_CH = (mean((Y-Y_fit).^2))^0.5;                            
Y_R_squared = 1- ( sum((Y-Y_fit).^2) ) / ( sum((Y-mean(Y)).^2) ); 


% now use LR to calculate the corresponding PLM - first on fine grid

X1mm=    BB_grid_2D ;
X2mm=    NN_grid_2D ;

X1m = reshape(X1mm,[nval_BB*nval_NN,1]);
X2m = reshape(X2mm,[nval_BB*nval_NN,1]);

Xm = f1_cheb2(X1m,X2m,Bmin,Bmax,Nmin,Nmax,order_cheb);
Y_mat2 = Xm*beta;
PLM_surf_CH = reshape(Y_mat2,size(X1mm));


% now use LR to calculate the corresponding PLM - now on HJB grid

X1mm=    squeeze(B(1,1,:,:)) ;
X2mm=    squeeze(N(1,1,:,:)) ;

X1m = reshape(X1mm,[nval_B*nval_N,1]);
X2m = reshape(X2mm,[nval_B*nval_N,1]);

Xm = f1_cheb2(X1m,X2m,Bmin,Bmax,Nmin,Nmax,order_cheb);
Y_mat2 = Xm*beta;
Y_mat2s= reshape(Y_mat2,size(X1mm));
for ia=1:nval_a
    for iz=1:nval_z
        PLM_CH(ia,iz,:,:)=Y_mat2s(:,:);
    end
end

