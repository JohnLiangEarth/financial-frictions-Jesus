% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

clear NN_all loss_final;

% prepare our observations (using only the observations after the initial pre-heat period of each simulation)

Y=[];
X1=[];
X2=[];

for it_sim=1:multi_sim
    Y = [Y  ; (Bsim((it_sim-1)*each_sim+delay_sim+1:it_sim*each_sim)-Bsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1))/dt];
    X1= [X1 ; Bsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1)];
    X2= [X2 ; Nsim((it_sim-1)*each_sim+delay_sim:it_sim*each_sim-1)];
end

X = [X1 X2];

% grid approximation, with LR on each grid knot

PLM_IP_finegrid = zeros(nval_BB,nval_NN)*NaN;

for iB=1:nval_BB
    for iN=1:nval_NN
        Y_=[];
        X1_=[];
        X2_=[];
        for it_sim=1:size(Y,1)
            if abs(X1(it_sim)-BB_grid(iB))<dBB/2
                if abs(X2(it_sim)-NN_grid(iN))<dNN/2
                    Y_=[Y_;Y(it_sim)];
                    X1_=[X1_;X1(it_sim)];
                    X2_=[X2_;X2(it_sim)];
                end
            end
        end
        if size(X1_,1)>5
            X0_=ones(size(X1_));
            X = [X0_ X1_ X2_];
            beta = (X'*X)^-1*X'*Y_;
            PLM_IP_finegrid(iB,iN) = [1 BB_grid(iB) NN_grid(iN)]*beta;
        end
    end
end

% we will use the knots in that grid of LR approximations to train the neural network
% (this helps because feeding less noisy data makes the NN training work faster, and takes random noise away from the full algorithm)

Y_=[];
X1_=[];
X2_=[];
for iB=1:nval_BB
    for iN=1:nval_NN
        if not(isnan(PLM_IP_finegrid(iB,iN)))
            Y_=[Y_;PLM_IP_finegrid(iB,iN)];
            X1_=[X1_;BB_grid_2D(iB,iN)];
            X2_=[X2_;NN_grid_2D(iB,iN)];
        end
    end
end

% an alternative to the neural network would be this built-in interpolant (IP)
% it's very fast but linear extrapolation on two dimmensions generates ugly ridges
F = scatteredInterpolant(X1_,X2_,Y_,'natural','linear');

Y_fit_IP = F(X1,X2);
Y_error_IP = Y-Y_fit_IP;
Y_RMSE_IP = (mean(Y_error_IP.^2))^0.5;

Y_fit=Y_fit_IP;
Y_error=Y_error_IP;

% now use F to calculate the corresponding PLM - first on fine grid

X1mm=    BB_grid_2D ;
X2mm=    NN_grid_2D ;

X1m = reshape(X1mm,[nval_BB*nval_NN,1]);
X2m = reshape(X2mm,[nval_BB*nval_NN,1]);

PLM_finegrid_2 = F(X1m,X2m);
PLM_finegrid_2 = reshape(PLM_finegrid_2,size(X1mm));


% now use F to calculate the corresponding PLM - now on HJB grid

X1mm=    squeeze(B(1,1,:,:)) ;
X2mm=    squeeze(N(1,1,:,:)) ;

X1m = reshape(X1mm,[nval_B*nval_N,1]);
X2m = reshape(X2mm,[nval_B*nval_N,1]);

Y_mat2s = F(X1m,X2m);
Y_mat2s = reshape(Y_mat2s,size(X1mm));

for ia=1:nval_a
    for iz=1:nval_z
        PLM_2(ia,iz,:,:)=Y_mat2s(:,:);
    end
end


PLM_IP = F(X1m,X2m);
PLM_IP = reshape(PLM_IP,size(X1mm));


Y_RMSE_NN2 = Y_RMSE_IP;
