% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

samplesize_quarters = size(Y_obs,1);          % number of quarters in the sample
samplesize_dt = samplesize_quarters/4/dt;     % number of dt-periods in the sample (if dt=1/12, number of months)
dt_per_q=1/dt/4;                              % if dt=1/12, this is months per quarter, i.e. 3
estim_likelihood=zeros(samplesize_quarters-1,1);

% we have to reformulate the solution of the model so that it is written in terms of Y instead of N

B20 = repmat(B,1,1,1,20); % we are going to use a 20-times finer grid for Y than we used for N (without this, approximation errors break everything)

% define Y grid

Ymin = Zeta*(Bmin+Nmin).^alpha;
Ymax = Zeta*(Bmax+Nmax).^alpha;

dY   = (Ymax-Ymin)/(nval_N*20-1);           % size of Y jump
Y_grid   = linspace(Ymin,Ymax,nval_N*20)';  % 1D - output
Y        = zeros(nval_a, nval_z, nval_B, nval_N*20);  % wouldn't need to be 4D but well, just so I can reuse old code below...
for ia=1:nval_a
    for iz=1:nval_z
        for iB=1:nval_B
            Y(ia,iz,iB,:)=Y_grid;
        end
    end
end

% initial values

Yt=Y_obs(1);

Yt_posD =floor((Yt-Ymin)/dY)+1;
Yt_posU = ceil((Yt-Ymin)/dY)+1;
wY=(Y_grid(Yt_posU)-Yt)/dY;     % Definition of the weights

f_matrix = zeros(nval_B,nval_N*20);
f_matrix(:,Yt_posD)=   wY ;
f_matrix(:,Yt_posU)=(1-wY);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constructing my very big sparse matrix A3 %
%       but now on the new space (B,Y)      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iN=1:nval_N*20
    A2_estim{iN}=sparse(nval_B,nval_B);
end

for iN=1:nval_N*20
    Xmat{iN}=sparse(nval_B,nval_B);
end

for iN=1:nval_N*20
    Pmat{iN}=sparse(nval_B,nval_B);
end

A3_estim=sparse(nval_B*nval_N*20,nval_B*nval_N*20);

% use NN1 to calculate the corresponding PLM - now using Y_grid

X1mm=    squeeze(B20(1,1,:,:)) ;
X2mm=    squeeze(Y(1,1,:,:)) ;

X2mm=X2mm.^(1/alpha)-X1mm;

X1m = reshape(X1mm,[nval_B*nval_N*20,1]);
X2m = reshape(X2mm,[nval_B*nval_N*20,1]);

X1mr=(X1m-X1r3)/X1r4;
X2mr=(X2m-X2r3)/X2r4;

Xmr = [X1mr X2mr];

Yr_mat2 = f2_NN_eval(Xmr,network_width,NN1_flat);
Y_mat2 = (Yr_mat2*Yr4)+Yr3;
Y_mat2s= reshape(Y_mat2,size(X1mm));
elem_h_Y=zeros(size(B20));
for ia=1:nval_a
    for iz=1:nval_z
        elem_h_Y(ia,iz,:,:)=Y_mat2s(:,:);
    end
end

% prepare the new elem_m and elem_s

elem_m_Y = alpha * (Y.^((alpha-1)/alpha)).*( elem_h_Y + alpha*Y + ((alpha-1)*sigma2/2-delta)*(Y.^(1/alpha)) - ( alpha*(Y.^((alpha-1)/alpha)) - delta - sigma2*(Y.^(1/alpha))./(Y.^(1/alpha)-B20) ).*B20 ) - (alpha*(Y.^((alpha-1)/alpha))*rhohat).*(Y.^(1/alpha)-B20);
elem_s_Y = alpha * sigma * Y;

elem_b_Y = - elem_h_Y/dB - elem_m_Y/dY - ((elem_s_Y/dY).^2);

% redefine elem_x and elem_r so that they refer to the new space (B,Y)

elem_r_Y       = ((alpha*sigma*Y/dY).^2)/2;
elem_x_Y       = elem_m_Y/dY + elem_r_Y;

% Fill A2 collection

jumper = 1;
for iN=1:nval_N*20
    for iB=1:nval_B
        A2_estim{iN}((iB-1)*jumper+1:iB*jumper,(iB-1)*jumper+1:iB*jumper) = elem_b_Y(1,1,iB,iN);
    end
end
for iN=1:nval_N*20
    for iB=1:nval_B-1
        A2_estim{iN}((iB-1)*jumper+1:iB*jumper,iB*jumper+1:(iB+1)*jumper) = (elem_h_Y(1,1,iB,iN)/dB);
    end
    A2_estim{iN}((nval_B-1)*jumper+1:nval_B*jumper,(nval_B-1)*jumper+1:nval_B*jumper)=A2_estim{iN}((nval_B-1)*jumper+1:nval_B*jumper,(nval_B-1)*jumper+1:nval_B*jumper)+(elem_h_Y(1,1,nval_B,iN)/dB);
end

% Prepare X2 and P2 collections

for iN=1:nval_N*20
    for iB=1:nval_B
        Xmat{iN}((iB-1)*jumper+1:iB*jumper,(iB-1)*jumper+1:iB*jumper) = elem_x_Y(1,1,iB,iN) ;
    end
end

for iN=1:nval_N*20
    for iB=1:nval_B
        Pmat{iN}((iB-1)*jumper+1:iB*jumper,(iB-1)*jumper+1:iB*jumper) = elem_r_Y(1,1,iB,iN) ;
    end
end

% Fill A3

jumper = nval_B;
for iN=1:nval_N*20
    A3_estim((iN-1)*jumper+1:iN*jumper,(iN-1)*jumper+1:iN*jumper) = A2_estim{iN};
end
for iN=1:nval_N*20-1
    A3_estim(iN*jumper+1:(iN+1)*jumper,(iN-1)*jumper+1:iN*jumper) = Pmat{iN+1};
end
for iN=1:nval_N*20-1
    A3_estim((iN-1)*jumper+1:iN*jumper,iN*jumper+1:(iN+1)*jumper) = Xmat{iN};
end
A3_estim(1:jumper,1:jumper) = A3_estim(1:jumper,1:jumper) + Pmat{1};
A3_estim((nval_N*20-1)*jumper+1:nval_N*20*jumper,(nval_N*20-1)*jumper+1:nval_N*20*jumper) = A3_estim((nval_N*20-1)*jumper+1:nval_N*20*jumper,(nval_N*20-1)*jumper+1:nval_N*20*jumper) + Xmat{nval_N*20};
    

% Calculate likelihood for the observation of each quarter

for it1=1:samplesize_quarters-1

    % probability propagation

    f_vector0=f_matrix(:);
    
    f_vector1=(speye(nval_B*nval_N*20)-A3_estim'*dt) \f_vector0;
    f_vector1=f_vector1/(sum(f_vector1)*dY);
    f_vector2=(speye(nval_B*nval_N*20)-A3_estim'*dt) \f_vector1;
    f_vector2=f_vector2/(sum(f_vector2)*dY);
    f_vector3=(speye(nval_B*nval_N*20)-A3_estim'*dt) \f_vector2;
    f_vector3=f_vector3/(sum(f_vector3)*dY);

    f_vector=(f_vector1+f_vector2+f_vector3)/3;
    
    f_matrix = reshape(f_vector,size(f_matrix));

    % now for the observed values

    Yt=Y_obs(1+it1);

    Yt_posD =floor((Yt-Ymin)/dY)+1;
    Yt_posU = ceil((Yt-Ymin)/dY)+1;
    wY=(Y_grid(Yt_posU)-Yt)/dY;     % Definition of the weights

    f_obs_matrix = zeros(nval_B,nval_N*20);

    f_obs_matrix(:,Yt_posD)=   wY *f_matrix(:,Yt_posD);
    f_obs_matrix(:,Yt_posU)=(1-wY)*f_matrix(:,Yt_posU);
    
    estim_likelihood(it1)=sum(sum(f_obs_matrix))*dY;

    if estim_likelihood(it1)<0
    	estim_likelihood(it1)=min(estim_likelihood(1:it1-1));    % sometimes estim_likelihood comes out negative; when that happens, replace it and show error message, but continue calculations
        disp(['ERROR: negative likelihood found at t=' num2str(it1)])
    end

    f_matrix = f_obs_matrix / (sum(sum(f_obs_matrix))*dY);

end

estim_likelihood_log=log(estim_likelihood);
estim_likelihood_logsum=sum(estim_likelihood_log);

disp('loglikelihood (estim_likelihood_logsum):');
disp(estim_likelihood_logsum);

