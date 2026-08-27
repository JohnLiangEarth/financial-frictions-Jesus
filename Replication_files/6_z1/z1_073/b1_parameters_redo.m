% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dt            = 1/12;                       % size of t jump
da  = (amax-amin)/(nval_a-1);               % size of a jump
dB  = (Bmax-Bmin)/(nval_B-1);               % size of B jump on the coarse grid
dN  = (Nmax-Nmin)/(nval_N-1);               % size of N jump on the coarse grid

dBB = (Bmax-Bmin)/(nval_BB-1);              % size of B jump on the fine grid
dNN = (Nmax-Nmin)/(nval_NN-1);              % size of B jump on the fine grid


% VARIABLES:
% the main ones will be 4-dimmensional matrices
% of size (nval_a, nval_z, nval_B, nval_N), with subscripts (ia, iz, iB, iN)
% there will be lots of repeated information in them (many of these matrices
% will be flat in all dimensions except one), but I can afford it
% (in terms of memory these matrices are not that big anyway)

a_grid   = linspace(amin,amax,nval_a)';           % 1D - assets
z_grid   = [z1,z2];                               % 1D - productivity
B_grid   = linspace(Bmin,Bmax,nval_B)';           % 1D - capital
N_grid   = linspace(Nmin,Nmax,nval_N)';           % 1D - TFP

BB_grid  = linspace(Bmin,Bmax,nval_BB)';          % finer grid, used only for determining visited range and convergence criteria
NN_grid  = linspace(Nmin,Nmax,nval_NN)';          % finer grid, used only for determining visited range and convergence criteria

a      = zeros(nval_a, nval_z, nval_B, nval_N);
z      = zeros(nval_a, nval_z, nval_B, nval_N);
B      = zeros(nval_a, nval_z, nval_B, nval_N);
N      = zeros(nval_a, nval_z, nval_B, nval_N);

BB_grid_2D = zeros(nval_BB, nval_NN);
NN_grid_2D = zeros(nval_BB, nval_NN);

for iz=1:nval_z    % repmat would be faster, but this is clearer
    for iB=1:nval_B
        for iN=1:nval_N
            a(:,iz,iB,iN)=a_grid;
        end
    end
end

for ia=1:nval_a
    for iB=1:nval_B
        for iN=1:nval_N
            z(ia,:,iB,iN)=z_grid;
        end
    end
end

for ia=1:nval_a
    for iz=1:nval_z
        for iN=1:nval_N
            B(ia,iz,:,iN)=B_grid;
        end
    end
end

for ia=1:nval_a
    for iz=1:nval_z
        for iB=1:nval_B
            N(ia,iz,iB,:)=N_grid;
        end
    end
end

for iN=1:nval_NN
    BB_grid_2D(:,iN)=BB_grid;
end
for iB=1:nval_BB
    NN_grid_2D(iB,:)=NN_grid;
end

a2=squeeze(a(:,:,1,1));    % this one is 2D instead of 4D, we need it for a simpler KFE algorithm

% Interest rates and wages (4D matrices that don't depend on anything but parameters) - WE ARE ASSUMING L=1
r =  alpha * Zeta * ((B+N).^(alpha-1)) - delta - sigma2*((B+N)./N);
w = (1-alpha) * Zeta * (B+N).^alpha;
