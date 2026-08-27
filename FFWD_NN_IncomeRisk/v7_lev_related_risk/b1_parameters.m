% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2023)
% Financial Frictions and the Wealth Distribution
% New feature, leverage-related labor-income risk




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



alpha         = 0.35;                       % Capital share   
Zeta          = 1;                          % Aggregate productivity
delta         = 0.1;                        % Depreciation
gamma         = 2;                          % CRRA utility with parameter s
rho           = 0.05;                       % discount rate
rhohat        = 0.04971;                    % discount rate for banks
la1           = 0.986;                      % transition probability
la2           = 0.052;
la            = [la1,la2];                  % vector of transition prob

amin          = 0.5;                          % borrowing constraint changed, previous 0.
amax          = 20;                         % max value of individual savings
z1            = 0.72;                       % labor productivity
z2            = 1 + la2/la1 * (1-z1);

% Parameters governing leverage-related labor-income risk
mean_logz = la2/(la1+la2)* log(z1) + la1/(la1+la2)*log(z2);
var_logz = la2/(la1+la2)*(log(z1)- mean_logz)^2 + la1/(la1+la2)*(log(z2)- mean_logz)^2;
V_ss = var_logz;                            % baseline variance of log productivity
eta_l = 0.01;                               % local response of Var(log labor income) to log leverage at leverage_ref
h_l = 0.20;                                 % smoothness of the leverage-risk transition
leverage_ref = 0.75;                        % leverage level at which theta(B,N) equals one

Bmin          = 0.7;                        % relevant range for aggregate savings
Bmax          = 2.7;
Nmin          = 1.2;                        % relevant range for aggregate equity
Nmax          = 3.2;
                       
nval_a        = 501;                        % number of points in amin-to-amax range (individual savings)
nval_z        = 2;                          % number of options for z (the idiosincratic shock)
nval_B        = 8;                          % number of points in Bmin-to-Bmax range (aggregate savings), on the coarse grid used for the HJB
nval_N        = 51;                         % number of points in Nmin-to-Nmax range (aggregate equity) , on the coarse grid used for the HJB

nval_BB       = 101;                        % finer grid, used for training the NN, for determining visited range and for the convergence criteria
nval_NN       = 101;                        % finer grid, used for training the NN, for determining visited range and for the convergence criteria

dt            = 1/12;                       % size of t jump
da  = (amax-amin)/(nval_a-1);               % size of a jump
dB  = (Bmax-Bmin)/(nval_B-1);               % size of B jump on the coarse grid
dN  = (Nmax-Nmin)/(nval_N-1);               % size of N jump on the coarse grid

dBB = (Bmax-Bmin)/(nval_BB-1);              % size of B jump on the fine grid
dNN = (Nmax-Nmin)/(nval_NN-1);              % size of B jump on the fine grid

sigma         = 0.0140;                    % original sigma
sigma2        = sigma^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters for the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

multi_sim    =   4    ;                     % number of simulation starts (i.e. one plus number of times we stop the simulation and start again from the ss)
delay_sim    =  500/dt;                     % number of initial periods in each simulation that will not be used (pre-heat)
used_sim     = 5000/dt;                     % number of periods in each simulation that will actually be used
each_sim     = delay_sim+used_sim;          % number of periods on each simulation
nval_sim     = multi_sim*each_sim;          % total number of periods in all simulations
rngseed1     = 123;                         % RNG seed for calculating the shocks for the simulation
                                            % there used to be a second seed for the mini-batch selection in the NN training, but we got rid of that by moving to batch gradient descent

maxitHJB = 100;                             % max number of iterations for the HJB
critHJB  = 10^(-6);                         % convergence crit for the HJB 
weHJB    = 0.5;                             % relaxation algorithm for HJB
Delta    = 1000;                            % Delta in HJB
maxitPLM = 200;                             % max number of iterations of the full algorithm
critPLM  = 0.00050;                         % convergence crit for determining that the PLM has converged
wePLM    = 0.3;                             % Initial weigth in the relaxation algorithm for PLM convergence
wePLM1   = 0.9;                             % reduction of the relaxation algorithm: wePLM = wePLM*wePLM1+wePLM2
wePLM2   = 0.005;                           % reduction of the relaxation algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural network parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
network_width  =    16;                     % Number of neurons in the hidden layer
mynoise        =     1;                     % Size of random initial NN parameters
lambda         =   0.1;                     % NN regularization parameter
NN_iters       = 10000;                     % Number of iterations to train the network
NN_starts      =    10;                     % Number of random restarts of the network training (only affects the first step, then it just starts from the previous NN)
learning_speed =  0.01;                     % Only affects the first step, then this becomes adaptative
reglimY        =     4;                     % NN input normalization
reglimX        =     4;                     % NN input normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

% Construct the leverage-related labor-income grid. The probability-
% weighted normalization keeps average labor income equal to w while
% theta_BN changes only its cross-sectional dispersion.
avg_zi = zeros(nval_a, nval_z, nval_B, nval_N);  % grid for cross-sectional average of effective productivity
leverage_grid = B./N;
theta_BN = exp((eta_l*h_l/(2*V_ss)).* ...
    tanh(log(leverage_grid./leverage_ref)./h_l));

for iB=1:nval_B
    for iN=1:nval_N
        leverage_temp = B_grid(iB)/N_grid(iN);
        theta_temp = exp((eta_l*h_l/(2*V_ss))* ...
            tanh(log(leverage_temp/leverage_ref)/h_l));
        avg_zi(:,:,iB,iN) = la2/(la1+la2)*z1^theta_temp + ...
                            la1/(la1+la2)*z2^theta_temp;
    end
end

yi = w.*z.^theta_BN./avg_zi;
yi_previous = w.*z;

disp("Difference between previous income grid and current income grid")
temp = (yi_previous-yi)./yi;
summary(temp(:));

% Report idiosyncratic labor-income risk and its covariance with leverage.
income_state_prob = reshape([la2,la1]./(la1+la2),[1,nval_z,1,1]);
log_yi_report = log(yi(1,:,:,:));
mean_log_yi_report = sum(income_state_prob.*log_yi_report,2);
sd_log_yi_grid = squeeze(sqrt(sum(income_state_prob.* ...
    (log_yi_report-mean_log_yi_report).^2,2)));
leverage_report_grid = squeeze(leverage_grid(1,1,:,:));
leverage_ss = B_ss/N_ss;
mean_sd_log_yi = mean(sd_log_yi_grid(:));
sd_income_risk_grid = std(sd_log_yi_grid(:),1);
sd_leverage_report_grid = std(leverage_report_grid(:),1);
if sd_income_risk_grid == 0 || sd_leverage_report_grid == 0
    corr_sd_log_yi_leverage = NaN;
    warning('IncomeRiskReport:ZeroVariation', ...
        'Income risk or leverage has zero grid variation; correlation is undefined.');
else
    corr_sd_log_yi_leverage = mean((sd_log_yi_grid(:)-mean_sd_log_yi).* ...
        (leverage_report_grid(:)-mean(leverage_report_grid(:))))/ ...
        (sd_income_risk_grid*sd_leverage_report_grid);
end

fprintf('Mean conditional SD of log idiosyncratic income: %.8f\n',mean_sd_log_yi);
fprintf('Grid correlation of income SD with leverage: %.8f\n', ...
    corr_sd_log_yi_leverage);

% Plot cross-sectional labor-income risk as a function of leverage.
[leverage_plot,sort_index] = sort(leverage_report_grid(:));
sd_log_yi_plot = sd_log_yi_grid(sort_index);

figure('Name','Labor-income risk vs leverage');
plot(leverage_plot,sd_log_yi_plot,'LineWidth',1.5)
hold on
xline(leverage_ref,'k--','$B/N=0.75$', ...
    'Interpreter','latex','LabelVerticalAlignment','bottom')
xline(leverage_ss,'r:','$B_{ss}/N_{ss}$', ...
    'Interpreter','latex','LabelVerticalAlignment','top')
xlabel('$B/N$','Interpreter','latex')
ylabel('$\mathrm{SD}_z(\log y_i\mid B/N)$','Interpreter','latex')
title('Cross-sectional labor-income risk and leverage', ...
    'Interpreter','latex')
grid on
box on
hold off

% Calculate the mu^N = 0 locus and the second zero locus on the fine B grid.
N_muN_zero = NaN(numel(BB_grid),numel(NN_grid));
N_second_zero = NaN(numel(BB_grid),numel(NN_grid));

for iB = 1:numel(BB_grid)
    B_temp = BB_grid(iB);

    muN_equation = @(N_temp) (...
        alpha.*(N_temp+B_temp).^(alpha-1)-delta-rhohat + ...
        sigma^2.*((N_temp+B_temp)./N_temp-1).* ...
        (N_temp+B_temp)./N_temp).*N_temp;
    muN_on_N_grid = muN_equation(NN_grid);
    N_roots = NN_grid(muN_on_N_grid == 0);
    root_intervals = find(muN_on_N_grid(1:end-1).* ...
        muN_on_N_grid(2:end) < 0);

    for iroot = 1:numel(root_intervals)
        root_interval = root_intervals(iroot);
        N_roots(end+1,1) = fzero(muN_equation, ...
            [NN_grid(root_interval),NN_grid(root_interval+1)]); %#ok<SAGROW>
    end

    if ~isempty(N_roots)
        N_roots = unique(N_roots);
        N_muN_zero(iB,1:numel(N_roots)) = N_roots;
    end

    second_equation = @(N_temp) ...
        (1-alpha).*(N_temp+B_temp).^alpha + ...
        (alpha.*(N_temp+B_temp).^(alpha-1)-delta- ...
        sigma^2.*(N_temp+B_temp)./N_temp).*B_temp - 1.1;
    second_on_N_grid = second_equation(NN_grid);
    N_roots = NN_grid(second_on_N_grid == 0);
    root_intervals = find(second_on_N_grid(1:end-1).* ...
        second_on_N_grid(2:end) < 0);

    for iroot = 1:numel(root_intervals)
        root_interval = root_intervals(iroot);
        N_roots(end+1,1) = fzero(second_equation, ...
            [NN_grid(root_interval),NN_grid(root_interval+1)]); %#ok<SAGROW>
    end

    if ~isempty(N_roots)
        N_roots = unique(N_roots);
        N_second_zero(iB,1:numel(N_roots)) = N_roots;
    end
end

if all(isnan(N_muN_zero),'all')
    warning('ZeroCurve:NoMuNRoot', ...
        'No mu^N = 0 root was detected between Nmin and Nmax.');
end
if all(isnan(N_second_zero),'all')
    warning('ZeroCurve:NoSecondRoot', ...
        'No root of the second equation was detected between Nmin and Nmax.');
end

figure('Name','Aggregate zero loci');
hold on
muN_line = plot(BB_grid,N_muN_zero,'b-','LineWidth',1.5);
second_line = plot(BB_grid,N_second_zero,'r--','LineWidth',1.5);
xlabel('$B_t$','Interpreter','latex')
ylabel('$N_t$','Interpreter','latex')
title('Aggregate zero loci','Interpreter','latex')
legend([muN_line(1),second_line(1)], ...
    {'$\mu^N=0$', ...
    '$(1-\alpha)(N+B)^\alpha+[\alpha(N+B)^{\alpha-1}-\delta-\sigma^2(N+B)/N]B=0$'}, ...
    'Interpreter','latex','Location','best')
grid on
close;
