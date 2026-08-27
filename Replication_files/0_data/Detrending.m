
% Detrends the data using a bandpass filter
Y = xlsread('GDPC1','FRED Graph','B160:B295');

fL = 1/60;
fH = 1/20;
K_BP = 12/1;
logY =      bandpass(fL,fH,K_BP,log(Y));
Y_obs = exp(logY);

% hold on
% load z_obs
% plot(Y_obs)
save '../z_obs_data_bandpass.mat' Y_obs