% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

% This solves some additional versions of the model, as a robustness check


cd 9_robustness_BNsupport  % this changes slightly the bounds of the (B,N) support
    a2_launch
cd ..


cd 9_robustness_IP         % this uses a natural interpolator in the visited area and NN extrapolation in the rest
    a2_launch
cd ..


cd 9_robustness_IPlin      % this uses a natural interpolator in the visited area and linear extrapolation in the rest
    a2_launch
cd ..


cd 9_robustness_rng        % this changes the seed for the random number generator used for the simulation
    a2_launch
cd ..

close all

