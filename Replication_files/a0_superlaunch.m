% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

cd 0_data
    a2_launch
cd ..

a2_launch_model             % This takes approximately  11 hours
a4_launch_sigma             % This takes approximately 188 hours (7.8 days)
a6_launch_z1                % This takes approximately 172 hours (7.2 days)
a8_launch_rhohat            % This takes approximately 112 hours (4.7 days)
a9_launch_robustness        % This takes approximately  24 hours



% TOTAL RUNNING TIME FOR ALL THIS:
% almost three weeks on an Intel i5-6500 (4 cores, 3.2GHz) with >12GB of RAM and the parallel toolbox
% twice that much if there's no parallel toolbox or not enough RAM to use it (with 8GB it struggles)


% And we're not launching the z1_sigma grid from here because it would take too long: another 630 hours (26.3 days)
% Better get a workstation for that, and launch several instances of Matlab in parallel

