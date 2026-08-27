% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

% This solves the model for different values of sigma, and then evaluates the likelihood of each one

cd 4_sigma
    for it_sigma=130:2:150
        if it_sigma ~= 140
            eval(['cd sigma_' num2str(it_sigma)])
            a2_launch
            cd ..
        end
    end
cd ..


cd 5_estimate_sigma_obsdata
    a2_launch
cd ..


cd 4_sigma
    for it_sigma = [120 160 170 180 190 200 210 220 230 240 260 280 300 320 340 360 380 400]
        if it_sigma ~= 140
            eval(['cd sigma_' num2str(it_sigma)])
            a2_launch
            cd ..
        end
    end
cd ..


cd 4_sigma
    b1_plot_phase6
    b2_plot_SSS
cd ..





