% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

sigma_grid =nan(21,1);
sss_grid=nan(21,3);

it_sigma_=0;

for it_sigma = [120 130 140 150 160 170 180 190 200 210 220 230 240 260 280 300 320 340 360 380 400]

    it_sigma
    it_sigma_=it_sigma_+1;
    
    if it_sigma ~= 140
        eval(['load ''./sigma_' num2str(it_sigma) '/z_FinalWorkspace.mat'';']);
        cd 'phase';
        B3=0;
        N3=0;
        b1_phase
        close all
        cd ..
    else
        load '../3_model_NN/z_FinalWorkspace.mat'
        cd 'phase';
        B3=0;
        N3=0;
        b1_phase
        close all
        cd ..
    end
    
    sigma_grid(it_sigma_)=it_sigma/10000;
    
    if size(SSS_points,1)==1
        sss_grid(it_sigma_,1)=NaN;
        sss_grid(it_sigma_,2)=NaN;
        sss_grid(it_sigma_,3)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
    elseif size(SSS_points,1)==2
        sss_grid(it_sigma_,1)=NaN;
        sss_grid(it_sigma_,2)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
        sss_grid(it_sigma_,3)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
    else
        sss_grid(it_sigma_,1)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
        sss_grid(it_sigma_,2)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
        sss_grid(it_sigma_,3)=(SSS_points(3,1)+SSS_points(3,2))/SSS_points(3,2);
    end
end

%%

myfig=figure(1);
set(myfig, 'Position', [0 0 600 600])
plot(sigma_grid,sss_grid(:,1),' *','Color',[1,0.1,0.1],'linewidth',2);
hold on
plot(sigma_grid,sss_grid(:,2),' *','Color',[0.5,0.5,0.5],'linewidth',2);
plot(sigma_grid,sss_grid(:,3),' *','Color',[0,0.5,0.1],'linewidth',2);
plot(sigma_grid,ones(size(sss_grid(:,3)))*(B_ss+N_ss)/N_ss,'-','Color',[0.1,0.3,1],'linewidth',2);
plot([0.0140,0.0140],[-1000,+1000],'k--','linewidth',1);
title('Stochastic steady states for different values of $\sigma$', 'interpreter','latex','FontSize',14);
xlabel('$\sigma$','interpreter','latex','FontSize',12)
ylabel('leverage ($K/N$)','interpreter','latex','FontSize',12)
legend({'LL-SSS','Unstable SSS','HL-SSS','DSS'},'Location','southeast', 'interpreter','latex','FontSize',10)
ylim([0,3])
ax = gca;
ax.YAxis.Exponent = 0;
grid

print -dpdf h72_SSS_sigma
savefig(myfig,'h72_SSS_sigma.fig');

title(' ', 'interpreter','latex','FontSize',14);
print -dpdf g72_SSS_sigma

