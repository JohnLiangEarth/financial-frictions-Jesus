% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

rhohat_grid =nan(12,1);
sss_rhohat_grid=nan(12,3);
ss_rhohat_grid=nan(12,1);

it_rhohat_=0;

rhohat_set = [4901 4911 4921 4931 4941 4951 4961 4971 4981 4986 4991 4996]

for it_rhohat = rhohat_set

    it_rhohat
    it_rhohat_=it_rhohat_+1;
    
    if it_rhohat ~= 4971
        eval(['load ''./rhohat_' num2str(it_rhohat) '/z_FinalWorkspace.mat'';']);
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
    
    rhohat_grid(it_rhohat_)=it_rhohat/100000;
    
    if size(SSS_points,1)==1
        sss_rhohat_grid(it_rhohat_,1)=NaN;
        sss_rhohat_grid(it_rhohat_,2)=NaN;
        sss_rhohat_grid(it_rhohat_,3)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
    elseif size(SSS_points,1)==2
        sss_rhohat_grid(it_rhohat_,1)=NaN;
        sss_rhohat_grid(it_rhohat_,2)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
        sss_rhohat_grid(it_rhohat_,3)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
    else
        sss_rhohat_grid(it_rhohat_,1)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
        sss_rhohat_grid(it_rhohat_,2)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
        sss_rhohat_grid(it_rhohat_,3)=(SSS_points(3,1)+SSS_points(3,2))/SSS_points(3,2);
    end
    
    cd rhohat_DSS
    a2_launch
    close all
    cd ..
    ss_rhohat_grid(it_rhohat_,1)=(B_ss+N_ss)/N_ss;
    
end

%%

myfig=figure(78);
set(myfig, 'Position', [0 0 600 600])
plot(rhohat_grid,sss_rhohat_grid(:,1),' *','Color',[1,0.1,0.1],'linewidth',2);
hold on
plot(rhohat_grid,sss_rhohat_grid(:,2),' *','Color',[0.5,0.5,0.5],'linewidth',2);
plot(rhohat_grid,sss_rhohat_grid(:,3),' *','Color',[0,0.5,0.1],'linewidth',2);
plot(rhohat_grid(1:11),ss_rhohat_grid(1:11,1),'-','Color',[0.1,0.3,1],'linewidth',2);
plot([0.72,0.72],[-1000,+1000],'k--','linewidth',1);
title('Stochastic steady states for different values of $z_1$', 'interpreter','latex','FontSize',14);
xlabel('$\hat{\rho}$','interpreter','latex','FontSize',12)
ylabel('leverage ($K/N$)','interpreter','latex','FontSize',12)
legend({'LL-SSS','Unstable SSS','HL-SSS','DSS'},'Location','southeast', 'interpreter','latex','FontSize',10)
xlim([0.048 0.05])
ylim([0,3])
ax = gca;
ax.YAxis.Exponent = 0;
grid

print -dpdf h78_SSS_rhohat
savefig(myfig,'h78_SSS_rhohat.fig');

title(' ', 'interpreter','latex','FontSize',14);
print -dpdf g78_SSS_rhohat

