% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

myfig_=figure(1);
set(myfig_, 'Position', [0 0 800 600])

it_rhohat_=0;
mytitle=['(a) $\hat{\rho}=0.04941$';'(b) $\hat{\rho}=0.04951$';'(c) $\hat{\rho}=0.04961$';'(d) $\hat{\rho}=0.04971$';'(e) $\hat{\rho}=0.04981$';'(f) $\hat{\rho}=0.04991$'];

for it_rhohat = [4951 4961 4971 4981 4991 4996]
    
    it_rhohat
    it_rhohat_=it_rhohat_+1;
    
    if it_rhohat ~= 4971
        eval(['load ''./rhohat_' num2str(it_rhohat) '/z_FinalWorkspace.mat'';']);
    else
        load '../3_model_NN/z_FinalWorkspace.mat'
    end
    
    cd phase
        b1_phase
    cd ..
    
    cd rhohat_DSS
        a2_launch
    cd ..

    subplot(2,3,it_rhohat_);
    plot(PLM0(:,1),PLM0(:,2),'-','Color',[0.1,0.3,1],'linewidth',2);
    hold on;
    plot(NPLM0(:,1),NPLM0(:,2),'--r','linewidth',2);
    plot(SSS_points(:,1),SSS_points(:,2),' Ok','linewidth',2)
    plot(B_ss,N_ss,' s','Color',[0,0.5,0.1],'linewidth',1,'markers',6)
    title(mytitle(it_rhohat_,:), 'interpreter','latex','FontSize',12);
    xlabel('debt ($B$)', 'interpreter','latex','FontSize',12);
    ylabel('equity ($N$)', 'interpreter','latex','FontSize',12);
    xlim([Bmin Bmax])
    ylim([Nmin Nmax])
    grid
    
end

legend({'$h(B,N)=0$','$\mu^N(B,N)=0$','SSS','DSS'},'Location','best', 'interpreter','latex','FontSize',10)

%%

print -dpdf h70_phase6
savefig(myfig_,'h70_phase6.fig');

print -dpdf g70_phase6
