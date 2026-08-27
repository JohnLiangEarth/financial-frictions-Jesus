% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

z1_grid =nan(22,1);
sss_z1_grid=nan(22,3);
ss_z1_grid=nan(22,1);

it_z1_=0;

z1_set = [67 68 70 72 73 75 77 78 80 82 83 85 87 88 90 91 92 93 94 95 96 97]

for it_z1 = z1_set

    it_z1
    it_z1_=it_z1_+1;
    
    if it_z1 ~= 72
        eval(['load ''./z1_0' num2str(it_z1) '/z_FinalWorkspace.mat'';']);
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
    
    z1_grid(it_z1_)=it_z1/100;
    
    if size(SSS_points,1)==1
        if it_z1 < 80
            sss_z1_grid(it_z1_,1)=NaN;
            sss_z1_grid(it_z1_,2)=NaN;
            sss_z1_grid(it_z1_,3)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
        else
            sss_z1_grid(it_z1_,1)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
            sss_z1_grid(it_z1_,2)=NaN;
            sss_z1_grid(it_z1_,3)=NaN;
        end
    elseif size(SSS_points,1)==2
        if it_z1 < 90
            sss_z1_grid(it_z1_,1)=NaN;
            sss_z1_grid(it_z1_,2)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
            sss_z1_grid(it_z1_,3)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
        else
            sss_z1_grid(it_z1_,1)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
            sss_z1_grid(it_z1_,2)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
            sss_z1_grid(it_z1_,3)=NaN;
        end
    else
        sss_z1_grid(it_z1_,1)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
        sss_z1_grid(it_z1_,2)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
        sss_z1_grid(it_z1_,3)=(SSS_points(3,1)+SSS_points(3,2))/SSS_points(3,2);
    end
    
    cd z1_DSS
    a2_launch
    close all
    cd ..
    ss_z1_grid(it_z1_,1)=(B_ss+N_ss)/N_ss;
    
end

%%

myfig=figure(78);
set(myfig, 'Position', [0 0 600 600])
plot(z1_grid,sss_z1_grid(:,1),' *','Color',[1,0.1,0.1],'linewidth',2);
hold on
plot(z1_grid,sss_z1_grid(:,2),' *','Color',[0.5,0.5,0.5],'linewidth',2);
plot(z1_grid,sss_z1_grid(:,3),' *','Color',[0,0.5,0.1],'linewidth',2);
plot(z1_grid,ss_z1_grid(:,1),'-','Color',[0.1,0.3,1],'linewidth',2);
plot([0.72,0.72],[-1000,+1000],'k--','linewidth',1);
title('Stochastic steady states for different values of $z_1$', 'interpreter','latex','FontSize',14);
xlabel('$z_1$','interpreter','latex','FontSize',12)
ylabel('leverage ($K/N$)','interpreter','latex','FontSize',12)
legend({'LL-SSS','Unstable SSS','HL-SSS','DSS'},'Location','southeast', 'interpreter','latex','FontSize',10)
ylim([0,3])
ax = gca;
ax.YAxis.Exponent = 0;
grid

print -dpdf h78_SSS_z1
savefig(myfig,'h78_SSS_z1.fig');

title(' ', 'interpreter','latex','FontSize',14);
print -dpdf g78_SSS_z1

