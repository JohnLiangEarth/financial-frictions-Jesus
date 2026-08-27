% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

sigma_grid =nan(5,1);
z1_grid =nan(7,1);
sigma_z1_sss_grid = nan(4,7);
sigma_z1_leverage_grid = nan(4,7,3);


it_sigma_=0;

z1_set = [67 72 77 82 87 92 97]
sigma_set = [120 140 160 180 200 220 240 260 280 300]

for it_sigma = sigma_set
    
    it_z1_=0;
    
    it_sigma_ = it_sigma_+1;
    sigma_grid(it_sigma_)=it_sigma/10000;
    if it_sigma ~=140
        eval(['cd sigma' num2str(it_sigma)]);
    else
        cd ..
        cd 6_z1
    end
    
    for it_z1 = z1_set

        it_z1
        it_z1_=it_z1_+1;
        z1_grid(it_z1_)=it_z1/100;

        if it_z1 ~= 72
            eval(['load ''./z1_0' num2str(it_z1) '/z_FinalWorkspace.mat'';']);
        else
            if it_sigma ~=140
                eval(['load ''../../4_sigma/sigma_' num2str(it_sigma) '/z_FinalWorkspace.mat'';']);
            else
                load '../3_model_NN/z_FinalWorkspace.mat';
            end
        end
        cd 'phase';
        B3=0;
        N3=0;
        b1_phase
        close all
        cd ..

        sigma_z1_sss_grid(it_sigma_,it_z1_)=size(SSS_points,1);
        
        if size(SSS_points,1)==1
            if it_z1 < 85  % just because I know how it's going to look and it seems that this is the correct cutting value
                sigma_z1_leverage_grid(it_sigma_,it_z1_,1)=NaN;
                sigma_z1_leverage_grid(it_sigma_,it_z1_,2)=NaN;
                sigma_z1_leverage_grid(it_sigma_,it_z1_,3)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
            else
                sigma_z1_leverage_grid(it_sigma_,it_z1_,1)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
                sigma_z1_leverage_grid(it_sigma_,it_z1_,2)=NaN;
                sigma_z1_leverage_grid(it_sigma_,it_z1_,3)=NaN;
            end
        elseif size(SSS_points,1)==2
            sigma_z1_leverage_grid(it_sigma_,it_z1_,1)=NaN;
            sigma_z1_leverage_grid(it_sigma_,it_z1_,2)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
            sigma_z1_leverage_grid(it_sigma_,it_z1_,3)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
        else
            sigma_z1_leverage_grid(it_sigma_,it_z1_,1)=(SSS_points(1,1)+SSS_points(1,2))/SSS_points(1,2);
            sigma_z1_leverage_grid(it_sigma_,it_z1_,2)=(SSS_points(2,1)+SSS_points(2,2))/SSS_points(2,2);
            sigma_z1_leverage_grid(it_sigma_,it_z1_,3)=(SSS_points(3,1)+SSS_points(3,2))/SSS_points(3,2);
        end

    end
    
    if it_sigma ~=140
        cd ..
    else
        cd ..
        cd 7_z1_sigma
    end

end

%%

sss_one_LL = [];
sss_one_HL = [];
sss_two = [];

for it_sigma_ = 1:size(sigma_grid)
    for it_z1_ = 1:size(z1_grid)
        if sigma_z1_sss_grid(it_sigma_,it_z1_)>2
            sss_two=[sss_two ; [sigma_set(it_sigma_)/10000 z1_set(it_z1_)/100]];
        else
            if z1_set(it_z1_)/100 > 0.85  % just because I know how it's going to look and it seems that this is the correct cutting value
                sss_one_LL=[sss_one_LL ; [sigma_set(it_sigma_)/10000 z1_set(it_z1_)/100]];
            else
                sss_one_HL=[sss_one_HL ; [sigma_set(it_sigma_)/10000 z1_set(it_z1_)/100]];
            end
        end
    end
end


myfig=figure(81);
set(myfig, 'Position', [0 0 600 600])

plot(sss_one_LL(:,1),sss_one_LL(:,2),' o','Color',[1,0.1,0.1],'linewidth',6)
hold on
plot(sss_two(:,1),sss_two(:,2),' o','Color',[0.3,0.6,1],'linewidth',6)
plot(sss_one_HL(:,1),sss_one_HL(:,2),' o','Color',[0,0.5,0.1],'linewidth',6)

patch([0.01, 0.01, 0.032, 0.032, 0.017, 0.017],[0.895, 1.000, 1.000, 0.845, 0.845, 0.895],[1,0.1,0.1],'FaceAlpha',0.4,'EdgeColor','none')
patch([0.01, 0.01, 0.017, 0.017, 0.021, 0.021],[0.695, 0.895, 0.895, 0.845, 0.845, 0.695],[0.3,0.6,1],'FaceAlpha',0.4,'EdgeColor','none')
patch([0.01, 0.01, 0.021, 0.021, 0.032, 0.032],[0.600, 0.695, 0.695, 0.845, 0.845, 0.600],[0,0.5,0.1],'FaceAlpha',0.4,'EdgeColor','none')

title('Number of stable SSS found at each point', 'interpreter','latex','FontSize',12);
legend({'Only LL-SSS','HL-SSS and LL-SSS','Only HL-SSS'},'Location','best', 'interpreter','latex','FontSize',10)
xlabel('Aggregate risk ($\sigma$)', 'interpreter','latex','FontSize',12);
ylabel('Idiosyncratic risk ($z_1$)', 'interpreter','latex','FontSize',12);
xlim([0.0100 0.0320])
ylim([0.60   1.00  ])
grid

print -dpdf h81_SSS_sigma_z1
savefig(myfig,'h81_SSS_sigma_z1.fig');

title(' ', 'interpreter','latex','FontSize',14);
print -dpdf g81_SSS_sigma_z1




myfig=figure(82);
set(myfig, 'Position', [0 0 800 600])

for it_z1_ = 1:6

    subplot(2,3,it_z1_)
    plot(sigma_grid,squeeze(sigma_z1_leverage_grid(:,it_z1_,1)),' *','Color',[1,0.1,0.1],'linewidth',2);
    hold on
    plot(sigma_grid,squeeze(sigma_z1_leverage_grid(:,it_z1_,2)),' *','Color',[0.5,0.5,0.5],'linewidth',2);
    plot(sigma_grid,squeeze(sigma_z1_leverage_grid(:,it_z1_,3)),' *','Color',[0,0.5,0.1],'linewidth',2);
    % plot(z1_grid,ss_z1_grid(:,1),'-','Color',[0.1,0.3,1],'linewidth',2);
    % plot([0.72,0.72],[-1000,+1000],'k--','linewidth',1);
    title(['$z_1$=' num2str(z1_grid(it_z1_))], 'interpreter','latex','FontSize',14);
    xlabel('$\sigma$','interpreter','latex','FontSize',12)
    ylabel('leverage ($K/N$)','interpreter','latex','FontSize',12)
    ylim([0,3])
    ax = gca;
    ax.YAxis.Exponent = 0;
    grid

end

legend({'LL-SSS','Unstable SSS','HL-SSS','DSS'},'Location','best', 'interpreter','latex','FontSize',10)

print -dpdf h82_SSS_sigma_z1_6panels
savefig(myfig,'h82_SSS_sigma_z1_6panels.fig');

print -dpdf g82_SSS_sigma_z1_6panels



