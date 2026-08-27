% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

z1_list = [72 77 82 87 92 97];

it_z1_=0;

myfig_=figure(1);
set(myfig_, 'Position', [0 0 800 600])

mycolors = [1.0 0.1 0.1 ; 1.0 0.3 0.3 ; 1.0 0.5 0.5 ; 0.7 0.7 1.0 ; 0.4 0.4 1.0 ; 0.1 0.1 1.0 ];

for it_z1 = z1_list
    
    it_z1
    it_z1_=it_z1_+1;
    
    if it_z1 ~= 72
        eval(['load ''./z1_0' num2str(it_z1) '/z_FinalWorkspace.mat'';']);
    else
        load '../../4_sigma/sigma_300/z_FinalWorkspace.mat'
    end
    
    subplot(3,2,1)
    plot(BB_grid,distrib_B_,'-','Color',mycolors(it_z1_,:)','Linewidth',2)
    hold on
    
    subplot(3,2,2)
    plot(NN_grid,distrib_N_,'-','Color',mycolors(it_z1_,:)','Linewidth',2)
    hold on

    Ksim = Bsim+Nsim;
    [y_ax,x_ax] = hist(Ksim,60);
    dx = x_ax(2)-x_ax(1);
    y_ax = y_ax / sum(y_ax*dx);
    subplot(3,2,3:6)
    plot(x_ax,y_ax,'-','Color',mycolors(it_z1_,:)','Linewidth',2)
    hold on

end

subplot(3,2,1)
xlabel('(a) Debt ($B$)', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid
    
subplot(3,2,2)
xlabel('(b) Equity ($N$)', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid
    
subplot(3,2,3:6)
xlabel('(c) Capital ($K$)', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid
legend({['$z_1=' num2str(z1_list(1)/100) '$'],['$z_1=' num2str(z1_list(2)/100) '$'],['$z_1=' num2str(z1_list(3)/100) '$'],['$z_1=' num2str(z1_list(4)/100) '$'],['$z_1=' num2str(z1_list(5)/100) '$'],['$z_1=' num2str(z1_list(6)/100) '$']},'Location','best', 'interpreter','latex','FontSize',12)
    
print -dpdf h79_z1_distr
savefig(myfig_,'h79_z1_distr.fig');

print -dpdf g79_z1_distr

