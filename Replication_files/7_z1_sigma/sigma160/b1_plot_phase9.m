% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

myfig_=figure(1);
set(myfig_, 'Position', [0 0 800 600])

it_z1_=0;

letters = ['(a) ';'(b) ';'(c) ';'(d) ';'(e) ';'(f) ';'(g) ';'(h) ';'(i) ']

z1_set = [67 67 67 72 77 82 87 92 97]

for it_z1 = z1_set
    
    it_z1
    it_z1_=it_z1_+1;
    
    if it_z1 ~= 72
        eval(['load ''./z1_0' num2str(it_z1) '/z_FinalWorkspace.mat'';']);
    else
        load '../../4_sigma/sigma_300/z_FinalWorkspace.mat'
    end
    
    cd phase
        b1_phase
    cd ..
    
    subplot(3,3,it_z1_);
    plot(PLM0(:,1),PLM0(:,2),'-','Color',[0.1,0.3,1],'linewidth',2);
    hold on;
    plot(NPLM0(:,1),NPLM0(:,2),'--r','linewidth',2);
    plot(SSS_points(:,1),SSS_points(:,2),' Ok','linewidth',2)
    plot(B_ss,N_ss,' s','Color',[0,0.5,0.1],'linewidth',2,'markers',6)
    mytitle=[letters(it_z1_,:) '$z_1=' num2str(it_z1/100) '$'];
    title(mytitle, 'interpreter','latex','FontSize',12);
    xlabel('debt ($B$)', 'interpreter','latex','FontSize',12);
    ylabel('equity ($N$)', 'interpreter','latex','FontSize',12);
    xlim([0 2.7])
    ylim([1.2 4.2])
    grid
    
end

legend({'$h(B,N)=0$','$\mu^N(B,N)=0$','SSS','DSS'},'Location','southeast', 'interpreter','latex','FontSize',10)

print -dpdf h74_phase9_z1
savefig(myfig_,'h74_phase6_z1.fig');

print -dpdf g74_phase9_z1

%%

myfig_=figure(2);
set(myfig_, 'Position', [0 0 800 600])

it_z1_=0;

for it_z1 = z1_set
    
    it_z1
    it_z1_=it_z1_+1;
    
    if it_z1 ~= 72
        eval(['load ''./z1_0' num2str(it_z1) '/z_FinalWorkspace.mat'';']);
    else
        load '../../4_sigma/sigma_300/z_FinalWorkspace.mat'
    end
    
    subplot(3,3,it_z1_);
    h=pcolor(squeeze(BB_grid),squeeze(NN_grid),PLM_visits');
    set(h, 'EdgeColor', 'none');
    mytitle=[letters(it_z1_,:) '$f(B,N)$ with $z_1=' num2str(it_z1/100) '$'];
    title(mytitle, 'interpreter','latex','FontSize',12);
    xlabel('debt ($B$)', 'interpreter','latex','FontSize',12);
    ylabel('equity ($N$)', 'interpreter','latex','FontSize',12);
    xlim([0 2.7])
    ylim([1.2 4.2])
    grid
    
end

%%

print -dpdf h76_zone9_z1
savefig(myfig_,'h76_zone6_z1.fig');

print -dpdf g76_zone9_z1
