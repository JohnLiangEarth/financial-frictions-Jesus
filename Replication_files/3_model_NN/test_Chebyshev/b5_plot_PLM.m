% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution



% visited part of the PLM
PLM_surf_CH_v = PLM_surf_CH;
for iN=1:nval_NN
    for iB=1:nval_BB
       if isnan(PLM_visits(iB,iN))
           PLM_surf_CH_v(iB,iN)=NaN;      % put NaN on non-visited points
       end
    end
end



%%
%%%%%%%%%%%%%%%%%%%%%%%
% Now plot everything %
%%%%%%%%%%%%%%%%%%%%%%%


% PLM: h(B,N)

myfig=figure(94);
set(myfig, 'Position', [0 0 800 800])

mesh(BB_grid,NN_grid,PLM_surf_CH');
hold on
surf(BB_grid,NN_grid,PLM_surf_CH_v');
title('The perceived law of motion, $h(B,N)$ (Chebyshev approximation)', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
ylim([Nmin Nmax])
zlim([-0.1 0.1])
view([-150, 10])

print -dpdf h94_PLM_CH
savefig(myfig,'h94_PLM_CH.fig');

title(' ', 'interpreter','latex','FontSize',14);
print -dpdf g94_PLM_CH


% PLM: full range

myfig=figure(95);
set(myfig, 'Position', [0 0 800 800])

mesh(BB_grid,NN_grid,PLM_surf_CH');
hold on
surf(BB_grid,NN_grid,PLM_surf_CH_v');
title('The perceived law of motion, $h(B,N)$ (Chebyshev approximation)', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
ylim([Nmin Nmax])
view([60, 15])

print -dpdf h95_PLM_CH_fullrange
savefig(myfig,'h95_PLM_CH_fullrange.fig');

title(' ', 'interpreter','latex','FontSize',14);
print -dpdf g95_PLM_CH_fullrange


% PLM: difference CH-NN

PLM_surf_CH_dif_v=PLM_surf_CH_v-PLM_surf_v;

myfig=figure(96);
set(myfig, 'Position', [0 0 800 800])

surf(BB_grid,NN_grid,PLM_surf_CH_dif_v');
title('$h(B,N)$: difference between Chebyshev and NN approximations', 'interpreter','latex','FontSize',14);
xlabel('debt ($B$)', 'interpreter','latex','FontSize',14);
ylabel('equity ($N$)', 'interpreter','latex','FontSize',14);
xlim([Bmin Bmax])
ylim([Nmin Nmax])
view([60, 15])

print -dpdf h96_PLM_CH_dif
savefig(myfig,'h96_PLM_CH_dif.fig');

