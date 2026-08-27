% Jesus Fernandez-Villaverde, Samuel Hurtado and Galo Nuno (2018)
% Financial Frictions and the Wealth Distribution

close all

load('./z1_097/z_FinalWorkspace.mat')

IRF_097_hlsss_Y = IRF_hlsss_Y;
IRF_097_hlsss_C = IRF_hlsss_C;
IRF_097_hlsss_Chat = IRF_hlsss_Chat;
IRF_097_hlsss_K = IRF_hlsss_K;
IRF_097_hlsss_B = IRF_hlsss_B;
IRF_097_hlsss_N = IRF_hlsss_N;
IRF_097_hlsss_w = IRF_hlsss_w;
IRF_097_hlsss_r = IRF_hlsss_r;
IRF_097_hlsss_rk = IRF_hlsss_rk;

load '../3_model_NN/z_FinalWorkspace.mat' IRF_hlsss_Y IRF_hlsss_C IRF_hlsss_Chat IRF_hlsss_K IRF_hlsss_B IRF_hlsss_N IRF_hlsss_w IRF_hlsss_r IRF_hlsss_rk IRF_llsss_Y IRF_llsss_C IRF_llsss_Chat IRF_llsss_K IRF_llsss_B IRF_llsss_N IRF_llsss_w IRF_llsss_r IRF_llsss_rk;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


myfig=figure(3);
set(myfig, 'Position', [0 0 800 800])

subplot(3,3,1);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_Y,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_Y,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(a) Output, $Y$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,2);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_C,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_C,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(b) HH consumption, $c$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,3);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_Chat,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_Chat,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(c) Expert consumption, $\hat{c}$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,4);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_K,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_K,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(d) Capital, $K$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,5);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_B,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_B,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(e) Debt, $B$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,6);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_N,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_N,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(f) Equity, $N$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,7);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_w,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_w,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(g) Wages, $w$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

legend({'$z_1$=0.72','$z_1$=0.97'},'Location','best', 'interpreter','latex','FontSize',10)

subplot(3,3,8);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_r,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_r,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(h) Interest rate, $r$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,9);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_rk,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_rk,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(i) Return to capital, $r_k$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

print -dpdf h52_IRF_072_097_ownSSS
savefig(myfig,'h52_IRF_072_097_ownSSS.fig');

print -dpdf g52_IRF_072_097_ownSSS






myfig=figure(4);
set(myfig, 'Position', [0 0 800 800])

subplot(3,3,1);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_Y,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_Y,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_Y,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(a) Output, $Y$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,2);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_C,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_C,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_C,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(b) HH consumption, $c$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,3);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_Chat,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_Chat,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_Chat,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(c) Expert consumption, $\hat{c}$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,4);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_K,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_K,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_K,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(d) Capital, $K$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,5);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_B,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_B,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_B,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(e) Debt, $B$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,6);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_N,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_N,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_N,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(f) Equity, $N$', 'interpreter','latex','FontSize',12);
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,7);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_w,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_w,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_w,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(g) Wages, $w$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

legend({'$z_1$=0.72, HL-SSS','$z_1$=0.72, LL-SSS','$z_1$=0.97'},'Location','best', 'interpreter','latex','FontSize',10)

subplot(3,3,8);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_r,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_r,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_r,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(h) Interest rate, $r$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

subplot(3,3,9);
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_hlsss_rk,'-','Color',[0,0.5,0.1],'linewidth',2);
hold on;
plot([1:size(IRF_hlsss_Y,1)]*dt, IRF_llsss_rk,'-.','Color',[1,0.1,0.1],'linewidth',1);
plot([1:size(IRF_097_hlsss_Y,1)]*dt, IRF_097_hlsss_rk,'--','Color',[0.3,0.3,1],'linewidth',1);
title('(i) Return to capital, $r_k$', 'interpreter','latex','FontSize',12);
xlabel('years','interpreter','latex','FontSize',12)
ax = gca;
ax.YAxis.Exponent = 0;
grid on

print -dpdf h53_IRF_072_097_all_SSS
savefig(myfig,'h53_IRF_072_097_all_SSS.fig');

print -dpdf g53_IRF_072_097_all_SSS


