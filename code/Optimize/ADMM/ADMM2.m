%CME307 HW3-P9, ADMM Part2
%% initialization
close all;clear all;clc;
rng(1);
N = 30; M = 10;
A = randn(M,N);
x_true = randn(N,1);
x_true = x_true - min(x_true) + 1e-5;
b = A * x_true;
V = orth(randn(N));
D = diag(abs(randn(N,1)));
Q = V\D*V;

beta = 1;
max_iter = 1e3;
x0 = abs(randn(N,1)); y10 = rand(M,1); y20 = rand(N,1);

%% generate true solutions, using intermal matlab function
x_min = quadprog(Q,zeros(N,1),-eye(N),zeros(N,1),A,b);

%% d) 5-blocks
[x,y1,y2,x_hist] = block_update_p(x0,y10,y20,Q,A,b,beta,max_iter);
x_norm = vecnorm(x_hist-x_min);
visualize_xnorm(x_norm);
legend('seq');
set(gca,'YScale','log');
plot_beautify;
saveas(gcf,'ADMMd.pdf');

%% e) randomly permute
[xp,y1p,y2p,xp_hist] = block_update_pperm(x0,y10,y20,Q,A,b,beta,max_iter);
xp_norm = vecnorm(xp_hist-x_min);
visualize_xnorm(xp_norm);
legend('permute');
set(gca,'YScale','log');
plot_beautify;
saveas(gcf,'ADMMe.pdf');

%% f) sample-without-replacement
[xr,y1r,y2r,xr_hist] = block_update_prand(x0,y10,y20,Q,A,b,beta,max_iter);
xr_norm = vecnorm(xr_hist-x_min);
visualize_xnorm(xr_norm);
legend('rswr');
set(gca,'YScale','log');
plot_beautify;
saveas(gcf,'ADMMf.pdf');

%% Summary
figure(4);hold on;
visualize_xnorm(x_norm);
visualize_xnorm(xp_norm);
visualize_xnorm(xr_norm);
legend('seq','permute','rswr');
set(gca,'YScale','log');
plot_beautify;
saveas(gcf,'ADMM_all.pdf');