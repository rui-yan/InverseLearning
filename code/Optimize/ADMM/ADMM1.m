% CME307 HW3-P9, ADMM Part1
%% initialization
close all;clear all;clc;
beta = [.1,1,10];
A = [1 1 1; 1 1 2; 1 2 2];
max_iter = 1e4;
rng(1);
x0 = randn(3,1); y0 = rand(3,1);

%% a) beta = 0.1, 1, 10, zero objective, no difference: all diverge
Fig_a = figure(1);hold on;
for i = 1:length(beta)
    [x,y,xnorm] = ADMM_div(x0,y0,beta(i),A, max_iter);
    visualize_xnorm(xnorm);
end
legend('0.1','1','10');
set(gca,'YScale','log');
plot_beautify;
saveas(gcf,'ADMMa.pdf');

%% b) beta = 0.1, 1, 10, objective is 2-norm of 0.5*||x||^2
Fig_b = figure(2);hold on;
for i = 1:length(beta)
    [x,y,xnorm] = ADMM_obj(x0,y0,beta(i),A, max_iter);
    visualize_xnorm(xnorm);
end
legend('0.1','1','10');
set(gca,'YScale','log');
plot_beautify;
saveas(gcf,'ADMMb.pdf');

%% c) random permute, repeat a and b
Fig_c1 = figure(3);hold on;
i = 2;
[x,y,xnorm] = ADMM_divrand(x0,y0,beta(i),A, max_iter);
visualize_xnorm(xnorm);
legend('9a, \beta = 1');
set(gca,'YScale','log');
plot_beautify;
saveas(gcf,'ADMMc1.pdf');

Fig_d = figure(4);hold on;
[x,y,xnorm] = ADMM_objrand(x0,y0,beta(i),A, max_iter);
visualize_xnorm(xnorm);
legend('9b, \beta = 1');
set(gca,'YScale','log');
plot_beautify;
saveas(gcf,'ADMMc2.pdf');
