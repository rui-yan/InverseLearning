clc;close all;clear all;
warning('off','all')
warning


a_t = 2.4947e3; b_t = 765; g_t = 2.5897e3; o_t = 858;

% generate data from forward model
[d_true, K_true, G_true] = forward_model(a_t, b_t, g_t, o_t);
param_helper = struct('d_true',d_true,...
                      'K_true',K_true,...
                      'G_true',G_true);

% fix b, o; Change a, g
% calculate df_da, df_dg.
pt=50;
temp=linspace(-2300, 2300, pt);
record_ag=zeros(pt);
gradient_a=zeros(pt,1);
gradient_g=zeros(pt,1);
for i=1:pt
    fprintf('1: current i is %d',i)
    a=a_t+temp(i);
    b=b_t;
    o=o_t;
    for j=1:pt
        fprintf('1: current j is %d',j)
        g=g_t+temp(j);
        C=[a,b,g,o]';
        [F, df]=cg_fchem(C, param_helper);
        record_ag(i,j)=F;
        gradient_a(i)=df(1);
        gradient_g(i)=df(3);
        fprintf('\n')
    end
     fprintf('\n')
end

format long;

save convex_plot1_pt50.mat temp record_ag


%contour(gradient_a/1e6',gradient_g/1e4',record_ag,'ShowText','on')
contour(temp,temp,log(record_ag),30,'ShowText','on')
title('fix o=o_t, b=b_t. log(norm(d_t-d^)^2)')
xlabel('a-a_t=linspace(-2300, 2300, pt)')
ylabel('g-g_t=linspace(-2300, 2300, pt)')



% fix b=b_t+24.3878, o=o_t+24.3878; Change a, g
% calculate df_da, df_dg.
record_ag2=zeros(pt);
gradient_a2=zeros(pt,1);
gradient_g2=zeros(pt,1);
for i=1:pt
    fprintf('2: current i is %d',i)
    a=a_t+temp(i);
    b=b_t+24.3878;
    o=o_t+24.3878;
    for j=1:pt
        fprintf('2: current j is %d',i)
        g=g_t+temp(j);
        C=[a,b,g,o]';
        [F, df]=cg_fchem(C, param_helper);
        record_ag2(i,j)=F;
        gradient_a2(i)=df(1);
        gradient_g2(i)=df(3);
        fprintf('\n')
    end
    fprintf('\n')
end

save convex_plot2_pt50.mat temp record_ag2

%contour(gradient_a/1e6',gradient_g/1e4',record_ag,'ShowText','on')
contour(temp,temp,log(record_ag2),30,'ShowText','on')
title('fix o=o_t+24.3878, b=b_t+24.3878. log(norm(d_t-d^)^2)')
xlabel('a-a_t=linspace(-2300, 2300, pt)')
ylabel('g-g_t=linspace(-2300, 2300, pt)')

%%%%%%%%%%%%%%%%%%%%
% fix b=b_t+300, o=o_t+300; Change a, g
% calculate df_da, df_dg.
record_ag3=zeros(pt);
gradient_a3=zeros(pt,1);
gradient_g3=zeros(pt,1);
for i=1:pt
    fprintf('3: current i is %d',i)
    a=a_t+temp(i);
    b=b_t+300;
    o=o_t+300;
    for j=1:pt
        fprintf('3: current j is %d',i)
        g=g_t+temp(j);
        C=[a,b,g,o]';
        [F, df]=cg_fchem(C, param_helper);
        record_ag3(i,j)=F;
        gradient_a3(i)=df(1);
        gradient_g3(i)=df(3);
        fprintf('\n')
    end
    fprintf('\n')
end

save convex_plot3_pt50.mat temp record_ag3

%contour(gradient_a/1e6',gradient_g/1e4',record_ag,'ShowText','on')
contour(temp,temp,log(record_ag3),30,'ShowText','on')
title('fix o=o_t+300, b=b_t+300. log(norm(d_t-d^)^2)')
xlabel('a-a_t=linspace(-2300, 2300, pt)')
ylabel('g-g_t=linspace(-2300, 2300, pt)')


%%%%%%%%%%%%%%%%%%%%%
% find whether a multiple of the true solns is also a minimizer.
mult=[0.5 1 1.5 2 5 10];
record_mult=zeros(0, size(mult,2));
for i=1:size(mult,2)
    C=mult(i)*[a_t,b_t,g_t,o_t]';
    [F, df]=cg_fchem(C, param_helper);
    record_mult(i)=F;
end
record_mult