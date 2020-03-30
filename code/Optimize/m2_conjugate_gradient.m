% test
clc;close all;clear all;
warning('off','all');
global print_level;
print_level = 1;
%% 1. generate real data d
% define range & params
a_t = 2.4947e3; b_t = 765; g_t = 2.5897e3; o_t = 858;

f_is_constant = 0; % IF f is variable, f_is_constant = 0; otherwise 1.

% generate data from forward model
[d_true, K_true, G_true] = forward_model(a_t, b_t, g_t, o_t, f_is_constant);

%% 2. optmization
% initialize
a0 = 1.0e4; b0 = 700; g0 = 2.0e3; o0 = 800;

% Try different a0, obviously more robust regarding to the intial values
% a0 = 5.0e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works, 215 seconds
% a0 = 3.0e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works, 117 seconds
% a0 = 2.5e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works, 122 seconds

param0 = [a0,b0,g0,o0];

%iter_max = 10;

%param = param0;
order=2;
%err = zeros(2*order*iter_max+1,1);
%err(1) = norm(param - [a_t, b_t, g_t, o_t])^2;
%param_track = zeros(2*order*iter_max+1,2*order);
%param_track(1,:) = param;

%% CG:
acc = 1e-6;
maxfn = 5000;
dfpred = 3e5; 
n = 2*order;
C0 = param0';
param_helper = struct('d_true',d_true,...
                      'K_true',K_true,...
                      'G_true',G_true);
tic;                 
[F,a] = cgrelax('cg_fchem',n,acc,maxfn,dfpred,C0,param_helper, f_is_constant);
toc;