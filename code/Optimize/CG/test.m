% test
clc;close all;clear all;
global print_level;
print_level = 1;
%% 1. generate real data
% define range & params
xrange = linspace(1,220,101);
[X,Y] = meshgrid(xrange);
querypoints = [X(:),Y(:)]';
order = 2;
param_in = [ 0,0.0250,0,-0.0100];
% generate data from forward model
e_true = forward_chem(querypoints, param_in, order);

%% 2. optmization
% initialize
param0 = [0,0.02,0,0];
iter_max = 10;
alpha = 1e-4;

param = param0;
err = zeros(2*order*iter_max+1,1);
err(1) = sum(param - param_in).^2;
param_track = zeros(2*order*iter_max+1,2*order);
param_track(1,:) = param;
% batch gradient descent
% for i = 1:iter_max
%     param_grad = dfda(e_true, querypoints, param_in, order);
%     param = param - alpha .* param_grad;
%     err(i+1) = sum(param - param_in).^2;
%     param_track(i+1,:) = param;
%     fprintf("error %d: %3.4f", i, err(i));
% end

% single gradient: it seems that the number of nodes is very important...
% for i = 1:iter_max
%     for orderi = 1:2*order
%         param_grad = dfda_single(e_true, querypoints, param_in, order, orderi);
%         param = param - alpha .* param_grad;
%         err(2*order*(i-1)+1+orderi) = sum(param - param_in).^2;
%         param_track(2*order*(i-1)+1+orderi,:) = param;
%         fprintf("error %d: %3.4f", i, log(err(i)));
%     end
% end

% plot(1:2*order*iter_max+1, log(err));


%% CG:
acc = 1e-4;
maxfn = 1; 
dfpred = 3e-3; 
n = 2*order;
a0 = param0';
param_helper = struct('e_true',e_true,...
                      'querypoints',querypoints,...
                      'order',order);
[F,a] = cgrelax('cg_fchem',n,acc,maxfn,dfpred,a0,param_helper);
