% test grad2F
clc
clear
format long

%% Obtain simulated "real" data
% Note that this should be changed to experiment data later
f_is_constant = 0;

% First obtain our simulated "real" data
% mechanical parameters
alpha = 2.4947e3; beta = 765; gamma = 2.5897e3; omega = 858;
[d1, K1, G1] = forward_model(alpha, beta, gamma, omega, f_is_constant);

[~, K_alpha_at_1, g1] =  forward_model(1, 0, 0, 0, f_is_constant);
[~, K_beta_at_1,  g2] =  forward_model(0, 1, 0, 0, f_is_constant);
[~, K_gamma_at_1, g3] =  forward_model(0, 0, 1, 0, f_is_constant);
[~, K_omega_at_1, g4] =  forward_model(0, 0, 0, 1, f_is_constant);

A = K_alpha_at_1; B = K_beta_at_1; G = K_gamma_at_1; O = K_omega_at_1;

I = eye(size(K1, 1)); % used to compute transpose

% Now assume given "real"/"experimental"/"simulation" data, we want to find
% estimated mechanical parameters that minimize our objective function

% Objective funtion (defined in file f.m)

%% Newton Method:
% Initialization
%==========================================================================
% initialize alpha-->a, beta-->b, gamma-->g and omege-->o
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works

% Try different a0 [2.2e3, 2.4e3]
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 800;   % This one works

% a0 = 1.0e5; b0 = 700; g0 = 2.0e3; o0 = 800;   % This one works
% a0 = 1.0e4; b0 = 700; g0 = 2.0e3; o0 = 800;   % This one works
% a0 = 5.0e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.5e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.4e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.1e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one does not work
a0 = 100; b0 = 700; g0 = 2.0e3; o0 = 800; % This one does not work

% Try different b0 [600, 1100]
% a0 = 2.2e3; b0 = 1500; g0 = 2.0e3; o0 = 800; % This one does not work
% a0 = 2.2e3; b0 = 1200; g0 = 2.0e3; o0 = 800; % This one does not work
% a0 = 2.2e3; b0 = 1100; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 1000; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 900; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 800; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 600; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 500; g0 = 2.0e3; o0 = 800; % This one does not work
% a0 = 2.2e3; b0 = 400; g0 = 2.0e3; o0 = 800; % This one does not work

% Try different g0 [1.8e3, 2.4e3]
% a0 = 2.2e3; b0 = 700; g0 = 2.6e3; o0 = 800; % This one does not work
% a0 = 2.2e3; b0 = 700; g0 = 2.5e3; o0 = 800; % This one does not work
% a0 = 2.2e3; b0 = 700; g0 = 2.4e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 2.3e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 2.2e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 1.8e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 1.7e3; o0 = 800; % This one does not work
% a0 = 2.2e3; b0 = 700; g0 = 1.5e3; o0 = 800; % This one does not work

% Try different o0 [700, 800]
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 1000; % This one does not work
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 900; % This one does not work
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 800; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 700; % This one works
% a0 = 2.2e3; b0 = 700; g0 = 2.0e3; o0 = 500; % This one does not work

% Try the combination of lowest range for each parameter
% a0 = 2.2e3; b0 = 600; g0 = 1.8e3; o0 = 700; % Doesn't work :/

x0 = [a0, b0, g0, o0];

maxiter = 1000; % maximum number of allowed iterations
tol = 1e-6; % tolerance

% storing and tracking values
a_values = zeros(maxiter,1);
b_values = zeros(maxiter,1);
g_values = zeros(maxiter,1);
o_values = zeros(maxiter,1);

f_values = zeros(maxiter,1);
df_values = zeros(maxiter,4);

a = a0; b = b0; g = g0; o = o0;

delta = 1e-4;
H1 = grad2F(a, b, g, o, f_is_constant);
dF1 = gradF(a+delta*a, b, g, o, f_is_constant );
dF2 = gradF(a-delta*a, b, g, o, f_is_constant );
dF3 = gradF(a, b+delta*b, g, o, f_is_constant );
dF4 = gradF(a, b-delta*b, g, o, f_is_constant );
dF5 = gradF(a, b, g+delta*g, o, f_is_constant );
dF6 = gradF(a, b, g-delta*g, o, f_is_constant );
dF7 = gradF(a, b, g, o+delta*o, f_is_constant );
dF8 = gradF(a, b, g, o-delta*o, f_is_constant );
dF2delta = (dF1-dF2)./(2*delta*a);
dF2deltb = (dF3-dF4)./(2*delta*b);
dF2deltg = (dF5-dF6)./(2*delta*g);
dF2delto = (dF7-dF8)./(2*delta*o);