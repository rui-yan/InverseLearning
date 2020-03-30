clc
clear
format long

%% Obtain simulated "real" data
% Note that this should be changed to experiment data later

% First obtain our simulated "real" data
% mechanical parameters
alpha = 2.4947e3; beta = 765; gamma = 2.5897e3; omega = 858;
[e_xx, e_xy, e_yy] = forward_model([alpha beta gamma omega]);

% Now assume given "real"/"experimental"/"simulation" data, we want to find 
% estimated mechanical parameters that minimize our objective function

%% Define objective funtion
% Our objective function is defined in the script titled objF.m
objF = @(e_xx, e_xy, e_yy, a, b, g, o) sumsqr(e_xx - ele_forward_model([a b g o], 1)) ...
    + sumsqr(e_xy - ele_forward_model([a b g o], 2)) ...
    + sumsqr(e_yy - ele_forward_model([a b g o], 3));

%% Steepest Descent Algorithm:
% Initialization
%==========================================================================
% Maximum number of allowed iterations
maxiter = 100;

objF_value = zeros(maxiter,1);
a_track = zeros(maxiter,1);
b_track = zeros(maxiter,1);
g_track = zeros(maxiter,1);
o_track = zeros(maxiter,1);
tol = 1e-4; % Tolerance
h = 0.5; % Step size (fixed) % required tuning

% tracking the values of alpha_hat, beta_hat, gamma_hat, omega_hat
a_track(1) = 2.0e3; b_track(1) = 500; g_track(1) = 2e3; o_track(1) = 1.5e3;

i = 1; % Iteration Counter

% Gradient Computation: since we cannot compute diff directly, compute it
% by the definition of gradient (not accurate though)

delta = 1e-4; % what delta should I chose?

df_da = ( objF(e_xx, e_xy, e_yy, a_track(i)+delta, b_track(i), g_track(i), o_track(i)) ...
    - objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)) ) / delta;
df_db = ( objF(e_xx, e_xy, e_yy, a_track(i), b_track(i)+delta, g_track(i), o_track(i)) ...
    - objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)) ) / delta;
df_dg = ( objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i)+delta, o_track(i)) ...
    - objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)) ) / delta;
df_do = ( objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)+delta) ...
    - objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)) ) / delta;

J = [df_da df_db df_dg df_do];

S = -(J); % Search Direction

% Minimization Condition:
while and(norm(J)>=tol, i <= maxiter)
    I = [a_track(i), b_track(i), g_track(i), o_track(i)]';
    
    % Compute objective function value and store it
    objF_value(i) = objF(e_xx, e_xy, e_yy,a_track(i),b_track(i),g_track(i),o_track(i));
    fprintf('The objective function has value: %d\n', objF_value(i))
   
%     update step size
%     g = subs(f, [a,b,g,o], [a_track(i)+h*S(1), b_track(i)+h*S(2), ...
%         g_track(i)+h*S(3), g_track(i)+h*S(4)]);
%     syms h; % step size
%     dg_dh = diff(g,h);
%     h = solve(dg_dh, h); % Optimal Step Length

    i = i+1;
    
    a_track(i) = I(1) + h * S(1); % Updated a value
    b_track(i) = I(2) + h * S(2); % Updated b value
    g_track(i) = I(3) + h * S(3); % Updated g value
    o_track(i) = I(4) + h * S(4); % Updated o value
    
    df_da = ( objF(e_xx, e_xy, e_yy, a_track(i)+delta, b_track(i), g_track(i), o_track(i)) ...
    - objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)) ) / delta;

    df_db = ( objF(e_xx, e_xy, e_yy, a_track(i), b_track(i)+delta, g_track(i), o_track(i)) ...
        - objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)) ) / delta;
    
    df_dg = ( objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i)+delta, o_track(i)) ...
        - objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)) ) / delta;
    
    df_do = ( objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)+delta) ...
        - objF(e_xx, e_xy, e_yy, a_track(i), b_track(i), g_track(i), o_track(i)) ) / delta;

    J = [df_da df_db df_dg df_do]; % Updated Gradient
    S = -(J); % New Search Direction
    
    
end

% Output:
fprintf('Initial Objective Function Value: %d\n\n', ...
    objF(e_xx, e_xy, e_yy,a_track(1),b_track(1),g_track(1),o_track(1)));
if (norm(J) < tol)
    fprintf('Minimum succesfully obtained...\n\n');
end
fprintf('Number of Iterations for Convergence: %d\n\n', i);
fprintf('Point of Minima: [%d,%d,%d,%d]\n\n\n\n', a_track(i),b_track(i),g_track(i),o_track(i));
fprintf('Objective Function Minimum Value Post-Optimization: %d\n\n', ...
    objF(e_xx, e_xy, e_yy, a_track(i),b_track(i),g_track(i),o_track(i)));
