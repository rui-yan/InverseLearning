clc
clear
format long

%% Obtain simulated "real" data
% Note that this should be changed to experiment data later
f_is_constant = 1;

% First obtain our simulated "real" data
% mechanical parameters
alpha = 2.4947e3; beta = 765; gamma = 2.5897e3; omega = 858;
[d1, K1, G1] = forward_model(alpha, beta, gamma, omega, f_is_constant);

[~, K_alpha_at_1, ~] = forward_model(1, 0, 0, 0, f_is_constant);
[~, K_beta_at_1, ~]  = forward_model(0, 1, 0, 0, f_is_constant);
[~, K_gamma_at_1, ~] = forward_model(0, 0, 1, 0, f_is_constant);
[~, K_omega_at_1, ~] = forward_model(0, 0, 0, 1, f_is_constant);

A = K_alpha_at_1; B = K_beta_at_1; G = K_gamma_at_1; O = K_omega_at_1;

I = eye(size(K1, 1)); % used to compute transpose

% Now assume given "real"/"experimental"/"simulation" data, we want to find
% estimated mechanical parameters that minimize our objective function

% Objective funtion (defined in file f.m)

%% Steepest Descent Algorithm:
% Initialization
%==========================================================================
% initialize alpha-->a, beta-->b, gamma-->g and omege-->o
a0 = 2.0e3; b0 = 700; g0 = 2.0e3; o0 = 800;
maxiter = 5000; % maximum number of allowed iterations
tol = 1e-6; % tolerance
step_size = 0.5; % step size
scale = 0.8; % scaling step size

% storing and tracking values
a_values = zeros(maxiter,1);
b_values = zeros(maxiter,1);
g_values = zeros(maxiter,1);
o_values = zeros(maxiter,1);

f_values = zeros(maxiter,1);
df_values = zeros(maxiter,4);

a = a0; b = b0; g = g0; o = o0;

for i = 1:maxiter
    fprintf('%d th iteration\n', i);

    % compute df
    [this_d1, this_K1, ~] = forward_model(a, b, g, o, f_is_constant);
    inv_this_K1 = this_K1 \ I;
    dd1_da = - inv_this_K1 * A * inv_this_K1 * G1;
    dd1_db = - inv_this_K1 * B * inv_this_K1 * G1;
    dd1_dg = - inv_this_K1 * G * inv_this_K1 * G1;
    dd1_do = - inv_this_K1 * O * inv_this_K1 * G1;

    df_da = -2 * d1.' * dd1_da + 2 * this_d1.' * dd1_da;
    df_db = -2 * d1.' * dd1_db + 2 * this_d1.' * dd1_db;
    df_dg = -2 * d1.' * dd1_dg + 2 * this_d1.' * dd1_dg;
    df_do = -2 * d1.' * dd1_do + 2 * this_d1.' * dd1_do;

    J = [df_da, df_db, df_dg, df_do];

    dir = -J; % search direction

%     % backward searching step-size
    while f(a - step_size * J(1), b - step_size * J(2), ...
            g - step_size * J(3), o - step_size * J(4), d1, f_is_constant) ...
            > f(a, b, g, o, d1, f_is_constant) - step_size/2 * norm(J)^2
        step_size = step_size * 0.8;
    end

    % update parameters
    a_new = a + step_size * dir(1);
    b_new = b + step_size * dir(2);
    g_new = g + step_size * dir(3);
    o_new = o + step_size * dir(4);

    % terminate if ||[a_new b_new g_new o_new] - [a b g o]|| < tol
%     if (norm([a_new, b_new, g_new, o_new] - [a, b, g, o]) < tol)
%       break;
%     end

    % terminate when gradient smaller than tol
    if (norm(step_size * dir) < tol)
      break;
    end

    a = a_new; b = b_new; g = g_new; o = o_new; % update

    result = f(a, b, g, o, d1, f_is_constant);

    a_values(i) = a;
    b_values(i) = b;
    g_values(i) = g;
    o_values(i) = o;
    f_values(i) = result;
    df_values(i,:) = J;

    fprintf('f(x) = %f\n', result);
end


% Output:
figure()
subplot(3,1,1)
y1 = f_values - result;
plot(y1, 'r', 'LineWidth', 2);
xlabel('Number of iteration');
title('Subplot 1: f(x^k) - f^*');

subplot(3,1,2)
y2 = vecnorm(df_values');
plot(y2, 'b', 'LineWidth', 2);
xlabel('Number of iteration');
title('Subplot 2: ||\nabla f(x^k)||');

subplot(3,1,3)
y3 = vecnorm([a_values, b_values, g_values, o_values]' - [a, b, g, o]');
plot(y3, 'g', 'LineWidth', 2);
xlabel('Number of iteration');
title('Subplot 3: ||x^k - x^*||');
sgtitle('Initialization using solutions of SDP')
