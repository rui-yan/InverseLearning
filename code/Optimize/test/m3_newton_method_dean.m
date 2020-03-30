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

for i = 1:maxiter
    fprintf('%d th iteration\n', i);

    % Compute first derivative and second derivative
    [this_d1, this_K1, this_G1] = forward_model(a, b, g, o, f_is_constant);
    inv_this_K1 = this_K1 \ I;

    % compute first derivative df
    if f_is_constant == 0
        dd1_da = - inv_this_K1 * A * inv_this_K1 * this_G1 + inv_this_K1 * g1;
        dd1_db = - inv_this_K1 * B * inv_this_K1 * this_G1 + inv_this_K1 * g2;
        dd1_dg = - inv_this_K1 * G * inv_this_K1 * this_G1 + inv_this_K1 * g3;
        dd1_do = - inv_this_K1 * O * inv_this_K1 * this_G1 + inv_this_K1 * g4;
    else
        dd1_da = - inv_this_K1 * A * inv_this_K1 * G1;
        dd1_db = - inv_this_K1 * B * inv_this_K1 * G1;
        dd1_dg = - inv_this_K1 * G * inv_this_K1 * G1;
        dd1_do = - inv_this_K1 * O * inv_this_K1 * G1;
    end

    df_da = -2 * d1.' * dd1_da + 2 * this_d1.' * dd1_da;
    df_db = -2 * d1.' * dd1_db + 2 * this_d1.' * dd1_db;
    df_dg = -2 * d1.' * dd1_dg + 2 * this_d1.' * dd1_dg;
    df_do = -2 * d1.' * dd1_do + 2 * this_d1.' * dd1_do;

    df = [df_da, df_db, df_dg, df_do];

    % compute second derivative H
    if f_is_constant == 0
        dd1_da_2 = - 2 * inv_this_K1 * A * dd1_da;
        dd1_db_2 = - 2 * inv_this_K1 * B * dd1_db;
        dd1_dg_2 = - 2 * inv_this_K1 * G * dd1_dg;
        dd1_do_2 = - 2 * inv_this_K1 * O * dd1_do;

        dd1_dadb = - inv_this_K1 * A * dd1_db - inv_this_K1 * B * dd1_da;
        dd1_dadg = - inv_this_K1 * A * dd1_dg - inv_this_K1 * G * dd1_da;
        dd1_dado = - inv_this_K1 * A * dd1_do - inv_this_K1 * O * dd1_da;
        dd1_dbdg = - inv_this_K1 * B * dd1_dg - inv_this_K1 * G * dd1_db;
        dd1_dbdo = - inv_this_K1 * B * dd1_do - inv_this_K1 * O * dd1_db;
        dd1_dgdo = - inv_this_K1 * G * dd1_do - inv_this_K1 * O * dd1_dg;

        df_da_2 = -2 * d1.' * dd1_da_2 + 2 * (dd1_da' * dd1_da) + 2 * this_d1' * dd1_da_2;
        df_db_2 = -2 * d1.' * dd1_db_2 + 2 * (dd1_db' * dd1_db) + 2 * this_d1' * dd1_db_2;
        df_dg_2 = -2 * d1.' * dd1_dg_2 + 2 * (dd1_dg' * dd1_dg) + 2 * this_d1' * dd1_dg_2;
        df_do_2 = -2 * d1.' * dd1_do_2 + 2 * (dd1_do' * dd1_do) + 2 * this_d1' * dd1_do_2;

        df_dadb = -2 * d1.' * dd1_dadb + 2 * (dd1_db' * dd1_da) + 2 * this_d1' * dd1_dadb;
        df_dadg = -2 * d1.' * dd1_dadg + 2 * (dd1_dg' * dd1_da) + 2 * this_d1' * dd1_dadg;
        df_dado = -2 * d1.' * dd1_dado + 2 * (dd1_do' * dd1_da) + 2 * this_d1' * dd1_dado;
        df_dbdg = -2 * d1.' * dd1_dbdg + 2 * (dd1_dg' * dd1_db) + 2 * this_d1' * dd1_dbdg;
        df_dbdo = -2 * d1.' * dd1_dbdo + 2 * (dd1_do' * dd1_db) + 2 * this_d1' * dd1_dbdo;
        df_dgdo = -2 * d1.' * dd1_dgdo + 2 * (dd1_do' * dd1_dg) + 2 * this_d1' * dd1_dgdo;
    else
        dd1_da_2 = - 2 * inv_this_K1 * A * dd1_da;
        dd1_db_2 = - 2 * inv_this_K1 * B * dd1_db;
        dd1_dg_2 = - 2 * inv_this_K1 * G * dd1_dg;
        dd1_do_2 = - 2 * inv_this_K1 * O * dd1_do;

        dd1_dadb = - inv_this_K1 * A * dd1_db - inv_this_K1 * B * dd1_da;
        dd1_dadg = - inv_this_K1 * A * dd1_dg - inv_this_K1 * G * dd1_da;
        dd1_dado = - inv_this_K1 * A * dd1_do - inv_this_K1 * O * dd1_da;
        dd1_dbdg = - inv_this_K1 * B * dd1_dg - inv_this_K1 * G * dd1_db;
        dd1_dbdo = - inv_this_K1 * B * dd1_do - inv_this_K1 * O * dd1_db;
        dd1_dgdo = - inv_this_K1 * G * dd1_do - inv_this_K1 * O * dd1_dg;

        df_da_2 = -2 * d1.' * dd1_da_2 + 2 * (dd1_da' * dd1_da) + 2 * this_d1' * dd1_da_2;
        df_db_2 = -2 * d1.' * dd1_db_2 + 2 * (dd1_db' * dd1_db) + 2 * this_d1' * dd1_db_2;
        df_dg_2 = -2 * d1.' * dd1_dg_2 + 2 * (dd1_dg' * dd1_dg) + 2 * this_d1' * dd1_dg_2;
        df_do_2 = -2 * d1.' * dd1_do_2 + 2 * (dd1_do' * dd1_do) + 2 * this_d1' * dd1_do_2;

        df_dadb = -2 * d1.' * dd1_dadb + 2 * (dd1_db' * dd1_da) + 2 * this_d1' * dd1_dadb;
        df_dadg = -2 * d1.' * dd1_dadg + 2 * (dd1_dg' * dd1_da) + 2 * this_d1' * dd1_dadg;
        df_dado = -2 * d1.' * dd1_dado + 2 * (dd1_do' * dd1_da) + 2 * this_d1' * dd1_dado;
        df_dbdg = -2 * d1.' * dd1_dbdg + 2 * (dd1_dg' * dd1_db) + 2 * this_d1' * dd1_dbdg;
        df_dbdo = -2 * d1.' * dd1_dbdo + 2 * (dd1_do' * dd1_db) + 2 * this_d1' * dd1_dbdo;
        df_dgdo = -2 * d1.' * dd1_dgdo + 2 * (dd1_do' * dd1_dg) + 2 * this_d1' * dd1_dgdo;
    end

    H = [df_da_2 df_dadb df_dadg df_dado;
         df_dadb df_db_2 df_dbdg df_dbdo;
         df_dadg df_dbdg df_dg_2 df_dgdo;
         df_dado df_dbdo df_dgdo df_do_2];

    % update parameters
    J = H\df';
    a_new = a - J(1); b_new = b - J(2);
    g_new = g - J(3); o_new = o - J(4);

    % terminate when gradient smaller than tol
    if (norm(J) < tol)
      break;
    end

    a = a_new; b = b_new; g = g_new; o = o_new; % update

    result = f(a, b, g, o, d1, f_is_constant);

    a_values(i) = a;
    b_values(i) = b;
    g_values(i) = g;
    o_values(i) = o;
    f_values(i) = result;
    df_values(i,:) = norm(H);

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
