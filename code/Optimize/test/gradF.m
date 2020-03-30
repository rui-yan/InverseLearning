function df = gradF(a, b, g, o, f_is_constant)
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
end