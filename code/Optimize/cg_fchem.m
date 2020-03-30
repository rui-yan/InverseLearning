function [F, df] = cg_fchem(C, param_true, f_is_constant)
    
    d_true = param_true.d_true;
    K_true = param_true.K_true;
    G_true = param_true.G_true;
    
    a=C(1); b=C(2); g=C(3); o=C(4);
    
    [this_d1, this_K1, this_G1]  = forward_model(a, b, g, o, f_is_constant);
    F = norm(d_true - this_d1)^2;
    
    % warning
    [~, K_alpha_at_1, g1] =  forward_model(1, 0, 0, 0, f_is_constant);
    [~, K_beta_at_1,  g2] =  forward_model(0, 1, 0, 0, f_is_constant);
    [~, K_gamma_at_1, g3] =  forward_model(0, 0, 1, 0, f_is_constant);
    [~, K_omega_at_1, g4] =  forward_model(0, 0, 0, 1, f_is_constant);

    % B G O all zero sparse
    A = K_alpha_at_1; B = K_beta_at_1; G = K_gamma_at_1; O = K_omega_at_1;

    I = eye(size(this_K1, 1));
    
    inv_this_K1 = this_K1 \ I;
    if f_is_constant == 0
        dd1_da = - inv_this_K1 * A * inv_this_K1 * this_G1 + inv_this_K1 * g1;
        dd1_db = - inv_this_K1 * B * inv_this_K1 * this_G1 + inv_this_K1 * g2;
        dd1_dg = - inv_this_K1 * G * inv_this_K1 * this_G1 + inv_this_K1 * g3;
        dd1_do = - inv_this_K1 * O * inv_this_K1 * this_G1 + inv_this_K1 * g4;
    else
        dd1_da = - inv_this_K1 * A * inv_this_K1 * G_true;
        dd1_db = - inv_this_K1 * B * inv_this_K1 * G_true;
        dd1_dg = - inv_this_K1 * G * inv_this_K1 * G_true;
        dd1_do = - inv_this_K1 * O * inv_this_K1 * G_true;
    end
    
    df_da = -2 * d_true.' * dd1_da + 2 * this_d1.' * dd1_da;
    df_db = -2 * d_true.' * dd1_db + 2 * this_d1.' * dd1_db;
    df_dg = -2 * d_true.' * dd1_dg + 2 * this_d1.' * dd1_dg;
    df_do = -2 * d_true.' * dd1_do + 2 * this_d1.' * dd1_do;
   
    df = [df_da, df_db, df_dg, df_do];
end