% Objective funtion
function f_val = f(alpha, beta, gamma, omega, d1, f_is_constant)
    [d1_, ~, ~] = forward_model(alpha, beta, gamma, omega, f_is_constant);
    f_val = norm(d1 - d1_)^2;
end
