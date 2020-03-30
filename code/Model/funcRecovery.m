function F = funcRecovery(basis, coeffs, x)
% @brief: constructs function with basis and coeffs, evaluated at x
% @param: basis, basis function handle, f(i,x), 
%               i: series degree, x: variable
% @param: coeffs, coefficients 
% @param: x, variable
% @return: F, function value

    B = zeros(length(x),length(coeffs));
    for i  = 1:length(coeffs)
        B(:,i) = basis(i,x)';
    end
    F = B*coeffs;
end

