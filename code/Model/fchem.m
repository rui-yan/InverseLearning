function f = fchem(c, params)
% calculates chemical strain
% @param: c, concentration of Li
% @param: params, struct of order and coefficients
% @returns: chemical strain evaluated at c
    f = 0;
    for iter = 1:params.order
        f = f + params.coefficients(iter) .* legendre0(iter, 2*c-1);
    end
end