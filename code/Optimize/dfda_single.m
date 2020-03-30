function gradf = dfda_single(e_true, querypoints,param_in,order, order_idx)
    delta = 1e-8;
    gradf = zeros(1,2*order);
    for i = order_idx
        param_p = param_in;
        param_p(i) = param_p(i) + delta;
        param_m = param_in;
        param_m(i) = param_m(i) - delta;
        e_p = forward_chem(querypoints,param_p,order);
        e_m = forward_chem(querypoints,param_m,order);
        f_p = nansum((e_true - e_p).^2);
        f_m = nansum((e_true - e_m).^2);
        gradf(i) = (f_p - f_m)./(2*delta);
    end
end