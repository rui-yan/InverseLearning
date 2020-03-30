function [F, dFda] = cg_fchem(a, param)
    e_true = param.e_true;
    querypoints = param.querypoints;
    order = param.order;
    e_cal = forward_chem(querypoints, a, order);
    F = nansum((e_true - e_cal).^2);
    dFda = dfda(e_true, querypoints,a,order)';
end