function L0 = legendre0(n,x)
% @brief: calculates legendre polynomial with order 0, P_iter^0
% @param: n, degree number
% @param: x, concentration, [-1,1]
% @return: L0, value Pn(2*c-1);

    L0 = legendreP(n,x);
end