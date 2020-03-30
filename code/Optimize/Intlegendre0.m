function L0 = Intlegendre0(iter,c)
% @brief: calculates legendre polynomial with order 0, P_iter^0
% @param: iter, degree number
% @param: c, concentration, [0,1]
% @return: L0, value Pn(2*c-1);
    if iter == 0
        L0 = ones(size(c));
    else
        L = legendre(iter,c);
        L0 = L(1,:);
        L0 = reshape(L0,size(c));
    end
end