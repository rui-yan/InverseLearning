% ADMM0
function [x,y,xnorm] = ADMM_obj(x0,y0,beta,A, max_iter)
% @brief: ADMM to solve the objective 0.5*(x1^2+x2^2+x3^2)
% @params: x0,y0, initial x0,y0 vector
% @params: beta, augmented constraint
% @params: A, linear constraint
% @params: max_iter, max iteration number
% @returns: x,y, final output
% @returns: x norm during the optmization process

    x = x0;
    y = y0;
    xnorm = zeros(max_iter,1);
    for iter = 1:max_iter
        b = beta * (x(2)*A(:,2) + x(3)*A(:,3))-y;
        x(1) = -(beta*A(:,1)'*A(:,1) + 1)\(A(:,1)'*b);
        b = beta * (x(1)*A(:,1) + x(3)*A(:,3))-y;
        x(2) = -(beta*A(:,2)'*A(:,2) + 1)\(A(:,2)'*b);
        b = beta * (x(2)*A(:,2) + x(1)*A(:,1))-y;
        x(3) = -(beta*A(:,3)'*A(:,3) + 1)\(A(:,3)'*b);
        y = y - beta*(A*x);
        xnorm(iter) =norm(x);
    end
end