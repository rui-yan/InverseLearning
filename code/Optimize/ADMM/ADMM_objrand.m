function [x,y,xnorm] = ADMM_objrand(x0,y0,beta,A, max_iter)
% @brief: ADMM to solve the objective 0.5*(x1^2+x2^2+x3^2), random
% order
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
        p = randperm(3);
        b = beta * (x(p(2))*A(:,p(2)) + x(p(3))*A(:,p(3)))-y;
        x(p(1)) = -(beta*A(:,p(1))'*A(:,p(1)) + 1)\(A(:,p(1))'*b);
        b = beta * (x(p(1))*A(:,p(1)) + x(p(3))*A(:,p(3)))-y;
        x(p(2)) = -(beta*A(:,p(2))'*A(:,p(2)) + 1)\(A(:,p(2))'*b);
        b = beta * (x(p(2))*A(:,p(2)) + x(p(1))*A(:,p(1)))-y;
        x(p(3)) = -(beta*A(:,p(3))'*A(:,p(3)) + 1)\(A(:,p(3))'*b);
        y = y - beta*(A*x);
        xnorm(iter) =norm(x);
    end
end