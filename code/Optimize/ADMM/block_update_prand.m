function [x,y1,y2,x_hist] = block_update_prand(x0,y10,y20,Q,A,b,beta,max_iter)
% @brief: ADMM to find minimum solution to 1/2xTQx 
%         using random sampling update,
% Transformation: Q = VDV^T, z = V^T*x, D is diag(N); OP becomes:
% min 1/2 z^TDz, s.t. (A*V)z = b, Vz>=0, now we can update separate blocks
% @param: x0,y0, initial guess
% @param: Q, quadratic matrix
% @param: A,b, constraint matrix and vector
% @param: beta, augmented lagrangian coefficient
% @param: max_iter, maximum iteration number
% @return: x,y1,y2, final output of contrl var, dual var1, dual var2
% @return: x_hist, x during iteration

    rng(0);
    
    % decomposition
    [M,N] = size(A);
    [V, D] = eig(Q);
    A = A*V;
    x0 = V'*x0;
    
    % matrix/vector blocks
    % x
    x = zeros(5,6);
    xx = x0;
    
    % y
    y1 = y10;
    y2 = y20;
    
    % s
    s = abs(randn(N,1));
    
    % A_block & D block & V block initialization
    A_block = cell(5,1);
    V_block = cell(5,1);
    p_block = cell(5,1);
    D_block = zeros(5,6);
    d = diag(D);

    x_hist = zeros(N,max_iter);
    for k = 1:max_iter
        % permute order 
        p = randperm(30);
        for i = 1:5
            p_block{i,1} = p(6*(i-1)+1:6*i);
            A_block{i,1} = A(:,p_block{i,1});
            V_block{i,1} = V(:,p_block{i,1});
            D_block(i,:) = d(p_block{i,1});
        end
    
        % x
        for i = 1:5
            x(i,:) = xx(p_block{i,1});
        end
        
        % update
        for i = 1:5
            r = 0;
            rv = 0;
            for j = 1:5
                if j ~= i
                    r = r + A_block{j,1} * x(j,:)';
                    rv = rv + V_block{j,1} * x(j,:)';
                end
            end
            
            x(i,:) = (diag(D_block(i,:)) + beta * A_block{i,1}' * A_block{i,1}...
                + beta * V_block{i,1}' * V_block{i,1})...
                \(A_block{i,1}' * y1 + V_block{i,1}' * y2 - beta * A_block{i,1}' * (r-b)...
                - beta * V_block{i,1}' * (rv-s));
        end
        
        % restore order
        xT = x';
        for i = 1:length(p)
            xx(p(i)) = xT(i);
        end
        
        % s
        s = V*xx - 1/beta * y2;
        s(s<0) = 0;
        
        % y
        y1 = y1 - beta * (A*xx - b);
        y2 = y2 - beta * (V*xx - s);
        x_hist(:,k) = V*xx;
    end
    x = V*xx;
end