function [x,y1,y2,x_hist] = block_update_p(x0,y10,y20,Q,A,b,beta,max_iter)
% @brief: ADMM to find minimum solution to 1/2xTQx using 5-block update
% Transformation: Q = VDV^T, z = V^T*x, D is diag(N); OP becomes:
% min 1/2 z^TDz, s.t. (A*V)z = b, Vz>=0, now we can update separate blocks
% @param: x0,y0, initial guess
% @param: Q, quadratic matrix
% @param: A,b, constraint matrix and vector
% @param: beta, augmented lagrangian coefficient
% @param: max_iter, maximum iteration number
% @return: x,y1,y2, final output of contrl var, dual var1, dual var2
% @return: x_hist, x during iteration

    
    % decomposition
    [M,N] = size(A);
    [V, D] = eig(Q);
    A = A*V;
    x0 = V'*x0;
    
    % matrix/vector blocks
    % x
    x = zeros(5,6);
    for k = 1:5
        x(k,:) = x0(6*(k-1)+1:6*k);       
    end
    
    % y
    y1 = y10;
    y2 = y20;
    
    % s
    s = abs(randn(N,1));
    
    % A_block & D block & V block
    A_block = cell(5,1);
    V_block = cell(5,1);
    D_block = zeros(5,6);
    d = diag(D);
    for k = 1:5
        A_block{k,1} = A(:,6*(k-1)+1:6*k);
        V_block{k,1} = V(:,6*(k-1)+1:6*k);
        D_block(k,:) = d(6*(k-1)+1:6*k);
    end
    
    % block updates
    x_hist = zeros(N,max_iter);
    for k = 1:max_iter
        % x
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
        
        % s
        xx = [x(1,:)';x(2,:)';x(3,:)';x(4,:)';x(5,:)'];
        s = V*xx - 1/beta * reshape(y2,[N,1]);
        s(s<0) = 0;
        
        % y
        y1 = y1 - beta * (A*xx - b);
        y2 = y2 - beta * (V*xx - s);
        x_hist(:,k) = V*xx;
    end
    x = V*xx;
end