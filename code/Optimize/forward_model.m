function [d1, K1, G1] = forward_model(alpha, beta, gamma, omega, f_is_constant)
    %% parameters
    % pde parameter
    N = 2;
    filename = 'LFP50_P2.mat';
%     hmax = 5e-2;
    hmax = 5;
%     xm_avg = 0.5;
    epsilon = 1e-2;
    stepsize = 1e-3;
    
    % mechanical parameters 
    % alpha = 2.4947e3; beta = 765; gamma = 2.5897e3; omega = 858;
    C = struct('alpha', alpha, 'beta', beta, 'gamma', gamma, 'omega', omega);
    paramx = struct('order',2, 'coefficients',[0, 0.025]);
    paramy = struct('order',2, 'coefficients',[0, -0.01]);
    rhs_param = struct('C', C, 'stepsize', stepsize,...
        'filename', filename, 'epsilon', epsilon, ...
        'paramx', paramx, 'paramy', paramy);
    
    % define c matrix for general pde
    if f_is_constant == 0
        bPDE = @(location, state) rhs(location, state, rhs_param);
        fPDE = @(location, state) -1.*rhs(location, state, rhs_param);
    else
        bPDE = -[1;0];
        fPDE = [1;0];
    end

    %% make pde model
    model = createpde(N);
    % geometry
    tri = loadGeometry(filename);
    geometryFromMesh(model, tri.Points', tri.ConnectivityList');
    msh = generateMesh(model, 'Hmax', 250*hmax);
    
    % general pde
    c = create_c(C);
    specifyCoefficients(model, 'm', 0, 'd', 0, 'a', 0, 'c', c, 'f', fPDE, 'Face', 1);
    
    % boundary condition
    force = applyBoundaryCondition(model, 'neumann', 'Edge',...
        1:model.Geometry.NumEdges, 'g', bPDE,'q',[0;0;0;0]);
    
    %% solve
%     tic
%     fprintf("pde formed, assemble FEM matrix\n");
    FEM = assembleFEMatrices(model, 'nullspace');
    K = FEM.Kc;
    K1 = K(3:end-1,3:end-1); 
    G = FEM.Fc;
    G1 = G(3:end-1);
    d1 = K1\G1;
%     toc
end