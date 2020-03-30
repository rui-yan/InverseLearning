function [d1, K1, G1] = forward_model(alpha, beta, gamma, omega)
    % pde parameter
    N = 2;
    model = createpde(N);

    % mechanical parameters 
    xm_avg = 0.5;
    c11 = [alpha, 0, omega/2];
    c12 = [0, beta; omega/2, 0];
    c22 = [omega/2, 0, gamma];
    c = [c11(:);c12(:);c22(:)];
    rhs = zeros(2,1);
    bc_rhs = zeros(2,1);

    %% create a rectangle geometry, to be modified later for actual shape
    rect= [3;4;-1;1;1;-1;-1;-1;1;1];
    gd = [rect];
    ns = char('rect');
    ns = ns';
    sf = 'rect';
    [dl,bt] = decsg(gd,sf,ns);
    geometryFromEdges(model,dl);

    %% set up the model
    % equation
    specifyCoefficients(model, 'm', 0, 'd', 0, 'a',0, 'c', c, 'f', rhs, 'Face', 1);

    % boundary condition
    force = applyBoundaryCondition(model, 'neumann', 'Edge',...
        1:model.Geometry.NumEdges, 'g', [10;0],'q',[0;0;0;0]);

    % mesh generation
    hmax = 50e-2;
    msh = generateMesh(model,'Hmax',hmax);

    %% solve
    tic
    fprintf("pde formed, assemble FEM matrix\n");
    FEM = assembleFEMatrices(model, 'nullspace');
    K = FEM.Kc;
    K1 = K(3:end-1,3:end-1); 
    G = FEM.Fc;
    G1 = G(3:end-1);
    d1 = K1\G1;

    toc

end