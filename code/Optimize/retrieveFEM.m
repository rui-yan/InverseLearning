function FEM = retrieveFEM(param_in, order, hmax)


    assert(length(param_in) == 2*order);
    %% parameters
    N = 2;
    filename = 'LFP50_P2.mat';
    epsilon = 1e-2;
    stepsize = 1e-3;
    C = struct('alpha',2.50e3, 'beta', 765, 'gamma', 2.59e3, 'omega',858); 
    paramx = struct('order', order, 'coefficients', param_in(1:order));
    paramy = struct('order', order, 'coefficients', param_in(order+1:2*order));
    rhs_param = struct('C',C, 'stepsize',stepsize,...
        'filename', filename, 'epsilon', epsilon, ...
        'paramx',paramx, 'paramy',paramy);   
    % define c matrix for general pde
    bPDE = @(location, state) rhs(location, state, rhs_param);
    fPDE = @(location, state) -1.*rhs(location, state, rhs_param);
    %% make pde model
    model = createpde(N);
    % geometry
    tri = loadGeometry(filename);
    geometryFromMesh(model, tri.Points', tri.ConnectivityList');
    msh = generateMesh(model, 'Hmax', 250*hmax);

    % general pde
    %c = create_c(C);
    c = 1;
    specifyCoefficients(model, 'm', 0, 'd', 0, 'a',0, 'c',c, 'f',fPDE,'Face',1);

    % boundary condition
    force = applyBoundaryCondition(model, 'neumann', 'Edge',...
        1:model.Geometry.NumEdges, 'g', bPDE,'q',[0;0;0;0]);
    
    % get FEMmatrix
    fprintf("assembly FEM matrix\n");
    FEM = assembleFEMatrices(model,'nullspace');
end