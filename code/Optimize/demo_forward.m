clear all; close all;
%% parameters
N = 2;
filename = 'LFP50_P2.mat';
hmax = 5e-2;
xm_avg = 0.5;
epsilon = 1e-2;
stepsize = 1e-3;
C = struct('alpha',2.50e3, 'beta', 765, 'gamma', 2.59e3, 'omega',858); 
paramx = struct('order',2, 'coefficients',[0, 0.025]);
paramy = struct('order',2, 'coefficients',[0, -0.01]);
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
c = create_c(C);
specifyCoefficients(model, 'm', 0, 'd', 0, 'a',0, 'c',c, 'f',fPDE,'Face',1);

% boundary condition
force = applyBoundaryCondition(model, 'neumann', 'Edge',...
    1:model.Geometry.NumEdges, 'g', bPDE,'q',[0;0;0;0]);

%% test solution
tic
fprintf("pde formed, solving pde...\n");
result = solvepde(model);
rs = result.NodalSolution;
toc

% %% Visualization Check
% % define nodes for visualization
% xrange = linspace(1,220,101);
% [X,Y] = meshgrid(xrange);
% querypoints = [X(:),Y(:)]';
% [e_xx, e_xy, e_yy] = strainFromDisplacement(result, querypoints); 
% % solution visualization
% rsFig = visualizeRS(rs, model);
% % strain visualization
% hFig = plot_strain(e_xx, e_yy, e_xy, X, Y, filename, epsilon);

%% calculate target function

%% add constraint