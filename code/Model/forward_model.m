clear all;close all;
%% project forward model0
% pde parameter
N = 2;
model = createpde(N);

% mechanical parameters
alpha = 2.4947e3; beta = 765; gamma = 2.5897e3; omega = 858;
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

%% visualize shape
% pdegplot(model,'EdgeLabels','on','FaceLabels','on');
% xlim([-1.5,1.5]);
% axis equal;


%% set up the model
% equation
specifyCoefficients(model, 'm', 0, 'd', 0, 'a',0, 'c',c, 'f',rhs,'Face',1);

% boundary condition
force = applyBoundaryCondition(model, 'neumann', 'Edge',...
    1:model.Geometry.NumEdges, 'g', [1;0],'q',[0;0;0;0]);

% mesh generation
hmax = 1e-2;
msh = generateMesh(model,'Hmax',hmax,'GeometricOrder','quadratic');


%% solve
tic
fprintf("pde formed, solving pde...\n");
result = solvepde(model);
toc

%% result visualization
rs = result.NodalSolution;
varsToPlot = char('u, m', 'v, m');
for i = 1:size(varsToPlot,1)
    figure;
    pdeplot(model, 'XYData', rs(:,i), 'Contour', 'on');
    title(varsToPlot(i,:))
    % scale the axes to make it easier to view the contours
    axis square
end
