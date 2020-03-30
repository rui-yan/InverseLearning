clear all
close all
clc

cmax = 2.29e4; % mol/m^3
kappa_scale = 0.022e-11/2.29e4;  % (Jm^2/mol)/(mol/m^3)
mu_scale = 298*8.3144598; % J/mol
F_scale = 96485.33289; % A.s/mol
c_scale = 2.29e3; % mol/m^3
k_scale = 1.e-2; % A/m^2
L_scale = sqrt(kappa_scale*c_scale/mu_scale); % m
T_scale = F_scale*c_scale*L_scale/k_scale; % s
iter = 100;
numberOfNodesInX = 201; % a axis
numberOfNodesInY = 501; % c axis
c_avg = zeros(iter,1);
time = zeros(iter,1);
boundaryPot = zeros(iter,1);
boundaryMu = zeros(iter,1);
mu_0 = 3.42;
totalNodes = numberOfNodesInY*numberOfNodesInX;
dof = 4;
strainDof = 8;


count = 100;
    
file = importdata(sprintf('ACRSimulationAvgConcentration_0_iteration_%d.txt',count));

data = file(2:totalNodes*dof+1);
boundaryPotential = file(totalNodes*dof+2);
totalTime = file(totalNodes*dof+3);
strains = file(totalNodes*dof+5:totalNodes*...
    (dof+strainDof)+4);

%% make nodal wise vectors

x = (0:1:numberOfNodesInX-1)*2/(numberOfNodesInX-1);
y = (0:1:numberOfNodesInY-1)*5/(numberOfNodesInY-1);
u = zeros(numberOfNodesInX,numberOfNodesInY);
v = zeros(numberOfNodesInX,numberOfNodesInY);
c = zeros(numberOfNodesInX,numberOfNodesInY);
mu = zeros(numberOfNodesInX,numberOfNodesInY);
eaa = zeros(numberOfNodesInX,numberOfNodesInY);
ecc = zeros(numberOfNodesInX,numberOfNodesInY);
eac = zeros(numberOfNodesInX,numberOfNodesInY);
time(count) = totalTime*T_scale;
boundaryPot(count) = boundaryPotential;

for k = 1:numberOfNodesInX
    for j = 1:numberOfNodesInY
        u(k,j) = data(dof*((j-1)*numberOfNodesInX+(k-1))+1);
        v(k,j) = data(dof*((j-1)*numberOfNodesInX+(k-1))+2);
        c(k,j) = data(dof*((j-1)*numberOfNodesInX+(k-1))+3)/(cmax/c_scale);
        mu(k,j) = data(dof*((j-1)*numberOfNodesInX+(k-1))+4);
        eaa(k,j) = strains(strainDof*((j-1)*numberOfNodesInX+(k-1))+1);
        ecc(k,j) = strains(strainDof*((j-1)*numberOfNodesInX+(k-1))+2);
        eac(k,j) = strains(strainDof*((j-1)*numberOfNodesInX+(k-1))+3);
    end
end

% plot contours of concentration, strain, etc
p = figure(1);
set(gcf,'color','w');
ax=subplot(3,2,2);
[C,h] = contourf(y,x,eaa-mean(mean(eaa)));
set(gca,'CLim', [-0.03,0.03],'fontsize',12);
set(h,'LineColor','none');
xlabel('c (\mum)')
ylabel('a (\mum)')
colormap(ax,polarmap)
title('e_{aa}')
colorbar('eastoutside')
axis equal


ax1=subplot(3,2,4)
[C1,h1] = contourf(y,x,ecc-mean(mean(ecc)));
set(gca,'CLim', [-0.03,0.03],'fontsize',12);
set(h1,'LineColor','none');
xlabel('c (\mum)')
ylabel('a (\mum)')
colormap(ax1,polarmap)
title('e_{cc}')
colorbar('eastoutside')
axis equal

ax2=subplot(3,2,6)
[C2,h2] = contourf(y,x,eac-mean(mean(eac)));
set(gca,'CLim', [-0.03,0.03],'fontsize',12);
set(h2,'LineColor','none');
xlabel('c (\mum)')
ylabel('a (\mum)')
colormap(ax2,polarmap)
title('e_{ac}')
colorbar('eastoutside')
axis equal

ax3=subplot(3,2,1)
[C3,h3] = contourf(y,x,c);
set(gca, 'CLim', [0.,1.],'fontsize',12);
set(h3,'LineColor','none');
xlabel('c (\mum)')
ylabel('a (\mum)')
title('x_{Li}')
colormap(ax3,multigradient([0 1 0; 1 1 0; 1 0 0]))
colorbar('eastoutside')
axis equal

% voltage relaxation
ax4 = subplot(3,2,[3,5])
plot(time(1:count),mu_0-boundaryPot(1:count)*mu_scale/F_scale,time(count)...
    ,mu_0-boundaryPot(count)*mu_scale/F_scale,'ro','linewidth',2)
axis([0,5531,3.405,3.425])
set(gca,'fontsize',12)
xlabel('t (s)')
ylabel('E (V vs Li^+/Li)')
