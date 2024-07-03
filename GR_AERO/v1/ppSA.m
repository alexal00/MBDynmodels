% Author: Alejandro Alvaro, 2023-2024
% This script is used to simply plot the results obtained from GRaeroSA. It
% requires the file 3DSA.mat ot exist and contain the corresponding damping
% ratios.
%
% CAVEAT: This file is created as a post-processing tool for the results as
% the main script takes several hours to run.
close all; clear all; clc;
%% Postprocessing Sensitivy Analysis
setPlot
fid = 10;
% Load the results
load("3DSA.mat")
% Labels for inputs
Nb = 4;
qNR = NRdof(Nb);
% Define the mesh
[X, Y] = meshgrid(INPUT_A_1_H.*180/pi,INPUT_B_1_H.*180/pi);
% Loop for plotting
for zz = 1 : length(INPUT_THETA_0)
    % Find the minimum values of all subplots for CB axis
    mn(1) = min(min(zeta1cib(:,:,zz)));
    mn(2) = min(min(zeta1sib(:,:,zz)));
    mn(3) = min(min(zeta1ci2b(:,:,zz)));
    mn(4) = min(min(zeta1si2b(:,:,zz)));
    minColorLimit = min(mn);
    % Find the maximum value of all subplots for CB axis
    Mx(1) = max(max(zeta1cib(:,:,zz)));
    Mx(2) = max(max(zeta1sib(:,:,zz)));
    Mx(3) = max(max(zeta1ci2b(:,:,zz)));
    Mx(4) = max(max(zeta1si2b(:,:,zz)));
    maxColorLimit = max(Mx);

    fig = figure(fid);
    fig.Name = ['Ratcyc' int2str(zz)]; fid = fid+1;
    tiledlayout(2,2)

    ax = nexttile; % IB \xi1c
    contourf(X,Y,zeta1cib(:,:,zz)')
    xlabel('$A_{1,H}$ [deg]')
    ylabel('$B_{1,H}$ [deg]')
    title(['ib ' qNR{2}])
    clim(ax,[minColorLimit,maxColorLimit]);
    colormap(ax,'jet')

    ax = nexttile; % IB \xi1s
    contourf(X,Y,zeta1sib(:,:,zz)')
    xlabel('$A_{1,H}$ [deg]')
    ylabel('$B_{1,H}$ [deg]')
    title(['ib ' qNR{3}])
    clim(ax,[minColorLimit,maxColorLimit]);
    colormap(ax,'jet')

    ax = nexttile; % I2B \xi1c
    contourf(X,Y,zeta1ci2b(:,:,zz)')
    xlabel('$A_{1,H}$ [deg]')
    ylabel('$B_{1,H}$ [deg]')
    title(['i2b ' qNR{2}])
    clim(ax,[minColorLimit,maxColorLimit]);
    colormap(ax,'jet')

    ax = nexttile; % I2B \xi1s
    contourf(X,Y,zeta1si2b(:,:,zz)')
    xlabel('$A_{1,H}$ [deg]')
    ylabel('$B_{1,H}$ [deg]')
    title(['i2b ' qNR{3}])
    clim(ax,[minColorLimit,maxColorLimit]);
    colormap(ax,'jet')

    % h = axes(fig,'visible','off');

    c = colorbar;  % attach colorbar to h
    c.Layout.Tile = 'east';
    c.TickLabelInterpreter = 'latex';
    c.Label.String = '\zeta/\zeta^{std}';
    colormap(c,'jet')
    % clim(h,[minColorLimit,maxColorLimit]);             % set colorbar limits
end

%% Influence of collective only
% Fix A1 and B1 to zero
fig = figure(fid); fid = fid+1;
fig.Name = 'RatcycTheta0';
plot(INPUT_THETA_0.*180/pi,squeeze(zeta1cib(1,1,:)),'DisplayName',['ib ' qNR{2}]); hold on
plot(INPUT_THETA_0.*180/pi,squeeze(zeta1sib(1,1,:)),'DisplayName',['ib ' qNR{3}]);
plot(INPUT_THETA_0.*180/pi,squeeze(zeta1ci2b(1,1,:)),'DisplayName',['i2b ' qNR{2}]);
plot(INPUT_THETA_0.*180/pi,squeeze(zeta1si2b(1,1,:)),'DisplayName',['i2b ' qNR{3}]);
xlabel('$\theta_0$ [deg]')
ylabel('$\zeta / \zeta^{std}$')
title('$A_{1,H}$=0, $B_{1,H}$=0')
grid on
legend('Location','best')
%% Influence of lateral cyclic
% Fix theta0 and B1 to zero
fig = figure(fid); fid = fid+1;
fig.Name = 'RatcycA1';
plot(INPUT_A_1_H.*180/pi,squeeze(zeta1cib(:,1,1)),'DisplayName',['ib ' qNR{2}]); hold on
plot(INPUT_A_1_H.*180/pi,squeeze(zeta1sib(:,1,1)),'DisplayName',['ib ' qNR{3}]);
plot(INPUT_A_1_H.*180/pi,squeeze(zeta1ci2b(:,1,1)),'DisplayName',['i2b ' qNR{2}]);
plot(INPUT_A_1_H.*180/pi,squeeze(zeta1si2b(:,1,1)),'DisplayName',['i2b ' qNR{3}]);
xlabel('$A_{1,H}$ [deg]')
ylabel('$\zeta / \zeta^{std}$')
title('$\theta_{0}$=0, $B_{1,H}$=0')
grid on
legend('Location','best')
%% Influence of longitudinal cyclic
% Fix theta0 and A1 to zero
fig = figure(fid); fid = fid+1;
fig.Name = 'RatcycB1';
plot(INPUT_B_1_H.*180/pi,squeeze(zeta1cib(1,:,1)),'DisplayName',['ib ' qNR{2}]); hold on
plot(INPUT_B_1_H.*180/pi,squeeze(zeta1sib(1,:,1)),'DisplayName',['ib ' qNR{3}]);
plot(INPUT_B_1_H.*180/pi,squeeze(zeta1ci2b(1,:,1)),'DisplayName',['i2b ' qNR{2}]);
plot(INPUT_B_1_H.*180/pi,squeeze(zeta1si2b(1,:,1)),'DisplayName',['i2b ' qNR{3}]);
xlabel('$B_{1,H}$ [deg]')
ylabel('$\zeta / \zeta^{std}$')
title('$\theta_{0}$=0, $A_{1,H}$=0')
grid on
legend('Location','best')

%% Influence of longitudinal cyclic in FF
% Fix theta0 and A1 to zero
fig = figure(fid); fid = fid+1;
fig.Name = 'RatcycT0B1FF';
tlo = tiledlayout('flow');
ax1 = nexttile(tlo);
hold(ax1,"on")
xlabel(ax1,'$B_{1,H}$ [deg]')
ylabel(ax1,'$\zeta_{1c}^{ib} / \zeta^{std}$')
title(ax1,['ib ' qNR{2}])
grid(ax1,"on")
% legend(ax1,'Location','best')

ax2 = nexttile(tlo);
hold(ax2,"on")
xlabel(ax2,'$B_{1,H}$ [deg]')
ylabel(ax2,'$\zeta_{1s}^{ib} / \zeta^{std}$')
title(ax2,['ib ' qNR{3}])
grid(ax2,"on")
% legend(ax2,'Location','best')

ax3 = nexttile(tlo);
hold(ax3,"on")
xlabel(ax3,'$B_{1,H}$ [deg]')
ylabel(ax3,'$\zeta_{1c}^{i2b} / \zeta^{std}$')
grid(ax3,"on")
title(ax3,['i2b ' qNR{2}])
% legend(ax3,'Location','best')

ax4 = nexttile(tlo);
hold(ax4,"on")
xlabel(ax4,'$B_{1,H}$ [deg]')
ylabel(ax4,'$\zeta_{1c}^{i2b} / \zeta^{std}$')
title(ax4,['i2b ' qNR{3}])
grid(ax4,"on")
% legend(ax4,'Location','best')
lg  = legend(ax4,'Orientation','Horizontal','NumColumns',6); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

for zz = 1 :length(INPUT_THETA_0)
plot(ax1,INPUT_B_1_H.*180/pi,squeeze(zeta1cib(1,:,zz)),'Color',Color{zz},...
    'DisplayName',['$\theta_0$=' num2str(INPUT_THETA_0(zz)*180/pi) ' deg']); hold on
plot(ax2,INPUT_B_1_H.*180/pi,squeeze(zeta1sib(1,:,zz)),'Color',Color{zz},...
    'DisplayName',['$\theta_0$=' num2str(INPUT_THETA_0(zz)*180/pi) ' deg']);
plot(ax3,INPUT_B_1_H.*180/pi,squeeze(zeta1ci2b(1,:,zz)),'Color',Color{zz},...
    'DisplayName',['$\theta_0$=' num2str(INPUT_THETA_0(zz)*180/pi) ' deg']);
plot(ax4,INPUT_B_1_H.*180/pi,squeeze(zeta1si2b(1,:,zz)),'Color',Color{zz},...
    'DisplayName',['$\theta_0$=' num2str(INPUT_THETA_0(zz)*180/pi) ' deg']);
end

%% Influence of lateral cyclic in FF
fig = figure(fid); fid = fid+1;
fig.Name = 'RatcycT0A1FF';
tlo = tiledlayout('flow');
ax1 = nexttile(tlo);
hold(ax1,"on")
xlabel(ax1,'$A_{1,H}$ [deg]')
ylabel(ax1,'$\zeta_{1c}^{ib} / \zeta^{std}$')
title(ax1,['ib ' qNR{2}])
grid(ax1,"on")
% legend(ax1,'Location','best')

ax2 = nexttile(tlo);
hold(ax2,"on")
xlabel(ax2,'$A_{1,H}$ [deg]')
ylabel(ax2,'$\zeta_{1s}^{ib} / \zeta^{std}$')
title(ax2,['ib ' qNR{3}])
grid(ax2,"on")
% legend(ax2,'Location','best')

ax3 = nexttile(tlo);
hold(ax3,"on")
xlabel(ax3,'$A_{1,H}$ [deg]')
ylabel(ax3,'$\zeta_{1c}^{i2b} / \zeta^{std}$')
grid(ax3,"on")
title(ax3,['i2b ' qNR{2}])
% legend(ax3,'Location','best')

ax4 = nexttile(tlo);
hold(ax4,"on")
xlabel(ax4,'$A_{1,H}$ [deg]')
ylabel(ax4,'$\zeta_{1s}^{i2b} / \zeta^{std}$')
title(ax4,['i2b ' qNR{3}])
grid(ax4,"on")
% legend(ax4,'Location','best')
lg  = legend(ax4,'Orientation','Horizontal','NumColumns',6); 
lg.Layout.Tile = 'South'; % <-- Legend placement with tiled layout

for zz = 1 :length(INPUT_THETA_0)
plot(ax1,INPUT_A_1_H.*180/pi,squeeze(zeta1cib(:,1,zz)),'Color',Color{zz},...
    'DisplayName',['$\theta_0$=' num2str(INPUT_THETA_0(zz)*180/pi) ' deg']); hold on
plot(ax2,INPUT_A_1_H.*180/pi,squeeze(zeta1sib(:,1,zz)),'Color',Color{zz},...
    'DisplayName',['$\theta_0$=' num2str(INPUT_THETA_0(zz)*180/pi) ' deg']);
plot(ax3,INPUT_A_1_H.*180/pi,squeeze(zeta1ci2b(:,1,zz)),'Color',Color{zz},...
    'DisplayName',['$\theta_0$=' num2str(INPUT_THETA_0(zz)*180/pi) ' deg']);
plot(ax4,INPUT_A_1_H.*180/pi,squeeze(zeta1si2b(:,1,zz)),'Color',Color{zz},...
    'DisplayName',['$\theta_0$=' num2str(INPUT_THETA_0(zz)*180/pi) ' deg']);
end
