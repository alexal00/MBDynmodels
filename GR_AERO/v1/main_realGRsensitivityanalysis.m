% Author: Alejandro Alvaro, 2023-2024
% execute a generic rotor model and evaluta the evolution of the damping
% ratio under various simulation conditions.
%
% Objective: Assess the influence of generic pitch inputs in the damping
% ratios of the cyclic MBC. The appearance of the kinematic couplings
% caused by the proposed arrengement are expected to modify the loads in
% the damper and potentially the GR stability margins by altering the
% damping ratio.
%
% NOTE: This script takes several hours to run (~10 hours) as it investigates a broad
% range of operation conditions (~320 simulations)
%
% CAVEAT: The file writes/overwrites the matlab file 3DSA.mat. This file
% contains the results of the main script execution. UNLESS something has
% been chaged in the MBDyn files of the overall scripts refer to ppSA to
% simply plot the results obtained.
%
% This code requires the use of the script pp
close all; clear all; clc
%% Plotting parameters
setPlot
fid = 1;
lastfid = fid;
%% Helicopter data and enviroment variables
deg2rad = pi/180;
%================= FLIGHT CONDITIONS ============================
% % set MBDyn input for hover
V_INF = 0.; % m/s
ALTITUDE = 0.; % m
% NOTE: for convenience, mount the rotor with vertical shaft and put gamma + tau in alpha_H
ALPHA_H = 0.*pi/180; % radian
x_0 = 0.05;
y_0 = 0.;
%======================= GENERAL HELICOPTER parameters=================
% Select number of blades
Nb = 4;         % [-], number of blades
dpsi = 2*pi/Nb;  % [rad], angle between blades
omega = 40;     %[rad/s], rotor angular speed
nper = 60;      % [-], number of periods to run
e = 0.3048;     % [m], blade hinge offset from Hammond [1974]
% Activate or deactivate additional blade-hub damper in blades
% 0 : inactive
% 1 : active
% C_b2h = b2hdamp*C_xi(gamma~0.8)
b2hdamp = 1.;

% Use swashplate or not
sw = 1;

% Activate or deactivate blade degrees of freedom
% Since the degrees of freedom are given by means of a total joint in MBDyn
% 0 means that the degrre of freedom is not constrained and viceversa for
% 1.
% 0: allowed
% 1: clamped
if sw ==1
    pitch = 0;
else
    pitch = 1;
end
flap = 0;
lag = 0;


dof = '';
if pitch ==0
    dof = [dof 'p'];
end
if flap ==0
    dof = [dof 'f'];
end
if lag ==0
    dof = [dof 'l'];
end
% Execute or no MBDyn
% 1 : execute MBDyn
% else: do not execute
exec = 1.;

% Geometrical consideration
opt = 1;

% CL consideration
cl = 'lve';
%% MBDyn filepaths
% Filepaths for MBDyn
% In preffix
pref_in = '/home/ale/thesis/';
folder = 'GR_AERO/v1/';
model = 'rotor_';

% Out preffix
% pref_out = '/home/ale/thesis/AERO_MBDyn/aerodynamic_sw/';
% Shell file for MBDyn
shfile = '/home/ale/mbdyn/mbdyn/mbdyn ';


% Select damping type
damp = {'std' 'ib' 'i2b'};
% damp = {'ib' 'i2b'};
% damp = {'std' 'i2b'};
% damp = {'i2b'};
% damp = {'std'};
% damp = {'ib'};
% damp = {'std' 'ib'};

%% MBDyn variables
geom = calculateKgeom(opt,e,damp,Nb);
%% Influence of collective
% Input the cyclic conditions for the sensititivy analysis
INPUT_THETA_0 = (0:1:10).*pi/180; % radian
INPUT_A_1_H = 0.*pi/180; % radian
INPUT_B_1_H = 0.*pi/180; % radian
% INPUT_A_1_H = (0:1:7).*pi/180; % radian
% INPUT_B_1_H = (0:1:7).*pi/180; % radian

% Define the meshgrid for 3D plots
% [X, Y] = meshgrid(INPUT_A_1_H,INPUT_B_1_H);

% Declare variables for 3D plots
% zeta1cib = zeros(length(INPUT_A_1_H),length(INPUT_B_1_H),length(INPUT_THETA_0));
% zeta1sib = zeros(length(INPUT_A_1_H),length(INPUT_B_1_H),length(INPUT_THETA_0));
% 
% zeta1ci2b = zeros(length(INPUT_A_1_H),length(INPUT_B_1_H),length(INPUT_THETA_0));
% zeta1si2b = zeros(length(INPUT_A_1_H),length(INPUT_B_1_H),length(INPUT_THETA_0));

% indexes where the ib and i2b are stored in pp script
idib = find(strcmp(damp,'ib'));
idi2b = find(strcmp(damp,'i2b'));

% Initialize the progress bar
totalIterations = length(INPUT_THETA_0) * length(INPUT_B_1_H) * length(INPUT_A_1_H);
% progressBar = waitbar(0, 'Processing...');
iteration = 0; % Counter to keep track of progress

for zz = 1 : length(INPUT_THETA_0)
    for yy = 1 : length(INPUT_B_1_H)
        for xx = 1 : length(INPUT_A_1_H)
            variables =['"' 'const real OMEGA_100 = ' num2str(omega) ...
            '; const real N_PER = ' num2str(nper) ...
            '; const real V_INF = ' num2str(V_INF) ...
            '; const real ALTITUDE = ' num2str(ALTITUDE) ...
            '; const real ALPHA_H = ' num2str(ALPHA_H) ...
            '; const real INPUT_THETA_0 = ' num2str(INPUT_THETA_0(zz)) ...
            '; const real INPUT_A_1_H = ' num2str(INPUT_A_1_H(yy)) ...
            '; const real INPUT_B_1_H = ' num2str(INPUT_B_1_H(xx)) ...
            '; const real ACTIVE = ' num2str(b2hdamp) ...
            '; const integer PITCH = ' num2str(pitch) ...
            '; const integer FLAP = ' num2str(flap) ...
            '; const integer LAG = ' num2str(lag) ...
            '; const real A_BAR = ' num2str(geom.a) ...
            '; const real B_BAR = ' num2str(geom.b) ...
            '; const real C_BAR = ' num2str(geom.cd) ...
            '; const real C_BARA = ' num2str(geom.ca) ...
            '; const real C_BARB = ' num2str(geom.cb) ...
            '; const real D_BAR = ' num2str(geom.d) ...
            '; const real F_BAR = ' num2str(geom.f) ...
            '; const real GEOM = ' num2str(geom.Kxidelta^2) ...
            '; const real GEOMIB = ' num2str(geom.Kxil^2) ...
            '; const real x_0 = ' num2str(x_0) ...
            '; const real y_0 = ' num2str(y_0) ...
            '"'];
        
            setenv('MBDYNVARS', variables);
            % Post-processing script
            pp
            
            % Calculate ratios
            % \zeta_{\xi_{1c}}
            zeta1cib(zz) = zeta_cyc(idib,1)./zeta_cyc(1,1);
            % \zeta_{\xi_{1s}}
            zeta1sib(zz) = zeta_cyc(idib,2)./zeta_cyc(1,2);

            % \zeta_{\xi_{1c}}
            zeta1ci2b(zz) = zeta_cyc(idi2b,1)./zeta_cyc(1,1);
            % \zeta_{\xi_{1s}}
            zeta1si2b(zz) = zeta_cyc(idi2b,2)./zeta_cyc(1,2);

            % Update iteration count
            iteration = iteration + 1;
            
            % Update the progress bar
            disp(sprintf('Progress: %d%%', (100 * iteration / totalIterations)));

        end
    end


    % Contour plots
    % Find the minimum values of all subplots for CB axis
    mn(1) = min(min(zeta1cib2(:,:)));
    mn(2) = min(min(zeta1sib2(:,:)));
    mn(3) = min(min(zeta1ci2b2(:,:)));
    mn(4) = min(min(zeta1si2b2(:,:)));
    minColorLimit = min(mn);
    % Find the maximum value of all subplots for CB axis
    Mx(1) = max(max(zeta1cib2(:,:)));
    Mx(2) = max(max(zeta1sib2(:,:)));
    Mx(3) = max(max(zeta1ci2b2(:,:)));
    Mx(4) = max(max(zeta1si2b2(:,:)));
    maxColorLimit = max(Mx);


    fig(fid) = figure(fid);
    fig(fid).Name = ['Ratcyct0' int2str(10)];
    tiledlayout(2,2); fid = fid+1;

    ax = nexttile; % IB \xi1c
    contourf(X*180/pi,Y*180/pi,zeta1cib2(:,:,1))
    xlabel('$A_{1,H}$ [deg]')
    ylabel('$B_{1,H}$ [deg]')
    title(['ib ' qNR{2}])
    clim(ax,[minColorLimit,maxColorLimit]);
    colormap(ax,'jet')

    ax = nexttile; % IB \xi1s
    contourf(X*180/pi,Y*180/pi,zeta1sib2(:,:,1))
    xlabel('$A_{1,H}$ [deg]')
    ylabel('$B_{1,H}$ [deg]')
    title(['ib ' qNR{3}])
    clim(ax,[minColorLimit,maxColorLimit]);
    colormap(ax,'jet')

    ax = nexttile; % I2B \xi1c
    contourf(X*180/pi,Y*180/pi,zeta1ci2b2(:,:,1))
    xlabel('$A_{1,H}$ [deg]')
    ylabel('$B_{1,H}$ [deg]')
    title(['i2b ' qNR{2}])
    clim(ax,[minColorLimit,maxColorLimit]);
    colormap(ax,'jet')

    ax = nexttile; % I2B \xi1s
    contourf(X*180/pi,Y*180/pi,zeta1si2b2(:,:,1))
    xlabel('$A_{1,H}$ [deg]')
    ylabel('$B_{1,H}$ [deg]')
    title(['i2b ' qNR{3}])
    clim(ax,[minColorLimit,maxColorLimit]);
    colormap(ax,'jet')

    c = colorbar;  % attach colorbar to h
    c.Layout.Tile = 'east';
    c.TickLabelInterpreter = 'latex';
    c.Label.String = '\zeta/\zeta^{std}';
    colormap(c,'jet')

end

fig = figure(fid);
fig.Name = 'RatcycTheta0';
plot(INPUT_THETA_0*180/pi,zeta1cib,'DisplayName',['ib ' qNR{2}]); hold on
plot(INPUT_THETA_0*180/pi,zeta1sib,'DisplayName',['ib ' qNR{3}]);
plot(INPUT_THETA_0*180/pi,zeta1ci2b,'DisplayName',['i2b ' qNR{2}]);
plot(INPUT_THETA_0*180/pi,zeta1si2b,'DisplayName',['i2b ' qNR{3}]);
xlabel('$\theta_0$ [deg]')
ylabel('$\zeta / \zeta^{std}$')
title('$A_{1,H}$=0, $B_{1,H}$=0')
grid on
legend

%Save matrix results
save('3DSA_2.mat',"zeta1cib","zeta1sib","zeta1ci2b","zeta1si2b")
%%
return
savefigures
