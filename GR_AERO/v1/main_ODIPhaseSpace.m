% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
% Initial verification steps for the LD hypothesis
% Stability investigation of the One Damper Inoperative condition for the
% realistic GR model
%
% This code requires the use of the script phasespaceplot
close all; clear all; clc
%% Plotting parameters
setPlot
set(0,"defaultLineLineWidth",1)
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
INPUT_THETA_0 = 0.*pi/180; % radian
INPUT_A_1_H = 0.*pi/180; % radian
INPUT_B_1_H = 0.*pi/180; % radian

%======================= GENERAL HELICOPTER parameters=================
% Select number of blades
Nb = 4;         % [-], number of blades
dpsi = 2*pi/Nb;  % [rad], angle between blades
omega_0 = 40;     %[rad/s], rotor angular speed
nper = 150;      % [-], number of periods to run
e = 0.3048;     % [m], blade hinge offset from Hammond [1974]
% Activate or deactivate additional blade-hub damper in blades
% 0 : inactive
% 1 : active
% C_b2h = b2hdamp*C_xi(gamma~0.8)
b2hdamp = 1.;

x_0 = 0.1;
y_0 = 0.;
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
exec = 0.;

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
% damp = {'std' 'ib' 'i2b'};
damp = {'ib' 'i2b'};
% damp = {'std' 'i2b'};
% damp = {'i2b'};
% damp = {'std'};
% damp = {'ib'};
% damp = {'std' 'ib'};

%% MBDyn variables
geom = calculateKgeom(opt,e,damp,Nb);

qNR = NRdof(Nb);
% MBC degrees of freedom
if mod(Nb,2)==0
    dof_col = [1 Nb];
else
    dof_col = 1;
end
dof_cyc = setdiff(1:Nb,dof_col);

% omega = linspace(0.1,1.2,12)*omega_0;
omega = [250]*2*pi/60;
omega_plot = [1];

lblx = cell(1,6);
lbly = cell(1,6);

lblx{1} = '$x$ [m]'; lblx{2} = '$y$ [m]';
lbly{1} = '$\dot{x}$ [m/s]'; lbly{2} = '$\dot{y}$ [m/s]';
for nb = 1 : Nb
    lblx{2+nb} = ['$\xi_' int2str(nb) '$ [rad]'];
    lbly{2+nb} = ['$\dot{\xi}_' int2str(nb) '$ [rad/s]'];
end
tc = 0;

suffix = 'LLO';
disp('Starting lead-lag only case')
phasespaceplot
disp('Lead-lag only done...')

INPUT_THETA_0 = 8*deg2rad; % radian 
INPUT_A_1_H = 0.*deg2rad; % radian 
INPUT_B_1_H = 0*deg2rad; % radian 
suffix = 'T0';
tc = tc+10;
disp('Starting colective only case')
phasespaceplot
disp('FF3 done...')

INPUT_THETA_0 = 8*deg2rad; % radian 
INPUT_A_1_H = 7.*deg2rad; % radian 
INPUT_B_1_H = 0*deg2rad; % radian 
suffix = 'T0A1';
tc = tc+10;
disp('Starting colective and A1 case')
phasespaceplot
disp('FF1 done...')

INPUT_THETA_0 = 8*deg2rad; % radian 
INPUT_A_1_H = 0.*deg2rad; % radian 
INPUT_B_1_H = 7*deg2rad; % radian 
suffix = 'T0B1';
tc = tc+10;
disp('Starting colective and B1 case')
phasespaceplot
disp('FF2 done...')

%%
savefigures