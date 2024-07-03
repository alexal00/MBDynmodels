% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
% Initial verification steps for the LD hypothesis
% Stability investigation of the One Damper Inoperative condition for the
% realistic GR model
%
% This code requires the use of the script ppODI
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
INPUT_THETA_0 = 0.*pi/180; % radian
INPUT_A_1_H = 0.*pi/180; % radian
INPUT_B_1_H = 0.*pi/180; % radian

%======================= GENERAL HELICOPTER parameters=================
% Select number of blades
Nb = 4;         % [-], number of blades
dpsi = 2*pi/Nb;  % [rad], angle between blades
omega_0 = 40;     %[rad/s], rotor angular speed
nper = 60;      % [-], number of periods to run
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
if opt ==1 % original case
    geom.a = e/2;
    geom.b = e/2;
    geom.ca = e/4;
    geom.cb = e/4;
    geom.d = e/2;
    geom.f = e/2;
    geom.cd = e/4;
    geom.cf = e/4;
elseif opt==2 % modified case
    geom.a = e/4;
    geom.b = e/4;
    geom.ca = e/4;
    geom.cb = e/4;
    geom.d = e/4;
    geom.f = e/4;
    geom.cd = e/2;
    geom.cf = e/2;
end

geom.Kxidelta = 1;
geom.Kxil = 1;

if any(strcmp(damp, 'i2b'))
[geom.Kxidelta,k2,~]=i2bfun(e,geom.a,geom.b,geom.ca,geom.cb,...
    geom.d,geom.f,geom.cd,geom.cf,Nb);
    if geom.Kxidelta~=k2
        warning('Non-symetrical rotor, pay attention to geometrical factors')
    end
end
if any(strcmp(damp, 'ib'))
    [geom.Kxil,k2,~] = ibfun(e,geom.a,geom.b,geom.ca,geom.cb,Nb);
    if geom.Kxil~=-k2
        warning('Non-symetrical rotor, pay attention to geometrical factors')
    end
end
qNR = NRdof(Nb);
% MBC degrees of freedom
if mod(Nb,2)==0
    dof_col = [1 Nb];
else
    dof_col = 1;
end
dof_cyc = setdiff(1:Nb,dof_col);

omega = linspace(0.1,1.2,12)*omega_0;
suffix = 'LLO';
% omega = 26;
disp('Starting lead-lag only case')
ppODI
disp('Lead-lag only done...')

%% Collective and longitudinal cyclic
disp('Starting FF case with T0 and B1')
INPUT_THETA_0 = 10.*pi/180; % radian
INPUT_A_1_H = 0.*pi/180; % radian
INPUT_B_1_H = 7.*pi/180; % radian
suffix = 'T0B1';
ppODI
disp('FF with T0 and B1 done...')
%% Collective and longitudinal cyclic
disp('Starting FF case with T0 and A1')
INPUT_THETA_0 = 10.*pi/180; % radian
INPUT_A_1_H = 7.*pi/180; % radian
INPUT_B_1_H = 0.*pi/180; % radian
suffix = 'T0A1';
ppODI
%%
% return
savefigures
