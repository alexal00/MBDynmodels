% Author: Alejandro Alvaro, 2023-2024
% Simple verification script for the proposed model used to verify
% the correct implementation of the non-linear constitutive law.
close all; clear all; clc
%% Plotting parameters
setPlot

%% ODE and parameter definition

% System parameters
omega_f = 2*pi;             % Forcing frequency, [rad/s]
K_xi = (10.*2*pi)^2;        % Stiffness, [Nm]
M_z = 1650;                 % Forcing amplitude, [Nm]
Jb = 1;                     % System inertia aroun z axis, [kgm^2]
omega_n = sqrt(K_xi)/Jb;    % Natural frequency of the system, [rad/s]

% NL viscous law =: f*tanh(alpha*eps_prime)
f = 100;                    % Non-linear viscoelastic coefficient
alpha = 10;                 % Multiplier


% Simulation conditions
T = 2*pi/omega_n;           % Oscilation period, [s]
steps = 256;                % Steps per period, [-]
dt = T/steps;               % Time-step, [s]
nper = 150;                 % Periods to run simulation, [-]
te = nper*T;                 % Final time, [s]
t_ode = 0:dt:te;            % Time vector for ODE

% Necessary parameters for the ODE stored inside a structure
params.Jb = Jb;
params.K_xi = K_xi;
params.f = f;
params.alpha = alpha;
params.M_z = M_z;
params.omega_f = omega_f;

%% ODE solution
[t_ode,theta_ode] = ode45(@(t,y) odefun(t,y,params),t_ode,[0. 0.]);

% Moments in the damper from the ODE
M_dode = K_xi*theta_ode(:,1)+f*tanh(alpha*theta_ode(:,2));

%% MBDyn execution

% Filepaths for MBDyn
% Folder of the file
pref_in = '/home/ale/thesis/FRICTION/';
% Model
model = 'rotsym_sf.mbd';

% Shell file for MBDyn
shfile = '/home/ale/mbdyn/mbdyn/mbdyn ';


% INPUT file
f_in = [pref_in model];

% Select the possible constitutive laws
% 30: Linear viscoelastic
% 40: Linear elastic, non-linear viscous
CL = [40];
for ii = 1:length(CL)
    % Environment variables for MBDyn
    variables = ['"' 'const integer CL = ' int2str(CL(ii)) ...
        '; const real OMEGA = ' num2str(omega_n) ...
        '; const real f = ' num2str(f) ...
        '; const real K_XI = ' num2str(K_xi) ...
        '; const real M_z = ' num2str(M_z) ...
        '; const real NPER = ' num2str(nper) ...
        '; const real OMEGA_F = ' num2str(omega_f) ...
        '"'];
    setenv('MBDYNVARS', variables);

    % Suffix for output file
    if CL(ii)==30
        suff = 'lve';
    else
        suff = 'lenlv';
    end
    % Output file name
    fn_base = ['hhh' suff];
    % Output file location
    f_out = [pref_in fn_base];

    % ----------------------------EXECUTE MBDYN
    disp('executing MBDyn...');
    [rc, errmsg] = system(['wsl ' 'MBDYNVARS=' variables ' ' shfile  f_in ' -o ' f_out]);
    
    % [rc, errmsg] = system(['./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
    if (rc ~= 0)
        error(errmsg);
    end
    disp('   ... done');
    % ----------------------------EXECUTE MBDYN
    % filename
    fn = [f_out '.nc'];
    disp(sprintf('reading output from file ''%s''', fn));
    

    time = ncread(fn, 'time');      % Time vector from MBDyn
    psi_nd = time*omega_n/(2*pi);   % N.D time in revs
    
    % Parameters for FFT
    nT_f = steps*omega_n/omega_f;           % Number of steps per one cycle of the forcing freq.
    revN_f = fix(te/(omega_f/(2*pi))) - 1;  % NUmber of periods in the forcing frequency, ~10
    
    % Extract data from .nc file
    def_hinge = 'elem.joint.31.';
    theta_h = ncread(fn, [def_hinge 'Theta']);
    thetap_h = ncread(fn, [def_hinge 'Omega']);
    M_d = ncread(fn,[def_hinge 'M'],[3, 1], [1, Inf])';

    % Figure 1: Theta and Theta_dot against time
    fh = figure(fid); fid = fid+1;
    fh.Name = 'thetavst';
    ax = subplot(1,2,1);
    functionplot(ax,psi_nd,theta_h,128,Color{ii},'o','none','MBDyn');
    plot(ax,t_ode/T,theta_ode(:,1),'Color','k','Marker','none','LineStyle',':','DisplayName','ODE')
    xlabel(ax,'revN [-]')
    ylabel(ax,'$\theta$')
    grid(ax,"on")

    ax = subplot(1,2,2);
    functionplot(ax,psi_nd,thetap_h,128,Color{ii},'o','none','MBDyn');
    plot(ax,t_ode/T,theta_ode(:,2),'Color','k','Marker','none','LineStyle',':','DisplayName','ODE')
    xlabel(ax,'revN [-]')
    ylabel(ax,'$\dot{\theta}$')
    grid(ax,"on")
    legend(ax)

    % Figure 2: Moments in the damper against time
    fh = figure(fid); fid = fid+1;
    fh.Name = 'Mvst';
    ax = gca;
    functionplot(ax,psi_nd,M_d,128,Color{ii},Marker{ii},'none','MBDyn');
    plot(ax,t_ode/T,M_dode,'Color','k','Marker','none','LineStyle',':','DisplayName','ODE')
    xlabel('revN [-]')
    ylabel('$M_d$ [Nm]')
    grid(ax,"on")
    legend(ax)

    % Truncate variables to only one period of the forcing frequency
    t_1f = time(nT_f*revN_f+[1:nT_f]);
    theta_1f = theta_h(nT_f*revN_f+[1:nT_f])*180/pi;
    Md_1f = M_d(nT_f*revN_f+[1:nT_f]);
    thetap_1f = thetap_h(nT_f*revN_f+[1:nT_f])*180/pi;


    % Constitutive Law
    % Figure 3: M_d vs theta and M_d vs theta_prime
    fh = figure(fid); fid = fid+1;
    fh.Name = 'constituetivelaw';
    ax = subplot(1,2,1);
    functionplot(ax,theta_1f,Md_1f,128,Color{ii},Marker{ii},'none','MBDyn');
    plot(ax,theta_ode(nT_f*revN_f+[1:nT_f],1)*180/pi,M_dode(nT_f*revN_f+[1:nT_f]),'Color','k','Marker','none','LineStyle',':','DisplayName','ODE')
    xlabel(ax,'$\theta$ [deg]')
    ylabel(ax,'$M_d$ [Nm]')
    grid(ax,"on")
    legend(ax)

    ax = subplot(1,2,2);
    functionplot(ax,thetap_1f,Md_1f,128,Color{ii},Marker{ii},'none','MBDyn');
    plot(ax,theta_ode(nT_f*revN_f+[1:nT_f],2)*180/pi,M_dode(nT_f*revN_f+[1:nT_f]),'Color','k','Marker','none','LineStyle',':','DisplayName','ODE')
    xlabel('$\dot{\theta}$ [deg/s]')
    ylabel('$M_d$ [Nm]')
    grid(ax,"on")
    legend(ax)

    % Computation of the FFT
    L = length(theta_1f);       % Sampled points
    Fs = 1/dt;                  % Sampling frequency
    df = Fs/L;                  % Increments of df for plotting
    f = df*(0:L/2);             % Generate frequency vector for FFT plot
    harm = f/(omega_n/(2*pi));
    [P1Theta,~] = computeFFT(theta_1f,t_1f);
    [P1ThetaP,~] = computeFFT(thetap_1f,t_1f);
    [P1MD,~] = computeFFT(Md_1f,t_1f);
    
    % Figure 4: FFT
    fh = figure(fid); fid = fid+1;
    fh.Name = 'FFT';
    tlo = tiledlayout('flow');
    nexttile;
    stem(harm,P1Theta);
    xlabel('$f$ [Hz]')
    ylabel('$|\theta|$ [deg]')
    xlim([0 2])
    nexttile;
    stem(harm,P1ThetaP);
    xlabel('$f$ [Hz]')
    ylabel('$|\dot{\theta}|$ [deg/s]')
    xlim([0 2])
    nexttile([1 2]);
    stem(harm,P1MD);
    xlabel('$f$ [Hz]')
    ylabel('$|M_d|$ [Nm]')
    xlim([0 2])


end