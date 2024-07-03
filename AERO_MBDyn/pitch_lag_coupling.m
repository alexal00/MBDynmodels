% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
%
% Initial estimate to obtain the classic pitch-lag and flap-lag coupling.
% CAVEAT: This script is incomplete.
close all; clear all; clc
%% Plotting parameters
setPlot
fid = 1;
lastfid = fid;
%% Helicopter data and enviroment variables
deg2rad = pi/180;
%================= FLIGHT CONDITIONS ============================
% % set MBDyn input for standard HOVER flight
V_INF = 0.; % m/s 
ALTITUDE = 0.;% m, for standard air 
ALPHA_H = 0*deg2rad; % radian 
INPUT_THETA_0 = 6.54609*deg2rad; % radian 
INPUT_A_1_H = 0.*deg2rad; % radian 
INPUT_B_1_H = 0*deg2rad; % radian 

% Air properties
rho = 1.225;    % [kg/m^3], air density
%======================= GENERAL HELICOPTER parameters=================
% Select number of blades
Nb = 4;         % [-], number of blades

omega = 40;     %[rad/s], rotor angular speed
nper = 50;      % [-], number of periods to run
R = 6;        % [m], rotor radius
e = 1/25*R;     % [m], hinge position from hub
nuxi = sqrt(3/2*e/(R-e));   % [-], Rotating lag frequency
nubeta = sqrt(1+nuxi^2);    % [-], Rotating flap frequency
sigma = 0.08;               % [-], Solidity
c = sigma*pi*R/Nb;     % [m], chord
x_pitchhorn = 0.;
y_pitchhorn = 0.2;
delta3 = atan2(x_pitchhorn,y_pitchhorn);
% * inertia properties
Mb = 50.; % [kg]; define
Sb = 0.5*Mb*(R-e); % [kg m]; define
Ib = (Mb*(R-e)^2)/3; % [kg m^2]


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
exec = 1;

% Geometrical consideration
opt = 1;


%% MBDyn filepaths
% Filepaths for MBDyn
% In preffix
pref_in = '/home/ale/thesis/AERO_MBDyn/';
folder = 'aerodynamic_sw/';
model = 'rotor_';

% Out preffix
% pref_out = '/home/ale/thesis/AERO_MBDyn/aerodynamic_sw/';
% Shell file for MBDyn
shfile = '/home/ale/mbdyn/mbdyn/mbdyn ';


% Select damping type
% damp = {'std' 'ib' 'i2b'};
% damp = {'ib' 'i2b'};
% damp = {'std' 'i2b'};
% damp = {'i2b'};
damp = {'std'};
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

variables =['"' 'const real R = ' num2str(R) ...
    '; const real X_BLADE_FLAP_OFFSET = ' num2str(e) ...
    '; const real M_BLADE = ' num2str(Mb) ...
    '; const real OMEGA_100 = ' num2str(omega) ...
    '; const real N_PER = ' num2str(nper) ...
    '; const real V_INF = ' num2str(V_INF) ...
    '; const real ALTITUDE = ' num2str(ALTITUDE) ...
    '; const real ALPHA_H = ' num2str(ALPHA_H) ...
    '; const real INPUT_THETA_0 = ' num2str(INPUT_THETA_0) ...
    '; const real INPUT_A_1_H = ' num2str(INPUT_A_1_H) ...
    '; const real INPUT_B_1_H = ' num2str(INPUT_B_1_H) ...
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
    '; const real X_PITCH_HORN = ' num2str(x_pitchhorn)...
    '; const real Y_PITCH_HORN = ' num2str(y_pitchhorn)...
    '"'];
setenv('MBDYNVARS', variables);
%% Execute MBDyn
for kk=1:length(damp)
    if sw ==1
        f_in = [pref_in folder model damp{kk} '_Nb' int2str(Nb) '_sw.mbd'];
    else
        f_in = [pref_in folder model damp{kk} '_Nb' int2str(Nb) '_nsw.mbd'];
    end
    if b2hdamp==0
        fn_base = ['hhh' damp{kk} 'ud'];
    else
        fn_base = ['hhh' damp{kk} int2str(Nb) int2str(opt)];
    end
    f_out = [pref_in folder fn_base];
    
    if exec==1
        disp('executing MBDyn...');
    
    
        [rc, errmsg] = system(['wsl ' 'MBDYNVARS=' variables ' ' shfile  f_in ' -o ' f_out]);
    
        % [rc, errmsg] = system(['./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
        if (rc ~= 0)
            error(errmsg);
        end
        disp('   ... done');
    end
    
    % filename
    fn = [f_out '.nc'];
    disp(sprintf('reading output from file ''%s''', fn));
    
    
    % extract number of blades from database
    finfo = ncinfo(fn);
    nb = 0;
    while (1)
        nb = nb + 1;
        label = int2str(10000 + 1000*nb);
        n = 1;
        while (n <= size(finfo.Variables, 2))
            if (strcmp(finfo.Variables(n).Name, ['node.struct.',label,'.X']))
                break;
            end
            n = n + 1;
        end
        if (n > size(finfo.Variables, 2))
            break;
        end
    end
    
    nb = nb - 1;
    blade_txt = ['blade ', int2str(nb)];

    deltapsi = 2*pi/nb;
    %% Post-processing
    % Element labels
    blade_label = int2str(10000 + 1000*nb);         % Final blade label
    joint_label = int2str(10000 + 1000*nb + 30);    % Blade hinge label
    inflow_label = '99';                            % induced velocity or aerodynamic element
    
    % Get Omega from motion of the hub, assuming it starts with nominal RPM
    Omega = ncread(fn, 'node.struct.10000.Omega', [3, 1], [1, 1]);
    % compute azimuth vector from database
    psi_nd = ncread(fn, 'time')*Omega/(2*pi);
    
    nT = 256;
    revN = fix(length(psi_nd)/nT) - 1;
    
    % Induced velocity
    % vi = ncread(fn, ['elem.inducedvelocity.',inflow_label,'.UMean']);
    % Y = fft(vi(nT*revN+[1:nT]));
    % Y = Y/(nT/2);
    % vi0 = Y(1)/2;
    % lambdai0 = vi0/(omega*R);
    % N.D Induced velocity
    lambda = ncread(fn, ['elem.inducedvelocity.',inflow_label,'.Lambda']);
    Y = fft(lambda(nT*revN+[1:nT]));
    Y = Y/(nT/2);
    lambda0 = Y(1)/2;
    
    % Blade degrees of freedom
    blade_angles = ncread(fn, ['elem.joint.',joint_label,'.Phi']);
    
    for ii=1:size(blade_angles,1)
        if ii==1
            fact = 1;
        else
            fact = -1;
        end
    
        Y = fft(fact*blade_angles(ii,nT*revN+[1:nT]));   % FFT of the angles
        Y = Y/(nT/2);                               % Values are multiplied by the number of samples/2, nT/2
        Col(ii) = Y(1)/2;                           % Collective
        Cyc(ii,:) = [-real(Y(2)), +imag(Y(2))];   % Cyclic
        if ii~=3
            if ii==1
                A1c = 'A';
                B1s = 'B';
            else
                A1c = 'a';
                B1s = 'b';
            end
            disp(sprintf([A1c '_1_H = %g deg'], -real(Y(2))*180/pi));
            disp(sprintf([B1s '_1_H = %g deg'], +imag(Y(2))*180/pi));
        end
    end
    
    theta0 = Col(1);                % Collective pitch, should be equal to INPUT_THETA_0
    beta0 = Col(2);                 % Coning angle
    xi0 = Col(3);                  % Lead-lag angle
    

    fig(fid).Name = [damp{kk} 'PFL'];
    % fid = fid+1;
    % transform in deg, opposite sign for flap and lead-lag
    plot(psi_nd, 180/pi*ncread(fn, ['elem.joint.',joint_label,'.Phi'], [1, 1], [1, Inf]), ...
        psi_nd, -180/pi*ncread(fn, ['elem.joint.',joint_label,'.Phi'], [2, 1], [1, Inf]), ...
        psi_nd, -180/pi*ncread(fn, ['elem.joint.',joint_label,'.Phi'], [3, 1], [1, Inf]));
    xlabel('revolutions [-]');
    ylabel([blade_txt, ' angles, [deg]']);
    legend('$\theta$', '$\beta$', '$\xi$');
    grid on
    %--------------------------- REVISE THIS OPTIMIZATION ---------------------%
    % cla = 7.076;
    % cd0 = 0.008;
    % params.nubeta = nubeta;
    % params.nuxi = nuxi;
    % params.lamdda0 =lambda0;
    % params.theta0 = theta0;
    % params.beta0 = beta0;
    % params.xi0 = xi0;
    % params.INPUT_THETA_0 = INPUT_THETA_0;
    % params.delta3 = delta3;
    % lock = rho*cla*c*R^4/Ib; params.lock = lock;
    % fb = lock/8*(INPUT_THETA_0-4/3*lambda0); params.fb = fb;
    % fx = lock/8*(cd0/cla+4/3*(INPUT_THETA_0-3/2*lambda0)*lambda0); params.fx = fx;
    % 
    % delta2_1 = 1/(lock/8*xi0)*(nubeta^2*beta0+lock/8*beta0*tan(delta3)-fb);
    % delta2_2 = 1/(lambda0*lock*xi0/6)*(nuxi^2*xi0+lock/6*lambda0*beta0*tan(delta3)-fx);
    % X = getdelta3delta2(params);
    % theta_r = theta0+tan(X(1))*beta0-X(2)*xi0;
    %--------------------------- REVISE THIS OPTIMIZATION ---------------------%

    delta2 = (theta0-(INPUT_THETA_0-beta0*tan(delta3)))/xi0
    aux = 1;
    %
end