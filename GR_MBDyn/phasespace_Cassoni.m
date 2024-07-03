%% Code start
% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
% Initial verification steps for the LD hypothesis
% Stability investigation of the One Damper Inoperative condition for the
% realistic GR model
%
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
nper = 80;      % [-], number of periods to run
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
folder = 'GR_MBDyn/cassoni_model/';
model = 'GR_rotor.mbd';

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
qNR = NRdof(Nb);
% MBC degrees of freedom
if mod(Nb,2)==0
    dof_col = [1 Nb];
else
    dof_col = 1;
end
dof_cyc = setdiff(1:Nb,dof_col);

% omega = linspace(0.1,1.2,12)*omega_0;
omega = [125 250]*2*pi/60;
omega_plot = [1 2];

lblx = cell(1,6);
lbly = cell(1,6);

lblx{1} = '$x$ [m]'; lblx{2} = '$y$ [m]';
lbly{1} = '$\dot{x}$ [m/s]'; lbly{2} = '$\dot{y}$ [m/s]';
for nb = 1 : Nb
    lblx{2+nb} = ['$\xi_' int2str(nb) '$ [rad]'];
    lbly{2+nb} = ['$\dot{\xi}_' int2str(nb) '$ [rad/s]'];
end
% Test case numbering for figures
tc = 0;

suffix = 'LLO';
disp('Starting lead-lag only case')
for kk=1:length(damp)

    index = 1;
    for om = omega
    disp(['Omega is at:' num2str(om)])
    variables =['"' 'const real OMEGA_100 = ' num2str(om) ...
        '"'];
    setenv('MBDYNVARS', variables);
   
    % Execution of MBDyn and PP
    
    
        % Preprocessing of file names and testcase chosen
        if sw == 1
            swash = 'sw';
        else
            swash = 'nsw';
        end
    
        f_in = [pref_in folder model];
    
        suff = [damp{kk} int2str(Nb) int2str(opt)];
    
        fn_base = ['hhh' suff];
    
        f_out = [pref_in folder fn_base 'ODI'];
    
        fignames = [damp{kk} int2str(Nb) num2str(round(om*60/(2*pi),0))];
        % Execution of MBDyn
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
        
        
        % compute azimuth vector from database
        time = ncread(fn, 'time');
        dt = time(2)-time(1);
        fs = 1/dt;
        psi_nd = time*om/(2*pi);
    
        % Pp for fft and truncation of signals
        nT = 256;                           % Steps per period
        revN = fix(length(psi_nd)/nT) - 1;  % Number of revs -1
        
        % Perturbation applied between T_PERT=40 periods and 41.
        % Damper deactivated between T_PERT+1 (41) and 41*T_REV+dt.

        % Extract information about the x and y displacements
        % # 4 hub motions
        AIRFRAME_X = 100; 
        AIRFRAME_Y = 200;
        hub_x = ncread(fn, ['node.struct.',int2str(AIRFRAME_X),'.X'],[1, 1], [1, Inf]);
        hub_y = ncread(fn, ['node.struct.',int2str(AIRFRAME_Y),'.X'],[2, 1], [1, Inf]);
        v_x = ncread(fn,['node.struct.' num2str(AIRFRAME_X) '.XP'],[1,1],[1,Inf]);
        v_y = ncread(fn,['node.struct.' num2str(AIRFRAME_Y) '.XP'],[2,1],[1,Inf]);

        % Extract lead-lag angles and azimut for each blade for MBC
        % transformation
        % Read psi vector from the hub joint
        psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
        %Allocate vectors for angular position and velocity
        xi = zeros(Nb,length(time));
        xid = zeros(Nb,length(time));
        for ii = 1:Nb
            blade = int2str(10000+1000*ii + 30);
            % SIGN CRITERIA: In MBDyn the lead-lag is positive when leading,
            % while normally it is considered the other way around
            xi(ii,:) = -ncread(fn, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);
            xid(ii,:) = -ncread(fn, ['elem.joint.',blade,'.Omega'],[3, 1], [1, Inf]);
        end

        if find(omega_plot==index)
        nm = ['$\Omega$=' num2str(round(om*60/(2*pi),0))];
        fig(100+kk+tc) = figure(100+kk+tc);
        fig(100+kk+tc).Name = ['PS3D' damp{kk} suffix];

        ax(1) = subplot(3,2,1);
        plot3(ax(1),hub_x,v_x,time,'Color',Color{1},'LineStyle',LineStyle{index},...
            'Marker',Marker{index},'MarkerEdgeColor','k','MarkerIndices',1:nT:length(psi_nd),'DisplayName',nm)
        xlabel(ax(1),lblx{1})
        ylabel(ax(1),lbly{1})
        zlabel(ax(1),'time [s]')
        hold(ax(1),"on")
        legend(ax(1))

        ax(1) = subplot(3,2,2);
        plot3(ax(1),hub_y,v_y,time,'Color',Color{2},'LineStyle',LineStyle{index},...
            'Marker',Marker{index},'MarkerEdgeColor','k','MarkerIndices',1:nT:length(psi_nd),'DisplayName',nm)
        xlabel(ax(1),lblx{2})
        ylabel(ax(1),lbly{2})
        zlabel(ax(1),'time [s]')
        hold(ax(1),"on")
        legend(ax(1))

        for ii = 1 : Nb
            ax(1) = subplot(3,2,2+ii);
            plot3(ax(1),xi(ii,:),xid(ii,:),time,'Color',Color{2+ii},'LineStyle',LineStyle{index},...
                'Marker',Marker{index},'MarkerEdgeColor','k','MarkerIndices',1:nT:length(psi_nd),'DisplayName',nm)
            xlabel(ax(1),lblx{2+ii})
            ylabel(ax(1),lbly{2+ii})
            zlabel(ax(1),'revolutions')
            hold(ax(1),"on")
            legend(ax(1))
        end

        fig(fid) = figure(fid);
        fig(fid).Name = [fignames 'PhaseSpace' suffix];
        tlo(fid) = tiledlayout('flow');

        % x vs v_x
        ax(1) = nexttile(tlo(fid));
        plot(hub_x,v_x,'Color',Color{1},'Marker','o','MarkerEdgeColor','k','MarkerIndices',1:nT:length(psi_nd))
        hold(ax(1),"on")
        grid(ax(1),"on")
        xlabel(ax(1),'$x$ [m]')
        ylabel(ax(1),'$\dot{x}$ [m/s]')

        % y vs v_y
        ax(1) = nexttile(tlo(fid));
        plot(hub_y,v_y,'Color',Color{2},'Marker','o','MarkerEdgeColor','k','MarkerIndices',1:nT:length(psi_nd))
        hold(ax(1),"on")
        grid(ax(1),"on")
        xlabel(ax(1),'$y$ [m]')
        ylabel(ax(1),'$\dot{y}$ [m/s]')

        % x vs v_x
        for ii = 1 : Nb
            ax(1) = nexttile(tlo(fid));
            plot(xi(ii,:),xid(ii,:),'Color',Color{2+ii},'Marker','o','MarkerEdgeColor','k','MarkerIndices',1:nT:length(psi_nd))
            hold(ax(1),"on")
            grid(ax(1),"on")
            xlabel(ax(1),['$\xi_' int2str(ii) '$ [rad]'])
            ylabel(ax(1),['$\dot{\xi}_' int2str(ii) '$ [rad/s]'])
        end
        fid = fid+1;


        end

        % PSmat = [hub_x;v_x;hub_y;v_y;xi;xid];
        % MLCE(index,kk) = calculate_MLCE(PSmat,100,10,15,2,fs);
        index = index+1;
    end
end
disp('Lead-lag only done...')