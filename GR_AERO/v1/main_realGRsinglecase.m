% Author: Alejandro Alvaro, 2023-2024
% execute generic rigid rotor model and obtain the damping ratio after
% applying a perturbation in the longitudinal direction to the base.
%
% Objective: Under the ideal case scenario where only lead-lag motion is
% allowed and there are no collective inputs the damping ratios yield
% values of exactly one, justifying the validity of the EGDC proposed.
%
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
INPUT_THETA_0 = 10.*pi/180; % radian
INPUT_A_1_H = 0.*pi/180; % radian
INPUT_B_1_H = 5.*pi/180; % radian
% Write the prescribed Pitch angles
writematrix([INPUT_THETA_0;INPUT_A_1_H;INPUT_B_1_H].*180/pi,"Aeromechanics.xlsx",'Sheet','Hoja1','Range','B6:B8')

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
suffix = '';

if INPUT_THETA_0~=0
suffix = [suffix 'col' num2str(round(INPUT_THETA_0*180/pi,0))];
end
if INPUT_A_1_H~=0
suffix = [suffix 'cycA1' num2str(round(INPUT_A_1_H*180/pi,0))];
end
if INPUT_B_1_H~=0
suffix = [suffix 'cycB1' num2str(round(INPUT_B_1_H*180/pi,0))];
end

%% MBDyn variables
geom = calculateKgeom(opt,e,damp,Nb);

% Modify the geomtrical factors to elliminate the deutsch criterion and see
% % if the theroy holds
% deutschib = 1/(2*(1-cos(dpsi)));
% deutschi2b = 1/(2*(1-cos(2*dpsi)));
deutschib = 1;
deutschi2b = 1;

variables =['"' 'const real OMEGA_100 = ' num2str(omega) ...
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
    '; const real GEOM = ' num2str(geom.Kxidelta^2*deutschi2b) ...
    '; const real GEOMIB = ' num2str(geom.Kxil^2*deutschib) ...
    '; const real x_0 = ' num2str(x_0) ...
    '; const real y_0 = ' num2str(y_0) ...
    '"'];
setenv('MBDYNVARS', variables);
%% Figure initialization
if any(strcmp(damp, 'std'))
    labels = [];
    fig = figure(100);
    fig.Name = ['DampRatioComp' suffix];
    tlo(100) = tiledlayout('flow');
    ax(101) = nexttile(tlo(100));
    grid(ax(101),'on')
    hold(ax(101),'on')
    ylabel(ax(101),'$\zeta/\zeta^{std}$ \%')
    title(ax(101),'$\xi_{1c}$')
    ylim(ax(101),[0 120])
    ax(102) = nexttile(tlo(100));
    grid(ax(102),'on')
    hold(ax(102),'on')
    ylabel(ax(102),'$\zeta/\zeta^{std}$ \%')
    title(ax(102),'$\xi_{1s}$')
    ylim(ax(102),[0 120])
    if any(strcmp(damp, 'ib'))
        labels= [labels 'ib'];
    end
    if any(strcmp(damp, 'i2b'))
        labels= [labels 'i2b'];
    end

end

%% Execution of MBDyn and PP
% Table to verify the relationship between pitch and flap angles
aermechtab = zeros(3,3*length(damp));
% Ranges to write the respective values of the pitch and flap obtained form
% the MBDyn model
range = {'E6:E8','G6:G8','I6:I8'};
% MBC degrees of freedom
if mod(Nb,2)==0
        dof_col = [1 Nb];
else
         dof_col = 1;
end
dof_cyc = setdiff(1:Nb,dof_col);
% Allocate damping ratio
zeta_cyc = zeros(length(damp),length(dof_cyc));

for kk=1:length(damp)
    % Preprocessing of file names and testcase chosen
    if sw == 1
        swash = 'sw';
    else
        swash = 'nsw';
    end

    f_in = [pref_in folder model damp{kk} '_Nb' int2str(Nb) '_' swash cl '.mbd'];

    % Outfiles
    if b2hdamp==0
        suff = [damp{kk} 'ud'];
    else
        suff = [damp{kk} int2str(Nb) int2str(opt)];
    end

    fn_base = ['hhh' suff];

    f_out = [pref_in folder fn_base];

    % Preffix for figures
    fignames = [damp{kk} int2str(Nb)];
    
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
    psi_nd = time*omega/(2*pi);

    % Pp for fft and truncation of signals
    nT = 256;                           % Steps per period
    revN = fix(length(psi_nd)/nT) - 1;  % Number of revs -1


    % Extract lead-lag angles and azimut for each blade for MBC
    % transformation
    psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
    xi = zeros(Nb,length(psi_nd));
    xid = zeros(Nb,length(psi_nd));
    psi_b = zeros(Nb,length(psi_nd));

    for ii = 1:Nb
        psi_b(ii,:) = psi+ii*dpsi;
        blade = int2str(10000+1000*ii + 30);
        % SIGN CRITERIA: In MBDyn the lead-lag is positive when leading,
        % while normally it is considered the other way around
        xi(ii,:) = -ncread(fn, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);
        xid(ii,:) = -ncread(fn, ['elem.joint.',blade,'.Omega'],[3, 1], [1, Inf]);
    end
    % Store some indexes to perform the logarithmic decrement estimate
    % Perturbation time is from 40*T_REV to 40*T_REV+5*dt
    idx = 1:nT*12;
    T = 2*pi/omega;
    xi_NR_MB = MBC(psi_b,xi);
    qNR = NRdof(Nb);

    % Logarithmic decrement
    
    for jj = 1:length(dof_cyc)
    zeta_cyc(kk,jj) = logarithmic_decrement(psi_nd(idx),xi_NR_MB(dof_cyc(jj),idx),'PlotFlag',false,'MinPeakDistance',1);
    end


    if strcmp(damp{kk},'ib')
        ratib = zeta_cyc(kk,:)./zeta_cyc(1,:)*100;
        bib = bar(ax(101),["ib"], ratib(1));
        xtips1 = bib.XEndPoints;
        ytips1 = bib.YEndPoints;
        labels1 = string(bib.YData);
        text(ax(101),xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
        bib = bar(ax(102),["ib"],ratib(2));
        xtips1 = bib.XEndPoints;
        ytips1 = bib.YEndPoints;
        labels1 = string(bib.YData);
        text(ax(102),xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')

    elseif strcmp(damp{kk},'i2b')
        rati2b = zeta_cyc(kk,:)./zeta_cyc(1,:)*100;
        bi2b = bar(ax(101),["i2b"], rati2b(1));
        xtips1 = bi2b.XEndPoints;
        ytips1 = bi2b.YEndPoints;
        labels1 = string(bi2b.YData);
        text(ax(101),xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
        bi2b = bar(ax(102),["i2b"],rati2b(2));
        xtips1 = bi2b.XEndPoints;
        ytips1 = bi2b.YEndPoints;
        labels1 = string(bi2b.YData);
        text(ax(102),xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end

    % #9, Pith-flap lag
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames 'PFL' suffix];
    beta = -180/pi*ncread(fn, ['elem.joint.',blade,'.Phi'], [2, 1], [1, Inf]);
    Beta = fft(beta(nT*revN+[1:nT]));
    a0 = Beta/(nT);
    dt = time(2)-time(1);
    Fs = 1/dt;                          % Sampling frequency
    df = Fs/nT;                         % Increments of df for plotting
    [a,b] = fourierFFT(Beta',nT,df,omega,3);
    aermechtab(:,kk) = [a(1);a(2);b(2)];

    % Write the obtained flap angles and comapre them with the predicted
    % values from theory
    writematrix(aermechtab(:,kk),"Aeromechanics.xlsx",'Sheet','Hoja1','Range',range{kk})

    % theta = 180/pi*ncread(fn, ['elem.joint.',blade,'.Phi'], [1, 1], [1, Inf]);
    % Theta = fft(theta(nT*revN+[1:nT]));
    % theta0 = Theta(1)/nT;
    % 
    % [A,B] = fourierFFT(Theta',nT,df,omega,3);
    

    % inflow_label = '99';
    % lambda = ncread(fn, ['elem.inducedvelocity.',inflow_label,'.Lambda']);
    % Y = fft(lambda(nT*revN+[1:nT]));
    % lambda0 = Y(1)/(nT);
    % sigma = 0.08;
    % chord = 0.345575;
    % lock = 1.225*2*pi*chord*5.5^4/1084.7;
    % nub = sqrt(1+e*289.1/1084.7);
    % a0 = lock/(8*nub^2)*theta0;
    % % transform in deg, opposite sign for flap and lead-lag
    % plot(psi_nd, 180/pi*ncread(fn, ['elem.joint.',blade,'.Phi'], [1, 1], [1, Inf]), ...
    %     psi_nd, -180/pi*ncread(fn, ['elem.joint.',blade,'.Phi'], [2, 1], [1, Inf]), ...
    %     psi_nd, -180/pi*ncread(fn, ['elem.joint.',blade,'.Phi'], [3, 1], [1, Inf]));
    % xlabel('revolutions [-]');
    % xlim([psi_nd(1) psi_nd(end)])
    % ylabel(['Blade ' int2str(ii), ' angles, [deg]']);
    % legend('$\theta$', '$\beta$', '$\xi$');
    % grid on
    fid = kk*10+1;
    lastfid = fid;
end
%%
return
savefigures
