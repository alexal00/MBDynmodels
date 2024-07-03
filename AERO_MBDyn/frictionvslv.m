% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
%
% This function requires the script rotorbasicpp to use the basic post-processing tools
% generally included in the rotor_netcdf script. It is only used in the first comparisson
% that contains the basic plots to verify the correct execution of the MBDyn codes.
%
% This script is used to obtain the comparisson figures between the non-linear and linear
% damper constitutive law from Chapter 5.
close all; clear all; clc
%% Plotting parameters
setPlot
fid = 1;
lastfid = fid;
%% Helicopter data and enviroment variables
deg2rad = pi/180;
%================= FLIGHT CONDITIONS ============================
% % set MBDyn input for standard cruise flight
V_INF = 75.; % m/s 
ALTITUDE = 0.;% m, for standard air 
ALPHA_H = 2.41172*deg2rad; % radian 
INPUT_THETA_0 = 6.54609*deg2rad; % radian 
INPUT_A_1_H = 0.*deg2rad; % radian 
INPUT_B_1_H = 5.70475*deg2rad; % radian 

%======================= GENERAL HELICOPTER parameters=================
% Select number of blades
Nb = 4;         % [-], number of blades

omega = 40;     %[rad/s], rotor angular speed
nper = 70;      % [-], number of periods to run
R = 5.5;        % [m], rotor radius
e = 0.3048;     % [m], hinge position from hub
nuxi = sqrt(3/2*e/(R-e));

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

% Select constitutive Law
% Constitutive law
% cl = {'lve'};
cl = {'lve' 'fr'};
% cl = {'fr'};

% Select damping type
% damp = {'std' 'ib' 'i2b'};
% damp = {'ib' 'i2b'};
% damp = {'ib'};
damp = {'std' 'i2b'};
% damp = {'i2b'};

%% MBDyn filepaths
% Filepaths for MBDyn
% In preffix
pref_in = '/home/ale/thesis/AERO_MBDyn/';
folder = 'aero_NLV/';
model = 'rotor_';

% Out preffix
% pref_out = '/home/ale/thesis/AERO_MBDyn/aerodynamic_sw/';
% Shell file for MBDyn
shfile = '/home/ale/mbdyn/mbdyn/mbdyn ';

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


dpsi = 2*pi/Nb;
f_star = zeros(1,length(damp));
alpha_star = zeros(1,length(damp));
fd = 4067.5;
alpha = 10.;
Kxi = 0.; % Linear damping term, [Nm];
for kk = 1:length(damp)
   

    if strcmp(damp{kk},'std')
        K = 1;
        E = 1;
    elseif strcmp(damp{kk},'ib')
        K = geom.Kxil;
        E=(2*(1-cos(dpsi)));
        % alpha = alpha*K*E;
    elseif strcmp(damp{kk},'i2b')
        K = geom.Kxidelta;
        E=(2*(1-cos(2*dpsi)));
    end
    f_star(kk) = fd/(E*alpha*K);
    alpha_star(kk) = alpha/K;
end

variables =['"' 'const real R = ' num2str(R) ...
    '; const real X_BLADE_FLAP_OFFSET = ' num2str(e) ...
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
    '; const real fd =' num2str(fd) ...
    '; const real ALPHA = ' num2str(alpha) ...
    '; const real K_XI = ' num2str(Kxi) ...
    '"'];
setenv('MBDYNVARS', variables);
%% Execution of MBDyn and PP

idx = 1;
for ll=1:length(cl)
    for kk=1:length(damp)

        if sw ==1
            swash = 'sw';
            
        else
            swash = 'nsw';
        end

        modelname =  [model damp{kk} '_Nb' int2str(Nb) '_' swash  cl{ll} '.mbd'];
        disp(['Input file: ' modelname])
        f_in = [pref_in folder modelname];
        
        fn_base = ['hhh' damp{kk} int2str(Nb) int2str(opt) swash cl{ll}];
        disp(['Out-file: ' fn_base])
        f_out{idx} = [pref_in folder fn_base];
    
        if exec==1
            disp('executing MBDyn...');
        
        
            [rc, errmsg] = system(['wsl ' 'MBDYNVARS=' variables ' ' shfile  f_in ' -o ' f_out{idx}]);
        
            % [rc, errmsg] = system(['./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
            if (rc ~= 0)
                error(errmsg);
            end
            disp('   ... done');
        end
    idx = idx +1;
    end

end
%% Linear viscoelastic pp
pplve = 1.;
clplt = 1;
if pplve ==1
    for ii = 1:length(cl)
    aux = 1;
    filenc = {};
        for jj = 1:length(f_out)
            % Find indexes of damp{ii} in f_out{jj}
            match_idx = strfind(f_out{jj}, cl{ii});
            if ~isempty(match_idx)
                % fidx = [fidx, jj];
                filenc{aux} = f_out{jj};
                aux = aux+1;
            end
        end
        cli = cl{ii};
        rotorbasicpp
        clplt = clplt+1;
    end
    
end
disp('End of initial postprocessing regarding harmonic content')
disp('Continue with comparisson between Constitutive Laws? PRESS ENTER')
pause
%% Comparisson between linear and non-linear
% Loop through each string in damp
if sw == 0
    rn1 = {'B5:B14' 'D5:D14'};
    rn2 = {'C5:C14' 'E5:E14'};
else
    rn1 = {'F5:F14' 'H5:H14'};
    rn2 = {'G5:G14' 'I5:I14'};
end

for ii = 1:length(damp)
    % Initialize an empty array to store indexes
    fidx = [];
    files = {};
    aux = 1;
    % Loop through each string in f_out
    for jj = 1:length(f_out)
        % Find indexes of damp{ii} in f_out{jj}
        match_idx = strfind(f_out{jj}, damp{ii});
        % If a match is found, store the index
        if ~isempty(match_idx)
            % fidx = [fidx, jj];
            files{aux} = f_out{jj};
            aux = aux+1;
        end
    end
    % % Figure 3: M_d vs theta and M_d vs theta_prime
    fhandle = figure(fid); fid = fid+1;
    fhandle.Name = [damp{ii} 'constituetivelaw' dof];

    ax(1) = subplot(1,2,1);
    grid(ax(1),"on")
    hold(ax(1),"on")
    legend(ax(1))

    ax(2) = subplot(1,2,2);
    xlabel(ax(2),'$\dot{\phi}$ [rad/s]')
    ylabel(ax(2),'$M_d$ [Nm]')
    grid(ax(2),"on")
    hold(ax(2),"on")
    legend(ax(2))
    
     if strcmp(damp{ii},'ib')
        xlabel(ax(1),'$l$ [m]')
        ylabel(ax(1),'$F_d$ [N]')
        xlabel(ax(2),'$\dot{l}$ [m/s]')
        ylabel(ax(2),'$F_d$ [N]')
    else
        xlabel(ax(1),'$\phi$ [rad]')
        ylabel(ax(1),'$M_d$ [Nm]')
        xlabel(ax(2),'$\dot{\theta}$ [rad/s]')
        ylabel(ax(2),'$M_d$ [Nm]')
     end


    fhandle = figure(fid); fid = fid+1;
    fhandle.Name = [damp{ii} 'HarmonicPP' swash dof 'fr'];
    tlo(1) = tiledlayout(3,2);
    ax(100) = nexttile(tlo(1),1,[1 2]);
    xlim(ax(100),[0,9])
    xticks(ax(100),0:9)
    xlabel(ax(100),'harmonic')
    ylabel(ax(100),'$\mathrm{PSD}(M_{\xi})$ [N or Nm]')
    grid(ax(100),"on")
    hold(ax(100),"on")
    legend(ax(100));
    ax(110) = nexttile(tlo(1),3,[1 2]);
    xlim(ax(110),[0,9])
    xticks(ax(110),0:9)
    xlabel(ax(110),'harmonic')
    ylabel(ax(110),'$\mathrm{PSD}(\dot{\Delta})$ [rad/s or m/s]')
    grid(ax(110),"on")
    hold(ax(110),"on")
    legend(ax(110));
    ax(120) = nexttile(tlo(1),5,[1 2]);
    hold(ax(120),"on")
    xlim(ax(120),[0,9])
    xticks(ax(120),0:9)
    xlabel(ax(120),'harmonic')
    ylabel(ax(120),'$E_d^n/E_d$ [\%]')
    grid(ax(120),"on")
    legend(ax(120));

    harmonic_energy = zeros(length(files),11);
    contribution_percent = zeros(length(files),11);
    for jj = 1:length(files)
        outpath = [files{jj} '.nc'];
        disp(sprintf('reading output from file ''%s''', outpath));

        def_hinge = ['elem.joint.' int2str(10000+Nb*1000+30+1) '.'];
        if strcmp(damp{ii},'ib')
            theta = ncread(outpath, [def_hinge 'l']);     % Rod length, [m]
            thetap = ncread(outpath, [def_hinge 'lP']);     % Rod lengthening velocity, [m/s]
            Md = ncread(outpath,[def_hinge 'f'],[1, 1], [1, Inf])'; % Reation force*geom.Kxil, [N,]
        else
            theta = ncread(outpath, [def_hinge 'E'],[3, 1], [1, Inf])';   % Relative orientation, [rad]
            thetap = ncread(outpath, [def_hinge 'Omega'],[3, 1], [1, Inf])';    % Angular velocity, [rad/s]
            if strcmp(damp{ii},'i2b')                                     
                fact =-1;
            else
                fact = 1;
            end
            Md = fact*ncread(outpath,[def_hinge 'M'],[3, 1], [1, Inf])';    % Reaction moment, [Nm]
        end
        

        functionplot(ax(1),theta,Md,1,Color{jj},Marker{jj},'none',cl{jj});

        functionplot(ax(2),thetap,Md,1,Color{jj},Marker{jj},'none',cl{jj});

        % compute azimuth vector from database
        time = ncread(outpath, 'time');
        psi_nd = time*omega/(2*pi);
        % Rescale to zero time vector
        psi_nd = psi_nd-psi_nd(1);        

        blade =  ['elem.joint.' int2str(10000+Nb*1000+30) '.'];

        % Compute the energy dissipated by the damper in one oscillation
        % cycle
        T = 2*pi/omega;
        nT = 256;
        revN = fix(length(psi_nd)/nT) - 1;
   
        dt = time(2)-time(1);               % time step
        Fs = 1/dt;                          % Sampling frequency
        df = Fs/nT;                         % Increments of df for plotting
        index = revN*nT+[1:nT];             
        % Truncate signal to one period
        t_T = time(index);
        % Compute one-sided PSD of the moments
        [P1Md,f] = computeFFT(Md(index),t_T);
        f_nd = f/(1/T);                     % Frequency vector
        % Plot the norm of the n-th harmonic
        stem(ax(100),f_nd,P1Md,'filled','DisplayName',cl{jj})
        % fft of the angular velocity in the damper
        fftxid = fft(thetap(index));
        [P1xid,f] = computeFFT(thetap(index),t_T);
        stem(ax(110),f_nd,P1xid,'filled','DisplayName',cl{jj})
        
        % Compute harmonic content
        nF = 10;            % Number of harmonics to obtain
        [Ed, harmonic_energy(jj,:), contribution_percent(jj,:)] = calculate_energy_dissipation(Md(index),thetap(index),t_T,10);

        % writematrix(harmonic_energy(jj,1:10)',"SpectraAermech.xlsx",'Sheet','Hoja2','Range',rn1{ii})
        % writematrix(contribution_percent(jj,1:10)',"SpectraAermech.xlsx",'Sheet','Hoja2','Range',rn2{ii})
    end

    bar(ax(120),0:nF,contribution_percent')
    legend(ax(120),cl);
    
    fhandle = figure(fid); 
    fhandle.Name = [damp{ii} 'EnergyHarmonic' dof 'fr'];
    tlo(fid) = tiledlayout('flow');
    ax(fid) = nexttile(tlo(fid));
    xlim(ax(fid),[0,nF])
    xticks(ax(fid),0:nF)
    xlabel(ax(fid),'harmonic')
    ylabel(ax(fid),'$E_d$ [J]')
    grid(ax(fid),"on")
    hold(ax(fid),"on")
    bar(ax(fid),0:nF,harmonic_energy')
    legend(ax(fid),cl)
    ax(fid) = nexttile(tlo(fid));
    xlim(ax(fid),[0,nF])
    xticks(ax(fid),0:nF)
    xlabel(ax(fid),'harmonic')
    ylabel(ax(fid),'$E_d^n/E_d$ [\%]')
    grid(ax(fid),"on")
    hold(ax(fid),"on")
    bar(ax(fid),0:nF,contribution_percent')
    legend(ax(fid),cl)
    fid = fid+1;
    

end
%%
return
savefigures
