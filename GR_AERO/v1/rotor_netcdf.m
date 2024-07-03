% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
%
close all; clear all; clc
%% Plotting parameters
setPlot
fid = 1;
lastfid = fid;
%% Helicopter data and enviroment variables
deg2rad = pi/180;
%================= FLIGHT CONDITIONS ============================
% set MBDyn input for 50 m/s forward flight
% V_INF = 50.; % m/s
% ALTITUDE = 0.; % m
% % NOTE: for convenience, mount the rotor with vertical shaft and put gamma + tau in alpha_H
% ALPHA_H = 0.*pi/180; % radian
% INPUT_THETA_0 = 0*pi/180; % radian
% INPUT_A_1_H = 0.*pi/180; % radian
% INPUT_B_1_H = 0.*pi/180; % radian
% % set MBDyn input for hover
V_INF = 0.; % m/s
ALTITUDE = 0.; % m
% NOTE: for convenience, mount the rotor with vertical shaft and put gamma + tau in alpha_H
ALPHA_H = 0.*pi/180; % radian
INPUT_THETA_0 = 7.66777*pi/180; % radian
INPUT_A_1_H = 0*pi/180; % radian
INPUT_B_1_H = 0*pi/180; % radian

% % set MBDyn input for standard cruise flight
% V_INF = 75.; % m/s 
% ALTITUDE = 0.;% m, for standard air 
% ALPHA_H = 2.41172*deg2rad; % radian 
% INPUT_THETA_0 = 6.54609*deg2rad; % radian 
% INPUT_A_1_H = 0.*deg2rad; % radian 
% INPUT_B_1_H = 5.70475*deg2rad; % radian 

%======================= GENERAL HELICOPTER parameters=================
% Select number of blades
Nb = 4;         % [-], number of blades

omega = 40;     %[rad/s], rotor angular speed
nper = 60;      % [-], number of periods to run
R = 5.5;        % [m], rotor radius
e = 1/25*R;     % [m], hinge position from hub
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
exec = 1;

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

variables =['"' 'const real R = ' num2str(R) ...
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
    '"'];
setenv('MBDYNVARS', variables);
%% Figure initialization

% Figure for harmonic comparisson between models
% fig(100) = figure(100);
% fig(100).Name = 'PSDDeltaxidotcompls';
% tlo(100) = tiledlayout('flow');
% for ii = 1: Nb
%     ax(100+ii) = nexttile(tlo(100));
%     yscale(ax((100+ii)),"log")
%     hold(ax(100+ii),"on")
%     grid(ax(100+ii),"on")
% end

% Figure for harmonic comparisson between models, without log-scale
fig(110) = figure(110);
fig(110).Name = 'PSDDeltaxidotcomp';
tlo(110) = tiledlayout('flow');
for ii = 1: Nb
    ax(110+ii) = nexttile(tlo(110));
    hold(ax(110+ii),"on")
    grid(ax(110+ii),"on")
    xlim(ax(110+ii),[0,9])
    xticks(ax(110+ii),0:9)
    xlabel(ax(110+ii),'harmonic')
    ylabel(ax(110+ii),['$PSD(\Delta\dot{\xi}_' num2str(ii) ')$'])
    grid(ax(110+ii),"on")
end

% Figure for harmonic comparisson between models, without log-scale
fig(210) = figure(210);
fig(210).Name = ['PSDmxicomp' 'Nb' int2str(Nb) '_' dof];
tlo(210) = tiledlayout('flow');
for ii = 1: Nb
    ax(210+ii) = nexttile(tlo(210));
    hold(ax(210+ii),"on")
    grid(ax(210+ii),"on")
    xlim(ax(210+ii),[0,9])
    xticks(ax(210+ii),0:9)
    xlabel(ax(210+ii),'harmonic')
    ylabel(ax(210+ii),['$PSD(m_{\xi_' num2str(ii) '})$'])
    grid(ax(210+ii),"on")
end

if any(strcmp(damp, 'std'))
    if any(strcmp(damp, 'i2b'))
        fig(300) = figure(300);
        fig(300).Name = ['rati2bNb' int2str(Nb) 'opt' int2str(opt) '_' dof];
        ax(310) = subplot(1,2,1);
        hold(ax(310),"on")
        grid(ax(310),"on")
        xlabel(ax(310),'harmonic')
        ylabel(ax(310),'$\Delta\dot{\xi}_{i2b}/\Delta\dot{\xi}_{std}$')
        xlim(ax(310),[0 10])
        xticks(ax(310),1:9)
        ax(320) = subplot(1,2,2);
        hold(ax(320),"on")
        grid(ax(320),"on")
        xlabel(ax(320),'harmonic')
        ylabel(ax(320),'$m_{\xi}^{i2b}/m_{\xi}^{std}$')
        % xticklabels(ax(320),1 2 3 4 5 6 7 8 9)
        xlim(ax(320),[0 10])
        xticks(ax(320),1:9)
    end

    if any(strcmp(damp, 'ib'))
        fig(350) = figure(350);
        fig(350).Name = ['ratibNb' int2str(Nb) 'opt' int2str(opt) '_' dof];
        ax(360) = subplot(1,2,1);
        hold(ax(360),"on")
        grid(ax(360),"on")
        xlabel(ax(360),'harmonic')
        ylabel(ax(360),'$\Delta\dot{\xi}_{ib}/\Delta\dot{\xi}_{std}$')
        xlim(ax(360),[0 10])
        xticks(ax(360),1:9)
        ax(370) = subplot(1,2,2);
        hold(ax(370),"on")
        grid(ax(370),"on")
        xlabel(ax(370),'harmonic')
        ylabel(ax(370),'$m_{\xi}^{ib}/m_{\xi}^{std}$')
        xlim(ax(370),[0 10])
        xticks(ax(370),1:9)
    end
end
%% Execution of MBDyn and PP
% Correction factor for the same equivalent damping
c_xi = 4067.5;
c_xic = zeros(1,length(damp));
dpsi = 2*pi/Nb;
for kk = 1:length(damp)
    if strcmp(damp{kk},'std')
        fact = 1;
    elseif strcmp(damp{kk},'ib')
        Kxil = geom.Kxil;

        fact=1/(2*(1-cos(dpsi)))/(Kxil^2);
        c_d = 1.5567e+05;
    elseif strcmp(damp{kk},'i2b')
        Kxidelta = geom.Kxidelta;
        fact=1/(2*(1-cos(2*dpsi)))/Kxidelta^2;
    end
    c_xic(kk) = c_xi*fact;
end


for kk=1:length(damp)
    if sw ==1
        f_in = [pref_in folder model damp{kk} '_Nb' int2str(Nb) '_sw' cl '_ODI.mbd'];
    else
        f_in = [pref_in folder model damp{kk} '_Nb' int2str(Nb) '_nsw' cl '_ODI.mbd'];
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
    
    % angular velocity
    %%% Omega = 40; % rad/s
    % we could get it from the motion of the hub, assuming it starts with nominal RPM
    Omega = ncread(fn, 'node.struct.10000.Omega', [3, 1], [1, 1]);
    
    % compute azimuth vector from database
    time = ncread(fn, 'time');
    psi_nd = time*Omega/(2*pi);
    
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
    deltapsi = 2*pi/nb;
    fignames = [damp{kk} int2str(nb)];
    blade_label = int2str(10000 + 1000*nb);
    joint_label = int2str(10000 + 1000*nb + 30);
    inflow_label = '99';
    
    disp(sprintf('  found %d blades; using blade %d, label=''%s''', nb, nb, blade_label));
    
    % extract number of Gauss points from database
    label = blade_label;
    blade_txt = ['blade ', int2str(nb)];
    ng = 0;
    while (1)
        n = 1;
        while (n <= size(finfo.Variables, 2))
            if (strcmp(finfo.Variables(n).Name, ['elem.aerodynamic.',label,'.alpha_',int2str(ng)]))
                break;
            end
            n = n + 1;
        end
        if (n > size(finfo.Variables, 2))
            break;
        end
        ng = ng + 1;
    end
    
    disp(sprintf('  found %d Gauss points', ng));
    
    nT = 256;
    revN = fix(length(psi_nd)/nT) - 1;

    T_xi = 2*pi/(nuxi*Omega);

    revNxi = (time(end)/(T_xi));
    nTxi = fix(length(psi_nd)/revNxi);
    revNxi = fix(length(psi_nd)/nTxi)-2;

    % % #1, Lead-lag of single blade
    % % Figure = kk*10+1
    % % Tiled layout = kk*10+1
    fig(lastfid) = figure(lastfid);
    fig(lastfid).Name = [fignames 'llsb'];
    tlo(lastfid) = tiledlayout('flow');
    fid = fid +1 ;
    % 
    % % #2, Lead-lag of all blades
    % % Figure = kk*10+1+1
    fig(lastfid+1) = figure(lastfid+1);
    fig(lastfid+1).Name = [fignames 'llallb'];
    ax(2) = gca;
    hold on;
    fid = fid +1 ;

    psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
    for ii = 1:nb
        psi_b(ii,:) = psi+ii*deltapsi;
        blade = int2str(10000+1000*ii + 30);
        xi(ii,:) = ncread(fn, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);
        xid(ii,:) = ncread(fn, ['elem.joint.',blade,'.Omega'],[3, 1], [1, Inf]);

        ax(1) = nexttile(tlo(lastfid));
        hold(ax(1),"on");
        functionplot(ax(1),psi_nd,-180/pi*xi(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\xi_' int2str(ii) '$'])
        grid(ax(1),"on")
        xlabel(ax(1),'revolutions [-]');
        ylabel(ax(1),'$\xi$ [deg]')

        functionplot(ax(2),psi_nd,-180/pi*xi(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\xi_' int2str(ii) '$'])
        legend(ax(2))
        grid(ax(2),"on")
        xlabel(ax(2),'revolutions [-]')
        ylabel(ax(2),'$\xi$ [deg]')
    end

    
    % % # 3, MBC
    % % Figure = kk*10++1+2
    xi_NR_MB = MBC(psi_b,xi);
    fig(lastfid+2)=figure(lastfid+2); 
    fig(lastfid+2).Name = [fignames 'qNR'];
    % fid = fid+1;
    qNR = NRdof(nb);
    ax(3)= subplot(1,2,1);
    % Plot collective DOF
    if mod(nb,2)==0
        dof_col = [1 nb];
        functionplot(ax(3),psi_nd,xi_NR_MB(1,:),step,Color{1},Marker{1},LineStyle{1},[qNR{1} ', MBDyn'])
        % plot(psi_nd, xi_NR_a(1,:),'DisplayName',qNR{1},'Color',Color{1}); hold on
        functionplot(ax(3),psi_nd,xi_NR_MB(end,:),step,Color{2},Marker{2},LineStyle{2},[qNR{end} ', MBDyn'])
        % plot(psi_nd, xi_NR_a(end,:),'DisplayName',qNR{end},'Color',Color{2}); hold on
    else
        dof_col = 1;
        functionplot(ax(3),psi_nd,xi_NR_MB(1,:),step,Color{1},Marker{1},LineStyle{1},[qNR{1} ', MBDyn'])
    end
    xlabel('revolutions [-]');
    ylabel('Collective $q_{NR}$ [rad]')
    legend('show','interpreter','latex')
    grid on

    dof_cyc = setdiff(1:nb,dof_col);
    ax(4)= subplot(1,2,2);
    for jj = 1:length(dof_cyc)
        functionplot(ax(4),psi_nd,xi_NR_MB(dof_cyc(jj),:),step,Color{jj},Marker{jj},LineStyle{jj},[qNR{dof_cyc(jj)} ', MBDyn'])
    end
    xlabel('revolutions [-]');
    ylabel('Cyclic $q_{NR}$ [rad]')
    legend
    grid on
    % # 4 hub motions
    AIRFRAME_X = 100; 
    AIRFRAME_Y = 200;
    hub_x = ncread(fn, ['node.struct.',int2str(AIRFRAME_X),'.X'],[1, 1], [1, Inf]);
    hub_y = ncread(fn, ['node.struct.',int2str(AIRFRAME_Y),'.X'],[2, 1], [1, Inf]);
    v_x = ncread(fn,['node.struct.' num2str(AIRFRAME_X) '.XP'],[1,1],[1,Inf]);
    v_y = ncread(fn,['node.struct.' num2str(AIRFRAME_Y) '.XP'],[2,1],[1,Inf]);
    fig(lastfid+3)= figure(lastfid+3);
    fig(lastfid+3).Name = [fignames 'Hubxy'];
    ax(1) = subplot(2,1,1);
    functionplot(ax(1),psi_nd,hub_x,step,Color{1},Marker{1},LineStyle{1}, '$x_{hub}$'); hold on
    functionplot(ax(1),psi_nd,hub_y,step,Color{2},Marker{2},LineStyle{2}, '$y_{hub}$')
    xlabel('revolutions [-]');
    ylabel('$x,y$ [m]')
    legend

    ax(2) = subplot(2,1,2);
    functionplot(ax(2),psi_nd,v_x,step,Color{1},Marker{1},LineStyle{1}, '$v_{x,hub}$'); hold on
    functionplot(ax(2),psi_nd,v_y,step,Color{2},Marker{2},LineStyle{2}, '$v_{y,hub}$')
    xlabel('revolutions [-]');
    ylabel('$v_x,v_y$ [m/s]')
    legend

    % % #4, Relative velocity at dampers
    % % Figure = kk*10+1+3
    % % Tiled layout = kk*10+1+3
    % fig(lastfid+3) = figure(lastfid+3); 
    % fig(lastfid+3).Name = [fignames 'DeltaXidot'];
    % tlo(lastfid+3) = tiledlayout('flow');
    % fid=fid+1;

    for ii = 1:nb
        xi_d = xid(ii,:) ;
        if strcmp(damp{kk},'std')
            nxt = 0;
            int_b = ii;                 % Intermediate blade
        elseif strcmp(damp{kk},'ib')
            nxt = mod(ii, Nb) + 1;      % Calculate next blade
            int_b = ii ;
        elseif strcmp(damp{kk},'i2b')
            nxt = mod(ii+1, Nb) + 1;    % Calculate next blade
            int_b = mod(ii, Nb) + 1;    % Intermediate blade
        end

        if nxt ==0
            deltaxidot(ii,:) = xi_d;
        else
            xi_dnxt = xid(nxt,:);
            deltaxidot(ii,:) = xi_dnxt-xi_d;
        end
        steps = nT;
        xid_m = movmean(deltaxidot(ii,:),steps);
        calculate_rms = @(funct) sqrt(mean(funct.^2));
        
        % Calculate RMS for each time step
        rms_values(ii,:) = arrayfun(@(n) calculate_rms(xid_m(1:n)), ...
            1:length(deltaxidot));
        defhing = int2str(10000+1000*int_b + 31);
        if strcmp(damp{kk},'i2b')
            k = -1;
            kappa = geom.Kxidelta;
        elseif strcmp(damp{kk},'ib')
            k = 1;
            kappa = geom.Kxil;
        else
            k= 1;
            kappa = 1;
        end
        if ~strcmp(damp{kk},'ib')
        dxidot_mbd(ii,:) = 1/kappa*k*ncread(fn, ['elem.joint.',defhing,'.Omega'],[3, 1], [1, Inf]);
        else
        dxidot_mbd(ii,:) = 1/kappa*k*ncread(fn, ['elem.joint.',defhing,'.lP']);
        end

        
        % ax(1) = nexttile(tlo(lastfid+3));
        % functionplot(ax(1),psi_nd,deltaxidot(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\Delta\dot{\xi}_' int2str(ii) '$'])
        % functionplot(ax(1),psi_nd,dxidot_mbd(ii,:),step,'k','none',':',['$\Delta\dot{\delta}_' int2str(ii) '/K_{\xi\delta}$, MBDyn'])
        % functionplot(ax(1),psi_nd,rms_values(ii,:),step,'r','none','-',['$\mathrm{RMS}(\dot{\xi}_' int2str(ii) ')$'])
        % grid(ax(1),"on")
        % legend(ax(1))
        % xlabel(ax(1),'revolutions [-]');
        % ylabel(ax(1),['$\Delta\dot{\xi}_' int2str(ii) '$ [rad/s]']) 
    end
    time_o = time(nT*revN+[1:nT]);
    [P1xi,f_o]=computeFFT(deltaxidot(:,nT*revN+[1:nT])',time_o);
    f_ond = f_o./(Omega/(2*pi));

    time_xi = time(nTxi*revNxi+[1:nTxi*2]);
    [P1xi_xi,f_xi]=computeFFT(deltaxidot(:,nTxi*revNxi+[1:nTxi*2])',time_xi);
    f_xind = f_xi./(nuxi*Omega/(2*pi));
    

    % Moments in the damper
    % #6
    % Figure = kk*10+1+5
    % Tiled layout = kk*10+1+5
    fig(lastfid+5)=figure(lastfid+5); 
    fig(lastfid+5).Name = [fignames 'Mxi'];
    tlo(lastfid+5) = tiledlayout('flow');
    fid = fid+1;

    % #7
    % Figure = kk*10+1+6
    % Tiled layout = kk*10+1+6
    % fig(lastfid+6)=figure(lastfid+6); 
    % fig(lastfid+6).Name = [fignames 'Mxicomp'];
    % tlo(lastfid+6) = tiledlayout('flow');
    % fid = fid+1;

    % #8
    % Figure = kk*10+1+7
    % Tiled layout = kk*10+1+7
    % fig(lastfid+7)=figure(lastfid+7); 
    % fig(lastfid+7).Name = [fignames 'PSDMxi'];
    % tlo(lastfid+7) = tiledlayout('flow');
    % fid = fid+1;
    for ii = 1:nb
        if ~strcmp(damp{kk},'i2b')
        nxt = ii;  % Calculate next blade
        elseif strcmp(damp{kk},'i2b')
        nxt = mod(ii, Nb) + 1;  % Calculate the intermediate blade
        end
        defhing = int2str(10000+1000*nxt + 31);
        if strcmp(damp{kk},'std')
            m_xi(ii,:) = ncread(fn, ['elem.joint.',defhing,'.m'],[3, 1], [1, Inf]);
            kappa = 1;
        elseif strcmp(damp{kk},'i2b')
            m_xi(ii,:) = -1*ncread(fn, ['elem.joint.',defhing,'.m'],[3, 1], [1, Inf]);
            kappa = geom.Kxidelta;
        else
            Force_xi(ii,:) = ncread(fn, ['elem.joint.',defhing,'.f'],[1, 1], [1, Inf]);
            m_xi(ii,:) = Force_xi(ii,:);%*geom;
            dV = ncread(fn, ['elem.joint.',defhing,'.lP']);
            kappa = geom.Kxil;
        end
        ax(2) = nexttile(tlo(lastfid+5));
        functionplot(ax(2),psi_nd,m_xi(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$m_{\xi,' int2str(ii) '}$'])
        grid(ax(2),"on")
        legend(ax(2))
        xlabel(ax(2),'revolutions [-]');
        ylabel(ax(2),['$m_{\xi,' int2str(ii) '}$ [Nm]'])
    
        % ax(3) = nexttile(tlo(lastfid+6));
        % functionplot(ax(3),psi_nd,c_xic(kk).*kappa*deltaxidot(ii,:),step,Color{ii},Marker{ii},'-','$C_{xi}\Delta\dot{\xi}$')
        % functionplot(ax(3),psi_nd,m_xi(ii,:),step,'k','none','-.','$m_{\xi}$, MBDyn')
        % grid(ax(3),"on")
        % legend(ax(3))
        % xlabel(ax(3),'revolutions [-]');
        % ylabel(ax(3),['$m_{\xi,' int2str(ii) '}$'])


        [P1m,f_o]=computeFFT(m_xi(:,nT*revN+[1:nT])',time_o);

        % ax(3) = nexttile(tlo(lastfid+7));
        % stem(ax(3),f_ond,P1m(:,ii));
        % xlim(ax(3),[0 10])
        % xlabel(ax(3),'f [Hz]')
        % ylabel(ax(3),'$\mathrm{PSD}(m_{\xi})$') 
    end

    
    for ii=1:nb
    %     % Plot one-sided espectrum of lead-lag angular speed
    %     % stem(ax(100+ii),f_ond',P1xi(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
    %     % xlabel(ax(100+ii),'harmonic')
    %     % ylabel(ax(100+ii),['$PSD(\Delta\dot{\xi}_' num2str(ii) ')$'])
    %     % xlim(ax(100+ii),[0 10])
    %     % % title(ax(100+ii),['Revolution=' num2str(revN)])
    %     % legend(ax(100+ii))
    %     % grid(ax(100+ii),"on")
    % 
    %     % Plot one-sided espectrum of lead-lag angular speed without log
    %     % scale
        stem(ax(110+ii),f_ond',P1xi(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
        % xlim(ax(110+ii),[0 10])
        % title(ax(110+ii),['Revolution=' num2str(revN)])
        legend(ax(110+ii))
        
    % 
    %     % Plot one-sided espectrum of lead-lag angular speed without log
    %     % scale
    %     % stem(ax(180+ii),f_xind',P1xi_xi(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
    %     % xlabel(ax(180+ii),'harmonics in $\xi$')
    %     % ylabel(ax(180+ii),['$PSD(\Delta\dot{\xi}_' num2str(ii) ')$'])
    %     % xlim(ax(180+ii),[0 30])
    %     % % title(ax(170+ii),['Revolution=' num2str(revN)])
    %     % legend(ax(180+ii))
    %     % grid(ax(180+ii),"on")
    % 
    %     % Plot one-sided espectrum of moments in damper
    %     % stem(ax(200+ii),f_ond',P1m(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
    %     % xlabel(ax(200+ii),'harmonic')
    %     % ylabel(ax(200+ii),['$PSD(m_{\xi_' num2str(ii) '})$'])
    %     % xlim(ax(200+ii),[0 10])
    %     % % title(ax(100+ii),['Revolution=' num2str(revN)])
    %     % legend(ax(200+ii))
    %     % grid(ax(200+ii),"on")
    % 
    %     % Plot one-sided espectrum of moments in damper without log
    %     % scale
        if strcmp(damp{kk},'ib')
            stem(ax(210+ii),f_ond',geom.Kxil*P1m(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
        else
            stem(ax(210+ii),f_ond',P1m(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
        end

        legend(ax(210+ii))
       
    end

    % Plot the evolution of the harmonic content with time
    % stp = floor(revN/5);
    
    % for ii=1:nb
    %     index = 1;
    %     for rev=0:stp:revN
    %         c = mod(index,length(Color))+1;
    %         [P1xin,f_o]=computeFFT(deltaxidot(:,nT*rev+[1:nT])',time_o);
    %         % Plot one-sided espectrum of lead-lag angular speed
    %         stem(ax(idx(kk)+ii),f_ond',P1xin(:,ii),'filled',Color=Color{c}, ...
    %             Marker=Marker{c},DisplayName=['n=' num2str(rev)]); hold on
    %         xlabel(ax(idx(kk)+ii),'harmonic')
    %         ylabel(ax(idx(kk)+ii),['$PSD(\Delta\dot{\xi}_' num2str(ii) ')$'])
    %         xlim(ax(idx(kk)+ii),[0 10])
    %         legend(ax(idx(kk)+ii))
    %         grid(ax(idx(kk)+ii),"on")
    % 
    %         % c = mod(index,length(Color))+1;
    %         % [P1xin,f]=computeFFT(deltaxidot(:,nT*rev+[1:nT])',time);
    %         % Plot one-sided espectrum of lead-lag angular speed
    %         stem(ax(idx2(kk)+ii),f_ond',P1xin(:,ii),'filled',Color=Color{c}, ...
    %             Marker=Marker{c},DisplayName=['n=' num2str(rev)]); hold on
    %         xlabel(ax(idx2(kk)+ii),'harmonic')
    %         ylabel(ax(idx2(kk)+ii),['$PSD(\Delta\dot{\xi}_' num2str(ii) ')$'])
    %         xlim(ax(idx2(kk)+ii),[0 10])
    %         legend(ax(idx2(kk)+ii))
    %         grid(ax(idx2(kk)+ii),"on")
    % 
    %         index = index+1;
    %     end
    % end
    % #9, Pith-flap lag
    % label = joint_label;
    % fig(lastfid+8)=figure(lastfid+8); 
    % fig(lastfid+8).Name = [fignames 'PFL'];
    % % fid = fid+1;
    % % transform in deg, opposite sign for flap and lead-lag
    % plot(psi_nd, 180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [1, 1], [1, Inf]), ...
    %     psi_nd, -180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [2, 1], [1, Inf]), ...
    %     psi_nd, -180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [3, 1], [1, Inf]));
    % xlabel('revolutions [-]');
    % ylabel([blade_txt, ' angles, [deg]']);
    % legend('$\theta$', '$\beta$', '$\xi$');
    % grid on

    if any(strcmp(damp,'std'))
        x = 1:9;
        if strcmp(damp{kk},'std')
            PSDxid_std = P1xi(2:10,1);
            PSDMxi_std = P1m(2:10,1);
        elseif strcmp(damp{kk},'ib')
            ratdxi =  P1xi(2:10,1)./PSDxid_std;
            % ratmxi =  P1m(2:10,1)./PSDMxi_std;
            ratmxi =  geom.Kxil*P1m(2:10,1)./PSDMxi_std;
    
            % Analytical values
            ratdxi_a = sqrt(2*(1-cos(x*dpsi)));
            ratmxi_a = sqrt(2*(1-cos(x*dpsi)))/(2*(1-cos(dpsi)));
            
            ratio_dxi = [ratdxi';ratdxi_a];
            ratio_m = [ratmxi';ratmxi_a];
            bar(ax(360),x,ratio_dxi)
            legend(ax(360),{'MBDyn' 'Analytical'})
            bar(ax(370),x,ratio_m)
            legend(ax(370),{'MBDyn' 'Analytical'})
        elseif strcmp(damp{kk},'i2b')
            ratdxi =  P1xi(2:10,1)./PSDxid_std;
            ratmxi =  P1m(2:10,1)./PSDMxi_std;
            
            % Analytical values
            ratdxi_a = sqrt(2*(1-cos(x*2*dpsi)));
            ratmxi_a = sqrt(2*(1-cos(x*2*dpsi)))/(2*(1-cos(2*dpsi)))/geom.Kxidelta;
            
            ratio_dxi = [ratdxi';ratdxi_a];
            ratio_m = [ratmxi';ratmxi_a];
    
            bar(ax(310),x,ratio_dxi)
            legend(ax(310),{'MBDyn' 'Analytical'})
            bar(ax(320),x,ratio_m)
            legend(ax(320),{'MBDyn' 'Analytical'})
        end
    end
    
    fid = kk*10+1;
    lastfid = fid;
end
%%
return
FolderName=strcat(pwd,'\images');
FoldeName=char(FolderName);
% FolderName = '\Images';   % Your destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  saveas(FigHandle, fullfile(FolderName, [FigName, '.eps']),'epsc');    %<---- 'Brackets'
end
