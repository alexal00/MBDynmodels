% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
%
% Post-processor to obtain the corresponding figures from Chapter 4 of the thesis.
% NOTE: This script writes outputs to SpectraAermech.
close all; clear all; clc
%% Plotting parameters
setPlot
% cmap = colormap(lines(8));
% close all
% Color =cell(1,8);
% for ii=1:8
%    Color{ii}= cmap(ii,:);
% end
% Marker = {'o','+','x','s','d','^','*','.'};
% LineStyle = {'-' '--' ':' '-.' '-' '--' ':' '-.'};
% plt.Color = Color; plt.Marker = Marker; plt.LineStyle;
% step = 128;
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
nper = 50;      % [-], number of periods to run
R = 5.5;        % [m], rotor radius
e = 0.3048;     % [m], hinge position from hub
nuxi = sqrt(3/2*e/(R-e));

% Activate or deactivate additional blade-hub damper in blades
% 0 : inactive
% 1 : active
% C_b2h = b2hdamp*C_xi(gamma~0.8)
b2hdamp = 1.;

% Use swashplate or not
% 0 : inactive
% 1 : active
sw = 0;

% Activate or deactivate blade degrees of freedom
% Since the degrees of freedom are given by means of a total joint in MBDyn
% 0 means that the degrre of freedom is not constrained and viceversa for
% 1.
% 0: allowed
% 1: clamped
flap = 1;
lag = 0;
% The pitch dof is constrained if no swashplate is present to prevent
% static divergence and therefore numerical issues.
if sw ==1
    pitch = 0;
else
    pitch = 1;
end


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
% 1 : Original combination of parameters
% 2 : Combination that makes Kxidelta = 1 for Nb=4
opt = 1;

% Plot Kinematic couplings for the i2b
% 1 : Plot results and display couplings
% 2 : Dont plot results
kincoup = 1 ;

% Number of harmonics to study
nF = 10 ;
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
damp = {'std' 'ib' 'i2b'};
% damp = {'ib' 'i2b'};
% damp = {'std' 'i2b'};
% damp = {'i2b'};
% damp = {'std'};
% damp = {'ib'};
% damp = {'std' 'ib'};


% Ranges for table wriing in excel
if sw ==1
    rng = {'B3:D7' 'E3:G7' 'H3:J7'};
    rng2 = {'D5:D8' 'E5:E8'};
else
    rng = {'B12:B16' 'C12:C16' 'D12:D16'};
    rng2 = {'B5:B8' 'C5:C8'};
end

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
fig(110).Name = ['PSDDeltaxidotcomp_' dof];
ax(110+Nb) = gca;
hold(ax(110+Nb),"on")
grid(ax(110+Nb),"on")
xlim(ax(110+Nb),[0,9])
xticks(ax(110+Nb),0:9)
xlabel(ax(110+Nb),'harmonic')
ylabel(ax(110+Nb),['$PSD(\Delta\dot{\xi}_' num2str(Nb) ')$'])
grid(ax(110+Nb),"on")


% Evolution of harmonic content with revolutions
% idx = zeros(length(damp),1);
% idx2 = zeros(length(damp),1);
% for kk = 1:length(damp)
%     % with log-scale
%     idx(kk) = 120+10*(kk-1);
%     idx(kk) = int64(idx(kk));
%     fig(idx(kk)) = figure(idx(kk));
%     fig(idx(kk)).Name = ['PSDDxidotcompRevs' damp{kk} 'ls'];
%     tlo(idx(kk)) = tiledlayout('flow');
%     for ii = 1: Nb
%         ax(idx(kk)+ii) = nexttile(tlo(idx(kk)));
%         yscale(ax((idx(kk)+ii)),"log")
%         hold(ax(idx(kk)+ii),"on")
%         grid(ax(idx(kk)+ii),"on")
%     end
%     % without log-scale
%     idx2(kk) = 150+10*(kk-1);
%     idx2(kk) = int64(idx2(kk));
%     fig(idx2(kk)) = figure(idx2(kk));
%     fig(idx2(kk)).Name = ['PSDDxidotcompRevs' damp{kk}];
%     tlo(idx2(kk)) = tiledlayout('flow');
%     for ii = 1: Nb
%         ax(idx2(kk)+ii) = nexttile(tlo(idx2(kk)));
%         hold(ax(idx2(kk)+ii),"on")
%         grid(ax(idx2(kk)+ii),"on")
%     end
% 
% end

% Figure for harmonics in xi comparisson between models.
% fig(180) = figure(180);
% fig(180).Name = 'PSDDeltaxidotcomp_xi';
% tlo(180) = tiledlayout('flow');
% for ii = 1: Nb
%     ax(180+ii) = nexttile(tlo(180));
%     hold(ax(180+ii),"on")
%     grid(ax(180+ii),"on")
% end

% Figure for harmonic comparisson between models
% fig(200) = figure(200);
% fig(200).Name = 'PSDDeltaxidotcompls';
% tlo(200) = tiledlayout('flow');
% for ii = 1: Nb
%     ax(200+ii) = nexttile(tlo(200));
%     yscale(ax((200+ii)),"log")
%     hold(ax(200+ii),"on")
%     grid(ax(200+ii),"on")
% end

% Figure for harmonic comparisson between models, without log-scale
fig(210) = figure(210);
fig(210).Name = ['PSDmxicomp' 'Nb' int2str(Nb) '_' dof];
ax(210+Nb) = gca;
hold(ax(210+Nb),"on")
grid(ax(210+Nb),"on")
xlim(ax(210+Nb),[0,9])
xticks(ax(210+Nb),0:9)
xlabel(ax(210+Nb),'harmonic')
ylabel(ax(210+Nb),['$PSD(M_{\xi_' num2str(Nb) '})$'])
grid(ax(210+Nb),"on")

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
        ylabel(ax(320),'$M_{\xi}^{i2b}/M_{\xi}^{std}$')
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
% Power over one oscillation cycle
fig(380) = figure(380);
fig(380).Name = ['PdNb' int2str(Nb) '_' dof];
ax(380) = gca;
hold(ax(380),"on")
grid(ax(380),"on")
xlabel(ax(380), 'revolutions')
ylabel(ax(380),'$P_d$ [W]')
legend(ax(380))

% Figure for comparisson between energy dissipated for all arrangements
fig(400) = figure(400);
fig(400).Name = ['Ed' 'Nb' int2str(Nb) '_' dof];
tlo(400) = tiledlayout('flow');
% Dissipated energy
ax(401) = nexttile(tlo(400));
hold(ax(401),"on")
grid(ax(401),"on")
xlim(ax(401),[-1,9])
xticks(ax(401),0:9)
xlabel(ax(401),'harmonic')
ylabel(ax(401),'$E_d$ [J]')
grid(ax(401),"on")
% Percentage of energy per harmonic
ax(402) = nexttile(tlo(400));
hold(ax(402),"on")
grid(ax(402),"on")
xlim(ax(402),[-1,9])
xticks(ax(402),0:9)
xlabel(ax(402),'harmonic')
ylabel(ax(402),'$E_d^n/E_d$ [\%]')
grid(ax(402),"on")

% Ratio of energy dissipated between non-conventional and standard
if any(strcmp(damp, 'std'))
    if any(strcmp(damp, 'i2b'))
        fig(500) = figure(500);
        fig(500).Name = ['Edrati2bNb' int2str(Nb) 'opt' int2str(opt) '_' dof];
        ax(510) = gca;
        hold(ax(510),"on")
        grid(ax(510),"on")
        xlabel(ax(510),'harmonic')
        ylabel(ax(510),'$E_d^{n,i2b}/E_d^{n,std}$')
        xlim(ax(310),[0 10])
        xticks(ax(310),1:9)
    end

    if any(strcmp(damp, 'ib'))
        fig(520) = figure(520);
        fig(520).Name = ['EdratibNb' int2str(Nb) 'opt' int2str(opt) '_' dof];
        ax(530) = gca;
        hold(ax(530),"on")
        grid(ax(530),"on")
        xlabel(ax(530),'harmonic')
        ylabel(ax(530),'$E_d^{n,ib}/E_d^{n,std}$')
        xlim(ax(360),[0 10])
        xticks(ax(360),1:9)
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

Ed = zeros(1,length(damp));
contribution_percent = zeros(length(damp),nF+1);
harmonic_energy = zeros(length(damp),nF+1);

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
    
    % angular velocity
    %%% Omega = 40; % rad/s
    % we could get it from the motion of the hub, assuming it starts with nominal RPM
    Omega = ncread(fn, 'node.struct.10000.Omega', [3, 1], [1, 1]);
    
    % compute azimuth vector from database
    time = ncread(fn, 'time');
    % ND time in revolutions of the rotor
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
    
    % Steps per period from MBDyn
    nT = 256;
    % Number of revultions
    revN = fix(length(psi_nd)/nT) - 1;
    % Indexes for the last period
    periodidx = revN*nT+[1:nT];
    
    T_xi = 2*pi/(nuxi*Omega);
    % revNxi = (time(end)/(T_xi));
    % nTxi = fix(length(psi_nd)/revNxi);
    % revNxi = fix(length(psi_nd)/nTxi)-2;

    % % #1, Lead-lag of single blade
    % % Figure = kk*10+1
    % % Tiled layout = kk*10+1
    % fig(lastfid) = figure(lastfid);
    % fig(lastfid).Name = [fignames 'llsb'];
    % tlo(lastfid) = tiledlayout('flow');
    % fid = fid +1 ;
    % 
    % % #2, Lead-lag of all blades
    % % Figure = kk*10+1+1
    % fig(lastfid+1) = figure(lastfid+1);
    % fig(lastfid+1).Name = [fignames 'llallb'];
    % ax(2) = gca;
    % hold(ax(2),'on')
    % fid = fid +1 ;

    psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
    xid = zeros(nb,length(time));
    for ii = 1:nb
        blade = int2str(10000+1000*ii + 30);
        xid(ii,:) = ncread(fn, ['elem.joint.',blade,'.Omega'],[3, 1], [1, Inf]);

        % psi_b(ii,:) = psi+ii*deltapsi;
        % xii(ii,:) = ncread(fn, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);
        

        % ax(1) = nexttile(tlo(lastfid));
        % hold(ax(1),"on");
        % functionplot(ax(1),psi_nd,-180/pi*xii(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\xi_' int2str(ii) '$'])
        % grid(ax(1),"on")
        % xlabel(ax(1),'revolutions [-]');
        % ylabel(ax(1),'$\xi$ [deg]')

        % functionplot(ax(2),psi_nd,-180/pi*xii(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\xi_' int2str(ii) '$'])
        % legend(ax(2))
        % grid(ax(2),"on")
        % xlabel(ax(2),'revolutions [-]')
        % ylabel(ax(2),'$\xi$ [deg]')
    end

    
    % % # 3, MBC
    % % Figure = kk*10++1+2
    % xi_NR_MB = MBC(psi_b,xii);
    % fig(lastfid+2)=figure(lastfid+2); 
    % fig(lastfid+2).Name = [fignames 'qNR'];
    % % fid = fid+1;
    % qNR = NRdof(nb);
    % ax(3)= subplot(1,2,1);
    % % Plot collective DOF
    % if mod(nb,2)==0
    %     dof_col = [1 nb];
    %     functionplot(ax(3),psi_nd,xi_NR_MB(1,:),step,Color{1},Marker{1},LineStyle{1},[qNR{1} ', MBDyn'])
    %     % plot(psi_nd, xi_NR_a(1,:),'DisplayName',qNR{1},'Color',Color{1}); hold on
    %     functionplot(ax(3),psi_nd,xi_NR_MB(end,:),step,Color{2},Marker{2},LineStyle{2},[qNR{end} ', MBDyn'])
    %     % plot(psi_nd, xi_NR_a(end,:),'DisplayName',qNR{end},'Color',Color{2}); hold on
    % else
    %     dof_col = 1;
    %     functionplot(ax(3),psi_nd,xi_NR_MB(1,:),step,Color{1},Marker{1},LineStyle{1},[qNR{1} ', MBDyn'])
    %     % plot(psi_nd, xi_NR_a(1,:),'DisplayName',qNR{1},'Color',Color{1}); hold on
    % end
    % xlabel('revolutions [-]');
    % ylabel('Collective $q_{NR}$ [rad]')
    % legend('show','interpreter','latex')
    % grid on
    % 
    % dof_cyc = setdiff(1:nb,dof_col);
    % ax(4)= subplot(1,2,2);
    % for jj = 1:length(dof_cyc)
    %     functionplot(ax(4),psi_nd,xi_NR_MB(dof_cyc(jj),:),step,Color{jj},Marker{jj},LineStyle{jj},[qNR{dof_cyc(jj)} ', MBDyn'])
    %     % plot(t', xi_NR_a(dof_cyc(jj),:),'DisplayName',qNR{dof_cyc(jj)},'Color',Color{jj},'LineStyle',LineStyle{jj}); hold on
    % end
    % xlabel('revolutions [-]');
    % ylabel('Cyclic $q_{NR}$ [rad]')
    % legend
    % grid on

    deltaxidot = zeros(nb,length(time));
    for ii = 1:nb
        % Select blade number
        ant = mod(ii - 2, nb) + 1;  % Calculate previous index

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

        % Plot Kinematic couplings
        if kincoup == 1
            % The model only supports the case of the i2b
            if strcmp(damp{kk},'i2b')
                % Relevant node positions
                D_NOD = 80;
                F_NOD = 90;
                HUB = 10000;
                BLADE_FLAP_OFFSET = 30;
                if ii == nb % Plot results only for the last blade
                    % Obtain relevant hinge number
                    prevarm = int2str(HUB+1000*ant + D_NOD);               % D arm
                    actual = int2str(HUB+1000*ii + BLADE_FLAP_OFFSET);  % blade hinge
                    nextarm = int2str(HUB+1000*int_b + F_NOD);             % F arm
            
                    % Legend labels
                    nm1=['xi_' num2str(ii)];
                    leg1 = latexfmt(nm1);
                    % Arm 1 of damper i-1
                    nm2=['delta_D' num2str(ant)];
                    leg2 = latexfmt(nm2);
            
                    nmt =['theta_' num2str(ii)];
                    legt = latexfmt(nmt);
                    % Arm 2 of damper i+1
                    nm3=['delta_F' num2str(nxt)];
                    leg3 = latexfmt(nm3);
            
                    nmb =['beta_' num2str(ii)];
                    legb = latexfmt(nmb);
            
                    % Orientation of damper arms
                    % * connection with previous blade
                    delta_d = ncread(fn, ['elem.joint.',prevarm,'.Phi'],[3, 1], [1, Inf]);
                    delta_dp = ncread(fn, ['elem.joint.',prevarm,'.Omega'],[3, 1], [1, Inf]);
                    % * connection with next blade
                    delta_f = ncread(fn, ['elem.joint.',nextarm,'.Phi'],[3, 1], [1, Inf]);
                    delta_fp = ncread(fn, ['elem.joint.',nextarm,'.Omega'],[3, 1], [1, Inf]);
                    % Lead lag angle
                    xii = ncread(fn, ['elem.joint.',actual,'.Phi'],[3, 1], [1, Inf]);
                    xiip = ncread(fn, ['elem.joint.',actual,'.Omega'],[3, 1], [1, Inf]);
                    % Pitch angle and velocity
                    theta = ncread(fn, ['elem.joint.',actual,'.Phi'],[1, 1], [1, Inf]);
                    thetap = ncread(fn, ['elem.joint.',actual,'.Omega'],[1, 1], [1, Inf]);
                    % Flap angle and velocity
                    beta = -ncread(fn, ['elem.joint.',actual,'.Phi'],[2, 1], [1, Inf]);
                    betap = ncread(fn, ['elem.joint.',actual,'.Omega'],[2, 1], [1, Inf]);
            
            
                    % Linearized values
                    [Ktt,Kbb,Ktb,Kxx] = kinconsti2b(e,geom.a,geom.ca,geom.f,geom.cf,Nb);
                    % Call the function
                    % * LS for delta_d
                    [K, contributions] = LSkinconsti2bMBDyn(xii, theta, beta, delta_d);
                    % * LS for delta_f
                    [K2, contributions2] = LSkinconsti2bMBDyn(xii, theta, beta, delta_f);

                    % Display the results
                    disp('Coefficients:');
                    disp(['K_{xd} = ', num2str(K(1))]);
                    disp(['K_{xi} = ', num2str(K(2))]);
                    disp(['K_{td} = ', num2str(K(3))]);
                    disp(['K_{bd} = ', num2str(K(4))]);
                    disp(['K_{tb} = ', num2str(K(5))]);
            
                    disp('Percentage contributions:');
                    disp(['Contribution of xi: ', num2str(contributions(1)), '%']);
                    disp(['Contribution of xi^2: ', num2str(contributions(2)), '%']);
                    disp(['Contribution of theta^2: ', num2str(contributions(3)), '%']);
                    disp(['Contribution of beta^2: ', num2str(contributions(4)), '%']);
                    disp(['Contribution of theta * beta: ', num2str(contributions(5)), '%']);
            
            
                    % Plot Kinematic couplings influence
                    fid = fid+1;
                    fhandle = figure(fid);
                    fhandle.Name = ['Nb' int2str(Nb) damp{kk} 'KINCOUP_' dof];
                    subplot(1,2,1)
                    plot(psi_nd(periodidx),xii(periodidx),'DisplayName','$\xi$'); hold on
                    plot(psi_nd(periodidx),xii(periodidx).^2,'DisplayName','$\xi^2$')
                    plot(psi_nd(periodidx),beta(periodidx).^2,'DisplayName','$\beta^2$')
                    plot(psi_nd(periodidx),theta(periodidx).^2,'DisplayName','$\theta^2$')
                    plot(psi_nd(periodidx),theta(periodidx).*beta(periodidx),'DisplayName','$\theta\beta$')
                    plot(psi_nd(periodidx),delta_d(periodidx),'DisplayName','$\delta_D$')
                    xlabel('revolutions')
                    ylabel('Angles [rad or rad\textsuperscript{2}]')
                    legend
                    grid on

                    subplot(1,2,2)
                    plot(psi_nd(periodidx),xii(periodidx),'DisplayName','$\xi$'); hold on
                    plot(psi_nd(periodidx),xii(periodidx).^2,'DisplayName','$\xi^2$')
                    plot(psi_nd(periodidx),beta(periodidx).^2,'DisplayName','$\beta^2$')
                    plot(psi_nd(periodidx),theta(periodidx).^2,'DisplayName','$\theta^2$')
                    plot(psi_nd(periodidx),theta(periodidx).*beta(periodidx),'DisplayName','$\theta\beta$')
                    plot(psi_nd(periodidx),delta_f(periodidx),'DisplayName','$\delta_F$')
                    xlabel('revolutions')
                    ylabel('Angles [rad or rad\textsuperscript{2}]')
                    legend
                    grid on
                    fid = fid +1;
                end
            end
        end

        % % #4, Relative velocity at dampers
        % Figure = kk*10+1+3
        % Tiled layout = kk*10+1+3
        % fig(lastfid+3) = figure(lastfid+3); 
        % fig(lastfid+3).Name = [fignames 'DeltaXidot'];
        % tlo(lastfid+3) = tiledlayout('flow');
        % fid=fid+1;
        %
        % steps = nT;
        % xid_m = movmean(deltaxidot(ii,:),steps);
        % calculate_rms = @(funct) sqrt(mean(funct.^2));
        % 
        % % Calculate RMS for each time step
        % rms_values(ii,:) = arrayfun(@(n) calculate_rms(xid_m(1:n)), ...
        %     1:length(deltaxidot));
        % defhing = int2str(10000+1000*int_b + 31);
        % if ~strcmp(damp{kk},'ib')
        %     ddeltadot = ncread(fn, ['elem.joint.',defhing,'.Omega'],[3, 1], [1, Inf]);
        % else
        %     ddeltadot = ncread(fn, ['elem.joint.',defhing,'.lP']);
        % end
        % 
        % if strcmp(damp{kk},'i2b')
        %     k = -1;
        %     kappa = geom.Kxidelta;
        % elseif strcmp(damp{kk},'ib')
        %     k = 1;
        %     kappa = geom.Kxil;
        % else
        %     k= 1;
        %     kappa = 1;
        % end
        % 
        % dxidot_mbd(ii,:) = 1/kappa*k*ddeltadot;
        % 
        % 
        % ax(1) = nexttile(tlo(lastfid+3));
        % functionplot(ax(1),psi_nd,deltaxidot(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\Delta\dot{\xi}_' int2str(ii) '$'])
        % functionplot(ax(1),psi_nd,dxidot_mbd(ii,:),step,'k','none',':',['$\Delta\dot{\delta}_' int2str(ii) '/K_{\xi\delta}$, MBDyn'])
        % functionplot(ax(1),psi_nd,rms_values(ii,:),step,'r','none','-',['$\mathrm{RMS}(\dot{\xi}_' int2str(ii) ')$'])
        % grid(ax(1),"on")
        % legend(ax(1))
        % xlabel(ax(1),'revolutions [-]');
        % ylabel(ax(1),['$\Delta\dot{\xi}_' int2str(ii) '$ [rad/s]']) 
        % aux = 1;
    end
    
    % Truncate signals to last period
    time_o = time(periodidx);
    [P1xid,f_o]=computeFFT(deltaxidot(:,periodidx)',time_o);
    f_ond = f_o./(Omega/(2*pi));

    % time_xi = time(nTxi*revNxi+[1:nTxi*2]);
    % [P1xi_xi,f_xi]=computeFFT(deltaxidot(:,nTxi*revNxi+[1:nTxi*2])',time_xi);
    % f_xind = f_xi./(nuxi*Omega/(2*pi));
    
    % % #5
    % Figure = kk*10+1+4
    % Tiled layout = kk*10+1+4
    % fig(lastfid+4)=figure(lastfid+4); 
    % fig(lastfid+4).Name = [fignames 'PSDDxid'];
    % tlo(lastfid+4) = tiledlayout('flow');
    % % fid = fid+1;
    % for ii = 1:size(P1xi,2)
    %     ax(2) = nexttile(tlo(lastfid+4));
    %     stem(ax(2),f_ond,P1xi(:,ii));
    %     xlim(ax(2),[0 10])
    %     xlabel(ax(2),'f [Hz]')
    %     ylabel(ax(2),'$\mathrm{PSD}(\dot{\xi})$')
    %     grid(ax(2),'on')
    % end

    % Moments in the damper
    % #6
    % Figure = kk*10+1+5
    % Tiled layout = kk*10+1+5
    % fig(lastfid+5)=figure(lastfid+5); 
    % fig(lastfid+5).Name = [fignames 'Mxi'];
    % tlo(lastfid+5) = tiledlayout('flow');
    % fid = fid+1;

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
        % Select blade number
        ant = mod(ii - 2, nb) + 1;  % Calculate previous index

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

        antblade = 10000+1000*ant + 30;
        nxtblade = 10000+1000*int_b + 30;

        if ~strcmp(damp{kk},'i2b')
            damper = ii;  % Calculate next blade
        elseif strcmp(damp{kk},'i2b')
            damper = mod(ii, Nb) + 1;  % Calculate the intermediate blade
        end

        defhing = int2str(10000+1000*damper + 31);
        if strcmp(damp{kk},'std')
            m_xi(ii,:) = ncread(fn, ['elem.joint.',defhing,'.m'],[3, 1], [1, Inf]);
            delta_xi(ii,:) = ncread(fn, ['elem.joint.',defhing,'.Omega'],[3, 1], [1, Inf]);
            kappa = 1;
        elseif strcmp(damp{kk},'i2b')
            m_xi(ii,:) = -1*ncread(fn, ['elem.joint.',defhing,'.m'],[3, 1], [1, Inf]);
            delta_xi(ii,:) = -1*ncread(fn, ['elem.joint.',defhing,'.Omega'],[3, 1], [1, Inf]);
            kappa = geom.Kxidelta;

        else
            Force_xi(ii,:) = ncread(fn, ['elem.joint.',defhing,'.f'],[1, 1], [1, Inf]);
            delta_xi(ii,:) = ncread(fn, ['elem.joint.',defhing,'.lP']);
            m_xi(ii,:) = Force_xi(ii,:);%*geom;
            kappa = geom.Kxil;
        end

        if ii == nb && strcmp(damp{kk},'i2b')
            deltadelta = ncread(fn, ['elem.joint.',int2str(10000+1000*ii + 80),'.Omega'],[3, 1], [1, Inf]);
            dxi = ncread(fn, ['elem.joint.',int2str(nxtblade),'.Omega'],[3, 1], [1, Inf]);

            dxixip = ncread(fn, ['elem.joint.',int2str(nxtblade),'.Phi'],[3, 1], [1, Inf]).*...
                    ncread(fn, ['elem.joint.',int2str(nxtblade),'.Omega'],[3, 1], [1, Inf]);

            dthetathetap = ncread(fn, ['elem.joint.',int2str(nxtblade),'.Phi'],[1, 1], [1, Inf]).*...
                ncread(fn, ['elem.joint.',int2str(nxtblade),'.Omega'],[1, 1], [1, Inf]);

            dbetathetap = -ncread(fn, ['elem.joint.',int2str(nxtblade),'.Phi'],[2, 1], [1, Inf]).*...
                ncread(fn, ['elem.joint.',int2str(nxtblade),'.Omega'],[1, 1], [1, Inf])-...
                ncread(fn, ['elem.joint.',int2str(nxtblade),'.Omega'],[2, 1], [1, Inf]).*...
                ncread(fn, ['elem.joint.',int2str(nxtblade),'.Phi'],[1, 1], [1, Inf]);

            dbetabetap = ncread(fn, ['elem.joint.',int2str(nxtblade),'.Phi'],[2, 1], [1, Inf]).*...
                ncread(fn, ['elem.joint.',int2str(nxtblade),'.Omega'],[2, 1], [1, Inf]);
            Mdxi = geom.Kxidelta*dxi;%+2*K(2)*dxixip;
            % disp(['K_{xd} = ', num2str(K(1))]);
            % disp(['K_{xi} = ', num2str(K(2))]);
            % disp(['K_{td} = ', num2str(K(3))]);
            % disp(['K_{bd} = ', num2str(K(4))]);
            % disp(['K_{tb} = ', num2str(K(5))]);
            Mtti = 2*abs(Ktt)*dthetathetap;
            Mtbi = abs(K(5))*dbetathetap;
            Mbbi = 2*Kbb*dbetabetap;
            taylor = Mdxi(periodidx)+Mtti(periodidx)+Mtbi(periodidx)+Mbbi(periodidx);
            fig(fid) = figure(fid);
            fig(fid).Name = ['AnalyticalTaylor' dof];
            fid = fid+1;
            plot(psi_nd(periodidx),deltadelta(periodidx),'DisplayName','$\dot{\delta}_{D}$'); hold on
            plot(psi_nd(periodidx),Mdxi(periodidx),'DisplayName','$K_{\xi\delta}\dot{\xi}$')
            plot(psi_nd(periodidx),Mtti(periodidx),'DisplayName','$2K_{\theta^2\delta}\theta\dot{\theta}$')
            plot(psi_nd(periodidx),Mtbi(periodidx),'DisplayName','$K_{\beta\theta\delta}(\beta\dot{\theta}+\dot{\beta}\theta)$')
            plot(psi_nd(periodidx),Mbbi(periodidx),'DisplayName','$K_{\beta^2\delta}\beta\dot{\beta}$')
            plot(psi_nd(periodidx),taylor,'DisplayName','$\dot{\delta}_{D}^{Taylor}$')
            xlabel('revolutions')
            ylabel('$\dot{\delta}_D$ [rad/s]')
            legend
            grid on
        end

        % ax(2) = nexttile(tlo(lastfid+5));
        % functionplot(ax(2),psi_nd,m_xi(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$m_{\xi,' int2str(ii) '}$'])
        % grid(ax(2),"on")
        % legend(ax(2))
        % xlabel(ax(2),'revolutions [-]');
        % ylabel(ax(2),['$m_{\xi,' int2str(ii) '}$ [Nm]'])
    
        % ax(3) = nexttile(tlo(lastfid+6));
        % functionplot(ax(3),psi_nd,c_xic(kk).*kappa*deltaxidot(ii,:),step,Color{ii},Marker{ii},'-','$C_{xi}\Delta\dot{\xi}$')
        % functionplot(ax(3),psi_nd,m_xi(ii,:),step,'k','none','-.','$m_{\xi}$, MBDyn')
        % grid(ax(3),"on")
        % legend(ax(3))
        % xlabel(ax(3),'revolutions [-]');
        % ylabel(ax(3),['$m_{\xi,' int2str(ii) '}$'])


        [P1m,~]=computeFFT(m_xi(:,periodidx)',time_o);
        [P1dD,~] = computeFFT(delta_xi(:,periodidx)',time_o);
        % ax(3) = nexttile(tlo(lastfid+7));
        % stem(ax(3),f_ond,P1m(:,ii));
        % xlim(ax(3),[0 10])
        % xlabel(ax(3),'f [Hz]')
        % ylabel(ax(3),'$\mathrm{PSD}(m_{\xi})$') 
    end

    % Obtaine total energy dissipated, and contribution of harmonics
    [Ed(kk), harmonic_energy(kk,:), contribution_percent(kk,:)] = calculate_energy_dissipation(m_xi(nb,periodidx)',delta_xi(nb,periodidx)',time_o,nF);
    
    disp('Energy dissipated:')
    disp(Ed)
    % Plot the power dissipated over one cycle
    plot(ax(380),psi_nd(periodidx),m_xi(nb,periodidx).*delta_xi(nb,periodidx),'DisplayName',damp{kk})


    % Plot the energy dissipated by the damper and the contribution of each
    % harmonic
    fig(fid) = figure(fid); 
    fig(fid).Name = [damp{kk} 'HarmonicPP' dof]; fid = fid+1;
    tlo(1) = tiledlayout(3,2);
    ax(100) = nexttile(tlo(1),1,[1 2]);
    % Plot the norm of the n-th harmonic
    stem(ax(100),f_ond,P1m(:,nb),'filled')
    xlim(ax(100),[0,9])
    xticks(ax(100),0:9)
    xlabel(ax(100),'harmonic')
    grid(ax(100),"on")
    hold(ax(100),"on")
    ax(110) = nexttile(tlo(1),3,[1 2]);
    stem(ax(110),f_ond,P1dD(:,nb),'filled')
    xlim(ax(110),[0,9])
    xticks(ax(110),0:9)
    xlabel(ax(110),'harmonic')
    grid(ax(110),"on")
    hold(ax(110),"on")

    if strcmp(damp{kk},'ib')
        ylbl1 = '$\mathrm{PSD}(F_{\xi}^{ib})$ [N]';
        ylbl2 = '$\mathrm{PSD}(\dot{\Delta})$ [m/s]';
    else
        ylbl1 = ['$\mathrm{PSD}(M_{\xi}^{' damp{kk} '})$ [Nm]'];
        ylbl2 = '$\mathrm{PSD}(\dot{\Delta})$ [rad/s]';
    end
    ylabel(ax(100),ylbl1)
    ylabel(ax(110),ylbl2)
    
    ax(120) = nexttile(tlo(1),5,[1 2]);
    hold(ax(120),"on")
    bar(ax(120),0:nF,contribution_percent(kk,:))
    xlim(ax(120),[0,9])
    xticks(ax(120),0:9)
    xlabel(ax(120),'harmonic')
    ylabel(ax(120),'$E_d^n/E_d$ [\%]')
    grid(ax(120),"on")
    % legend(ax(120));

     
    % Plot one-sided espectrum of lead-lag angular speed without log
    % scale
    stem(ax(110+Nb),f_ond',P1xid(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
    legend(ax(110+Nb))
    if strcmp(damp{kk},'ib')
        stem(ax(210+Nb),f_ond',geom.Kxil*P1m(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
    else
        stem(ax(210+Nb),f_ond',P1m(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
    end

    legend(ax(210+Nb))
       
    label = joint_label;
    % transform in deg, opposite sign for flap and lead-lag
    theta = 180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [1, 1], [1, Inf]);
    beta = -180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [2, 1], [1, Inf]);
    xi = -180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [3, 1], [1, Inf]);
    % % #9, Pith-flap lag
    % fig(lastfid+8)=figure(lastfid+8); 
    % fig(lastfid+8).Name = [fignames 'PFL'];
    % % fid = fid+1;
    % plot(psi_nd, theta, ...
    %     psi_nd, beta, ...
    %     psi_nd, xii);
    % xlabel('revolutions [-]');
    % ylabel([blade_txt, ' angles, [deg]']);
    % legend('$\theta$', '$\beta$', '$\xi$');
    % grid on

    % Write the corresponding angles in the excel sheet
    [P1theta,~] = computeFFT(theta(periodidx)',time_o);
    [P1beta,~] = computeFFT(beta(periodidx)',time_o);
    [P1xi,~] = computeFFT(xi(periodidx)',time_o);
    if sw == 1
        M = [P1theta(1:5),P1beta(1:5),P1xi(1:5)];
    else
        M = P1xi(1:5);
    end
    writematrix(M,"SpectraAermech.xlsx",'Sheet','Hoja1','Range',rng{kk})
    
    % Plot the damping ratios
    if any(strcmp(damp,'std'))
        x = 1:9;
        if strcmp(damp{kk},'std')
            PSDxid_std = P1xid(2:10,1);
            PSDMxi_std = P1m(2:10,1);
            Edn_std = harmonic_energy(kk,:);
            if Edn_std(1) <= 1e-2
                Edn_std(1) = 1;
            end
        elseif strcmp(damp{kk},'ib')
            Edn_ib = harmonic_energy(kk,:);
            ratdxi =  P1xid(2:10,1)./PSDxid_std;
            % ratmxi =  P1m(2:10,1)./PSDMxi_std;
            ratmxi =  geom.Kxil*P1m(2:10,1)./PSDMxi_std;
    
            % Analytical values
            ratdxi_a = sqrt(2*(1-cos(x*dpsi)));
            ratmxi_a = sqrt(2*(1-cos(x*dpsi)))/(2*(1-cos(dpsi)));
            
            ratio_dxi = [ratdxi';ratdxi_a];
            ratio_m = [ratmxi';ratmxi_a];
            writematrix(ratmxi(1:4),"SpectraAermech.xlsx",'Sheet','Hoja3','Range',rng2{1})
            bar(ax(360),x,ratio_dxi)
            legend(ax(360),{'MBDyn' 'Analytical'})
            bar(ax(370),x,ratio_m)
            legend(ax(370),{'MBDyn' 'Analytical'})

            bar(ax(530),0:nF,Edn_ib./Edn_std)
        elseif strcmp(damp{kk},'i2b')
            Edn_i2b = harmonic_energy(kk,:);
            ratdxi =  P1xid(2:10,1)./PSDxid_std;
            ratmxi =  P1m(2:10,1)./PSDMxi_std;
            
            % Analytical values
            ratdxi_a = sqrt(2*(1-cos(x*2*dpsi)));
            ratmxi_a = sqrt(2*(1-cos(x*2*dpsi)))/(2*(1-cos(2*dpsi)))/geom.Kxidelta;
            
            ratio_dxi = [ratdxi';ratdxi_a];
            ratio_m = [ratmxi';ratmxi_a];
            writematrix(ratmxi(1:4),"SpectraAermech.xlsx",'Sheet','Hoja3','Range',rng2{2})
            bar(ax(310),x,ratio_dxi)
            legend(ax(310),{'MBDyn' 'Analytical'})
            bar(ax(320),x,ratio_m)
            legend(ax(320),{'MBDyn' 'Analytical'})

            bar(ax(510),0:nF,Edn_i2b./Edn_std)
        end
    end
    
    fid = kk*10+1;
    lastfid = fid;
end

bar(ax(401),0:nF,harmonic_energy')
legend(ax(401),damp)

bar(ax(402),0:nF,contribution_percent')
legend(ax(402),damp)
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
