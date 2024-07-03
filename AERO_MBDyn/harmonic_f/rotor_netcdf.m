% Author: Alejandro Alvaro, 2023-2024
% execute MBDyn generic rigid rotor model and visualize output
%
% The rotor is excited by means of a harmonic load applied at the lead-lag
% hinges. NO aerodynamics are involved.
%
% Preliminary verification step to validate the expected frequency content
% of relative blade velocity and loads in the damper.
close all; clear all; clc
%% Plotting parameters
setPlot
lastfid = fid;
%%
deg2rad = pi/180;
N_b = 4;
% Environment variables for MBDyn
omega = 40;  %[rad/s]
variables =['"' 'const real OMEGA_100 = ' num2str(omega) '"'];
setenv('MBDYNVARS', variables);

% Filepaths for MBDyn
% In preffix
pref_in = '/home/ale/thesis/AERO_MBDyn/harmonic_f/rotor_aero_';
% Out preffix
pref_out = '/home/ale/thesis/AERO_MBDyn/harmonic_f/';
% Shell file for MBDyn
shfile = '/home/ale/mbdyn/mbdyn/mbdyn ';


% Select damping type
damp = {'std' 'ib' 'i2b'};
% damp = {'i2b'};


% Select number of blades
Nb = 4;


% Figure for harmonic comparisson between models
% #5
fig(100) = figure(100);
fig(100).Name = 'PSDDeltaxidotcomp';
for ii = 1: Nb
    ax(100+ii) = subplot(2,2,ii);
    hold on
end

% Figure for harmonic comparisson between models
% #6
fig(110) = figure(110);
fig(110).Name = 'PSDDeltaxidotcompLS';
for ii = 1: Nb
    ax(110+ii) = subplot(2,2,ii);
    yscale(ax((110+ii)),"log")
    hold on
end

% Correction factor for the same equivalent damping
c_xi = 4067.5;
c_xic = zeros(1,length(damp));
dpsi = 2*pi/Nb;
for kk = 1:length(damp)
    if strcmp(damp{kk},'std')
        fact = 1;
    elseif strcmp(damp{kk},'ib')
        fact=1/(2*(1-cos(dpsi)));
    elseif strcmp(damp{kk},'i2b')
        fact=1/(2*(1-cos(2*dpsi)));
    end
    c_xic(kk) = c_xi*fact;
end


for kk=1:length(damp)
    f_in = [pref_in damp{kk} '_Nb' int2str(N_b) '.mbd'];
    fn_base = ['hhh' damp{kk}];
    f_out = [pref_out fn_base];
        
    disp('executing MBDyn...');


    [rc, errmsg] = system(['wsl ' 'MBDYNVARS=' variables ' ' shfile  f_in ' -o ' f_out])
    % [rc, errmsg] = system(['./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
    if (rc ~= 0),
        error(errmsg);
    end
    disp('   ... done');

    % filename
    fn = [f_out '.nc'];
    disp(sprintf('reading output from file ''%s''', fn));
    
    % angular velocity
    %%% Omega = 40; % rad/s
    % we could get it from the motion of the hub, assuming it starts with nominal RPM
    Omega = ncread(fn, 'node.struct.10000.Omega', [3, 1], [1, 1]);
    
    % compute azimuth vector from database
    psi_nd = ncread(fn, 'time')*Omega/(2*pi);
    
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
    
    % datum = 'alpha';
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames datum];
    % fid = fid+1;
    % legend_txt = {};
    % for p = 0:ng-1,
    %     plot(psi_nd, ncread(fn, ['elem.aerodynamic.',label,'.',datum,'_',int2str(p)]));
    %     legend_txt{p+1} = ['point ', int2str(p+1)];
    %     hold on
    % end
    % hold off
    % xlabel('revolutions [-]');
    % ylabel('$\alpha$ [deg]');
    % legend(legend_txt);
    % 
    % datum = 'gamma';
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames datum];
    % fid = fid+1;
    % for p = 0:ng-1
    %     plot(psi_nd, ncread(fn, ['elem.aerodynamic.',label,'.',datum,'_',int2str(p)]));
    %     hold on
    % end
    % hold off
    % xlabel('revolutions  [-]');
    % ylabel('$\gamma$ [deg]');
    % legend(legend_txt);
    % 
    % datum = 'Mach';
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames datum];
    % fid = fid+1;
    % for p = 0:ng-1
    %     plot(psi_nd, ncread(fn, ['elem.aerodynamic.',label,'.',datum,'_',int2str(p)]));
    %     hold on
    % end
    % hold off
    % xlabel('revolutions [-]');
    % ylabel('Mach [-]');
    % legend(legend_txt);
    % 
    % datum = 'cl';
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames datum];
    % fid = fid+1;
    % for p = 0:ng-1
    %     plot(psi_nd, ncread(fn, ['elem.aerodynamic.',label,'.',datum,'_',int2str(p)]));
    %     hold on
    % end
    % hold off
    % xlabel('revolutions [-]');
    % ylabel('$C_L$ [-]');
    % legend(legend_txt);
    % 
    % datum = 'cd';
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames datum];
    % fid = fid+1;
    % for p = 0:ng-1
    %     plot(psi_nd, ncread(fn, ['elem.aerodynamic.',label,'.',datum,'_',int2str(p)]));
    %     hold on
    % end
    % hold off
    % xlabel('revolutions [-]');
    % ylabel('$C_D$ [-]');
    % legend(legend_txt);
    % 
    % datum = 'cm';
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames datum];
    % fid = fid+1;
    % for p = 0:ng-1,
    %     plot(psi_nd, ncread(fn, ['elem.aerodynamic.',label,'.',datum,'_',int2str(p)]));
    %     hold on
    % end
    % hold off
    % xlabel('revolutions [-]');
    % ylabel('$C_M$ [-]');
    % legend(legend_txt);
    
    % label = inflow_label;
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames 'Forces'];
    % fid = fid+1;
    % plot(psi_nd, ncread(fn, ['elem.inducedvelocity.',label,'.f']));
    % xlabel('revolutions [-]');
    % ylabel('$H_H$, $Q_H$, $T_H$ [N]');
    % legend('$H_H$', '$Q_H$', '$T_H$');
    
    nT = 256;
    revN = fix(length(psi_nd)/nT) - 1;
    
    % F_H = ncread(fn, ['elem.inducedvelocity.',label,'.f']);
    % T_H = F_H(3, :);
    % Y = fft(T_H(nT*revN+[1:nT]));
    % Y = Y/(nT/2);
    % Y(1) = Y(1)/2;
    % disp(sprintf('    T_H(0) = %g N', Y(1)));
    % disp(sprintf('    T_H(%d) = %g %+gi N', nb, real(Y(nb+1)), -imag(Y(nb+1))));
    
    % H_H = F_H(1, :);
    % Y = fft(H_H(nT*revN+[1:nT]));
    % Y = Y/(nT/2);
    % Y(1) = Y(1)/2;
    % disp(sprintf('    H_H(0) = %g N', Y(1)));
    % disp(sprintf('    H_H(%d) = %g %+gi N', nb, real(Y(nb+1)), -imag(Y(nb+1))));
    
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames 'Moments'];
    % fid = fid+1;
    % plot(psi_nd, ncread(fn, ['elem.inducedvelocity.',label,'.m']));
    % xlabel('revolutions [-]');
    % ylabel('$\mathcal{L}$, $\mathcal{M}$, $\mathcal{N}$, [Nm]');
    % legend('Roll', 'Pitch', 'Torque');
    
    % M_H = ncread(fn, ['elem.inducedvelocity.',label,'.m']);
    % C_H = M_H(3, :);
    % Y = fft(C_H(nT*revN+[1:nT]));
    % Y = Y/(nT/2);
    % Y(1) = Y(1)/2;
    % disp(sprintf('    C_H(0) = %g Nm', Y(1)));
    % disp(sprintf('    C_H(%d) = %g %+gi Nm', nb, real(Y(nb+1)), -imag(Y(nb+1))));
    % 
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames 'vimean'];
    % fid = fid+1;
    % plot(psi_nd, ncread(fn, ['elem.inducedvelocity.',label,'.UMean']));
    % xlabel('revolutions [-]');
    % ylabel('$\bar{u}$ [m/s]');
    % % legend('u');
    % 
    % UMean = ncread(fn, ['elem.inducedvelocity.',label,'.UMean']);
    % Y = fft(UMean(nT*revN+[1:nT]));
    % Y = Y/(nT/2);
    % Y(1) = Y(1)/2;
    % disp(sprintf('    UMean(0) = %g m/s', Y(1)));
    % disp(sprintf('    UMean(%d) = %g %+gi m/s', nb, real(Y(nb+1)), -imag(Y(nb+1))));
    
    % fig(fid)=figure(fid); 
    % fid = fid+1;
    % plot(psi_nd, ncread(fn, 'elem.inducedvelocity.99.Alpha'));
    % xlabel('revolutions [-]');
    % ylabel('\alpha_H, deg');
    % legend('disk angle');
    
    % fig(fid)=figure(fid); 
    % fid = fid+1;
    % fig(fid).Name = [fignames 'mu'];
    % plot(psi_nd, ncread(fn, ['elem.inducedvelocity.',label,'.Mu']));
    % xlabel('revolutions [-]');
    % ylabel('$\mu$, [-]');
    % % legend('advance ratio');
    % 
    % fig(fid)=figure(fid); 
    % fid = fid+1;
    % fig(fid).Name = [fignames 'lambda'];
    % plot(psi_nd, ncread(fn, ['elem.inducedvelocity.',label,'.Lambda']));
    % xlabel('revolutions [-]');
    % ylabel('$\lambda$, [-]');
    % % legend('inflow ratio');
    % 
    % label = blade_label;
    % fig(fid)=figure(fid);
    % fid = fid+1;
    % plot(psi_nd, ncread(fn, ['node.struct.',label,'.X'], [3, 1], [1, Inf]));
    % xlabel('revolutions [-]');
    % ylabel([blade_txt, ' $z_{tip}$ [m]']);
    % % legend('$z_{tip}$');
    % 
    label = joint_label;
    % fig(fid)=figure(fid); 
    % fig(fid).Name = [fignames 'PFL'];
    % fid = fid+1;
    % % transform in deg, opposite sign for flap and lead-lag
    % plot(psi_nd, 180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [1, 1], [1, Inf]), ...
    %     psi_nd, -180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [2, 1], [1, Inf]), ...
    %     psi_nd, -180/pi*ncread(fn, ['elem.joint.',label,'.Phi'], [3, 1], [1, Inf]));
    % xlabel('revolutions [-]');
    % ylabel([blade_txt, ' angles, [deg]']);
    % legend('$\theta$', '$\beta$', '$\xi$');
    % #1
    fig(fid) = figure(fid);
    fig(fid).Name = [fignames 'llsb'];
    tlo(fid) = tiledlayout('flow');
    fid = fid +1 ;
    % #2
    fig(fid) = figure(fid);
    fig(fid).Name = [fignames 'llallb'];
    ax(2) = gca;
    hold on;
    fid = fid +1 ;

    psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
    for ii = 1:nb
        ax(1) = nexttile(tlo(lastfid));
        hold on;
        psi_b(ii,:) = psi+ii*deltapsi;
        blade = int2str(10000+1000*ii + 30);
        xi(ii,:) = ncread(fn, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);
        xid(ii,:) = ncread(fn, ['elem.joint.',blade,'.Omega'],[3, 1], [1, Inf]);

        functionplot(ax(1),psi_nd,-180/pi*xi(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\xi_' int2str(ii) '$'])
        grid(ax(1),"on")
        xlabel(ax(1),'revolutions [-]');
        ylabel(ax(1),'$\xi$ [deg]')

        functionplot(ax(2),psi_nd,-180/pi*xi(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\xi_' int2str(ii) '$'])
        % plot(psi_nd,-180/pi*xi(ii,:),'DisplayName',['$\xi_' int2str(ii) '$']); 
    end
    legend(ax(2))
    grid(ax(2),"on")
    xlabel(ax(2),'revolutions [-]')
    ylabel(ax(2),'$\xi$ [deg]')
    
    xi_NR_MB = MBC(psi_b,xi);
    % # 3
    fig(fid)=figure(fid); 
    fig(fid).Name = [fignames 'qNR'];
    fid = fid+1;
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
        % plot(psi_nd, xi_NR_a(1,:),'DisplayName',qNR{1},'Color',Color{1}); hold on
    end
    xlabel('revolutions [-]');
    ylabel('Collective $q_{NR}$ [rad]')
    legend('show','interpreter','latex')
    grid on
    
    dof_cyc = setdiff(1:nb,dof_col);
    ax(4)= subplot(1,2,2);
    for jj = 1:length(dof_cyc)
        functionplot(ax(4),psi_nd,xi_NR_MB(dof_cyc(jj),:),step,Color{jj},Marker{jj},LineStyle{jj},[qNR{dof_cyc(jj)} ', MBDyn'])
        % plot(t', xi_NR_a(dof_cyc(jj),:),'DisplayName',qNR{dof_cyc(jj)},'Color',Color{jj},'LineStyle',LineStyle{jj}); hold on
    end
    xlabel('revolutions [-]');
    ylabel('Cyclic $q_{NR}$ [rad]')
    legend
    grid on
        
    % #4
    fig(fid) = figure(fid); 
    fig(fid).Name = [fignames 'DeltaXidot'];
    tlo(fid) = tiledlayout('flow');
    fid=fid+1;
    
    for ii = 1:nb
        xi_d = xid(ii,:) ;
        if strcmp(damp{kk},'std')
            nxt = 0;
        elseif strcmp(damp{kk},'ib')
            nxt = mod(ii, Nb) + 1;  % Calculate next blade
        elseif strcmp(damp{kk},'i2b')
            nxt = mod(ii+1, Nb) + 1;  % Calculate next blade
        end

        if nxt ==0
            deltaxidot(ii,:) = xi_d;
        else
            xi_dnxt = xid(nxt,:);
            deltaxidot(ii,:) = xi_dnxt-xi_d;
        end

        ax(1) = nexttile(tlo(lastfid+3));
        functionplot(ax(1),psi_nd,deltaxidot(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$\Delta\dot{\xi}_' int2str(ii) '$'])
        grid(ax(1),"on")
        legend(ax(1))
        xlabel(ax(1),'revolutions [-]');
        ylabel(ax(1),['$\Delta\dot{\xi}_' int2str(ii) '$ [rad/s]'])

    end
    % psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
    time = ncread(fn, 'time');
    time = time(nT*revN+[1:nT]);
    [P1xi,f]=computeFFT(deltaxidot(:,nT*revN+[1:nT])',time);
    f_nd = f./(Omega/(2*pi));
    % #5
    fig(fid)=figure(fid); 
    fig(fid).Name = [fignames 'PSDDxid'];
    tlo(fid) = tiledlayout('flow');
    fid = fid+1;
    for ii = 1:size(P1xi,2)
        ax(2) = nexttile(tlo(lastfid+4));
        stem(ax(2),f,P1xi(:,ii));
        xlim(ax(2),[0 ceil(10*omega/(2*pi))])
        xlabel(ax(2),'f [Hz]')
        ylabel(ax(2),'$\mathrm{PSD}(\dot{\xi})$')
    end

    if ~strcmp(damp{kk},'ib')
        % #6
        fig(fid)=figure(fid); 
        fig(fid).Name = [fignames 'Mxi'];
        tlo(fid) = tiledlayout('flow');
        fid = fid+1;
        % #7
        fig(fid)=figure(fid); 
        fig(fid).Name = [fignames 'Mxicomp'];
        tlo(fid) = tiledlayout('flow');
        fid = fid+1;

        for ii = 1:nb
            psi_b(ii,:) = psi+ii*deltapsi;
            defhing = int2str(10000+1000*ii + 31);
            m_xi(ii,:) = ncread(fn, ['elem.joint.',defhing,'.m'],[3, 1], [1, Inf]);
        
            ax(2) = nexttile(tlo(lastfid+5));
            functionplot(ax(2),psi_nd,m_xi(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},['$m_{\xi,' int2str(ii) '}$'])
            grid(ax(2),"on")
            legend(ax(2))
            xlabel(ax(2),'revolutions [-]');
            ylabel(ax(2),['$m_{\xi,' int2str(ii) '}$ [Nm]'])
    
            ax(3) = nexttile(tlo(lastfid+6));
            functionplot(ax(3),psi_nd,c_xic(kk).*deltaxidot(ii,:),step,Color{ii},Marker{ii},LineStyle{ii},'$C_{xi}\Delta\dot{\xi}$')
            functionplot(ax(3),psi_nd,m_xi(ii,:),step,'k','none',LineStyle{ii},'$m_{\xi}$, MBDyn')
            grid(ax(3),"on")
            legend(ax(3))
            xlabel(ax(3),'revolutions [-]');
            ylabel(ax(3),['$m_{\xi,' int2str(ii) '}$'])
        end

        [P1m,f]=computeFFT(m_xi(:,nT*revN+[1:nT])',time);

        % #8
        fig(fid)=figure(fid); 
        fig(fid).Name = [fignames 'PSDMxi'];
        tlo(fid) = tiledlayout('flow');
        fid = fid+1;
        for ii=1:nb
            ax(3) = nexttile(tlo(lastfid+7));
            stem(ax(3),f,P1m(:,ii));
            xlim(ax(3),[0 ceil(10*omega/(2*pi))])
            xlabel(ax(3),'f [Hz]')
            ylabel(ax(3),'$\mathrm{PSD}(m_{\xi})$')
        end
    
    end

    
    for ii=1:nb
        % Plot one-sided espectrum of lead-lag angular speed
        stem(ax(100+ii),f_nd',P1xi(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
        xlabel(ax(100+ii),'$\omega/\Omega$ [-]')
        ylabel(ax(100+ii),['$PSD(\Delta\dot{\xi}_' num2str(ii) ')$'])
        xticks(ax(100+ii),[0:7])
        xticklabels(ax(100+ii),{'0','1','2','3','4','5','6','7'})
        xlim(ax(100+ii),[0 7])
        title(ax(100+ii),['Revolution=' num2str(kk)])
        legend(ax(100+ii),'Location','best')
        grid(ax(100+ii),"on")

        stem(ax(110+ii),f_nd',P1xi(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
        xlabel(ax(110+ii),'$\omega/\Omega$ [-]')
        ylabel(ax(110+ii),['$PSD(\Delta\dot{\xi}_' num2str(ii) ')$'])
        xticks(ax(110+ii),[0:7])
        xticklabels(ax(110+ii),{'0','1','2','3','4','5','6','7'})
        xlim(ax(110+ii),[0 7])
        title(ax(110+ii),['Revolution=' num2str(kk)])
        legend(ax(110+ii),'Location','best')
        grid(ax(110+ii),"on")
    end

    % theta = psi_nd(1, :);
    % Y = fft(theta(nT*revN+[1:nT]));
    % Y = Y/(nT/2);
    % Y(1) = Y(1)/2;
    % disp(sprintf('    theta_0 = %g deg', Y(1)*180/pi));
    % disp(sprintf('    A_1_H = %g deg', -real(Y(2))*180/pi));
    % disp(sprintf('    B_1_H = %g deg', +imag(Y(2))*180/pi));
    
    % beta = -psi_nd(2, :);
    % Y = fft(beta(nT*revN+[1:nT]));
    % Y = Y/(nT/2);
    % Y(1) = Y(1)/2;
    % disp(sprintf('    a_0 = %g deg', Y(1)*180/pi));
    % disp(sprintf('    a_1_H = %g deg', -real(Y(2))*180/pi));
    % disp(sprintf('    b_1_H = %g deg', +imag(Y(2))*180/pi));
    
    % fig(fid)=figure(fid); 
    % fid = fid+1;
    % lbl = int2str(10000 + 1000*nb+1);
    % plot3(ncread(fn, ['node.struct.',lbl,'.X'], [1, 1], [1, Inf]), ...
    %     ncread(fn, ['node.struct.',lbl,'.X'], [2, 1], [1, Inf]), ...
    %     ncread(fn, ['node.struct.',lbl,'.X'], [3, 1], [1, Inf]));
    % axis equal
    % title('blade tip path')
    fid = kk*10;
    lastfid = fid;
end
