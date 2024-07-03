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
% fig(110) = figure(110);
% fig(110).Name = 'PSDDeltaxidotcomp';
% tlo(110) = tiledlayout('flow');
% for ii = 1: Nb
%     ax(110+ii) = nexttile(tlo(110));
%     hold(ax(110+ii),"on")
%     grid(ax(110+ii),"on")
%     xlim(ax(110+ii),[0,9])
%     xticks(ax(110+ii),0:9)
%     xlabel(ax(110+ii),'harmonic')
%     ylabel(ax(110+ii),['$PSD(\Delta\dot{\xi}_' num2str(ii) ')$'])
%     grid(ax(110+ii),"on")
% end

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
fig(210) = figure(210+clplt);
fig(210).Name = ['PSDmxicomp' 'Nb' int2str(Nb) '_' dof swash cli ];
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
        fig(300) = figure(300+clplt);
        fig(300).Name = ['rati2bNb' int2str(Nb) 'opt' int2str(opt) '_' dof swash cli];
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
        fig(350) = figure(350+clplt);
        fig(350).Name = ['ratibNb' int2str(Nb) 'opt' int2str(opt) '_' dof swash cli];
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
%% Postprocessing
% Correction factor for the same equivalent damping
% c_xi = 4067.5;
% c_xic = zeros(1,length(damp));
% dpsi = 2*pi/Nb;
% for kk = 1:length(damp)
%     if strcmp(damp{kk},'std')
%         fact = 1;
%     elseif strcmp(damp{kk},'ib')
%         Kxil = geom.Kxil;
% 
%         fact=1/(2*(1-cos(dpsi)))/(Kxil^2);
%         c_d = 1.5567e+05;
%     elseif strcmp(damp{kk},'i2b')
%         Kxidelta = geom.Kxidelta;
%         fact=1/(2*(1-cos(2*dpsi)))/Kxidelta^2;
%     end
%     c_xic(kk) = c_xi*fact;
% end
for kk = 1:length(damp)
    fout = filenc{kk};
    % filename
    fn = [fout '.nc'];
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
    % Preprocessing complete
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
    % fig(lastfid+2)=figure(lastfid+2); 
    % fig(lastfid+2).Name = [fignames 'qNR'];
    % fid = fid+1;
    qNR = NRdof(nb);
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
    

    % % #4, Relative velocity at dampers
    % Figure = kk*10+1+3
    % Tiled layout = kk*10+1+3
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
        % stem(ax(110+ii),f_ond',P1xi(:,ii),'filled',Color=Color{kk},DisplayName=damp{kk}); hold on
        % xlim(ax(110+ii),[0 10])
        % title(ax(110+ii),['Revolution=' num2str(revN)])
        % legend(ax(110+ii))
        
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