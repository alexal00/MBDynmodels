% Author: Alejandro Alvaro, 2023-2024
% This script segment executes MBDyn and plots the phase-space
% corresponding phase-space plot for the conditions investigated.
%
% The embedding dimension was assumed to be 6x2 = 12, only including the
% dof related to the lead-lag angles of the blades and in-plane support/hub
% movements.
for kk=1:length(damp)

    index = 1;
    for om = omega
    disp(['Omega is at:' num2str(om)])
    variables =['"' 'const real OMEGA_100 = ' num2str(om) ...
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
        '; const real x_0 = ' num2str(x_0) ...
        '; const real y_0 = ' num2str(y_0) ...
        '"'];
    setenv('MBDYNVARS', variables);
   
    % Execution of MBDyn and PP
    
    
        % Preprocessing of file names and testcase chosen
        if sw == 1
            swash = 'sw';
        else
            swash = 'nsw';
        end
    
        f_in = [pref_in folder model damp{kk} '_Nb' int2str(Nb) '_' swash cl '_ODI.mbd'];
    
        % Outfiles
        if b2hdamp==0
            suff = [damp{kk} 'ud'];
        else
            suff = [damp{kk} int2str(Nb) int2str(opt)];
        end
    
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

        % Truncate time series to instants after the perturbation is applied
        idx = psi_nd>=41;
        psi_nd = psi_nd(idx);
        % Extract information about the x and y displacements
        % # 4 hub motions
        AIRFRAME_X = 100; 
        AIRFRAME_Y = 200;
        hub_x = ncread(fn, ['node.struct.',int2str(AIRFRAME_X),'.X'],[1, 1], [1, Inf]);
        hub_x = hub_x(idx);
        hub_y = ncread(fn, ['node.struct.',int2str(AIRFRAME_Y),'.X'],[2, 1], [1, Inf]);
        hub_y = hub_y(idx);
        v_x = ncread(fn,['node.struct.' num2str(AIRFRAME_X) '.XP'],[1,1],[1,Inf]);
        v_x = v_x(idx);
        v_y = ncread(fn,['node.struct.' num2str(AIRFRAME_Y) '.XP'],[2,1],[1,Inf]);
        v_y = v_y(idx);

        % Extract lead-lag angles and azimut for each blade for MBC
        % transformation
        % Read psi vector from the hub joint
        psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
        %Allocate vectors for angular position and velocity
        xi = zeros(Nb,length(time(idx)));
        xid = zeros(Nb,length(time(idx)));
        for ii = 1:Nb
            blade = int2str(10000+1000*ii + 30);
            % SIGN CRITERIA: In MBDyn the lead-lag is positive when leading,
            % while normally it is considered the other way around
            auxxi = -ncread(fn, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);
            xi(ii,:) = auxxi(idx);
            auxxid = -ncread(fn, ['elem.joint.',blade,'.Omega'],[3, 1], [1, Inf]);
            xid(ii,:) = auxxid(idx);
        end

        if find(omega_plot==index)
        nm = ['$\Omega$=' num2str(round(om*60/(2*pi),0))];
        fig(100+kk+tc) = figure(100+kk+tc);
        fig(100+kk+tc).Name = ['PS3D' damp{kk} suffix];

        ax(1) = subplot(3,2,1);
        plot3(ax(1),hub_x,v_x,psi_nd,'Color',Color{1},'LineStyle',LineStyle{index},...
            'Marker',Marker{index},'MarkerEdgeColor','k','MarkerIndices',1:nT:length(psi_nd),'DisplayName',nm)
        xlabel(ax(1),lblx{1})
        ylabel(ax(1),lbly{1})
        zlabel(ax(1),'revolutions')
        hold(ax(1),"on")
        legend(ax(1))

        ax(1) = subplot(3,2,2);
        plot3(ax(1),hub_y,v_y,psi_nd,'Color',Color{2},'LineStyle',LineStyle{index},...
            'Marker',Marker{index},'MarkerEdgeColor','k','MarkerIndices',1:nT:length(psi_nd),'DisplayName',nm)
        xlabel(ax(1),lblx{2})
        ylabel(ax(1),lbly{2})
        zlabel(ax(1),'revolutions')
        hold(ax(1),"on")
        legend(ax(1))

        for ii = 1 : Nb
            ax(1) = subplot(3,2,2+ii);
            plot3(ax(1),xi(ii,:),xid(ii,:),psi_nd,'Color',Color{2+ii},'LineStyle',LineStyle{index},...
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
        plot(hub_x,v_x,'Color',Color{1},'Marker','o','MarkerEdgeColor','k','MarkerIndices',1:nT:length(time(idx)))
        hold(ax(1),"on")
        grid(ax(1),"on")
        xlabel(ax(1),'$x$ [m]')
        ylabel(ax(1),'$\dot{x}$ [m/s]')

        % y vs v_y
        ax(1) = nexttile(tlo(fid));
        plot(hub_y,v_y,'Color',Color{2},'Marker','o','MarkerEdgeColor','k','MarkerIndices',1:nT:length(time(idx)))
        hold(ax(1),"on")
        grid(ax(1),"on")
        xlabel(ax(1),'$y$ [m]')
        ylabel(ax(1),'$\dot{y}$ [m/s]')

        % x vs v_x
        for ii = 1 : Nb
            ax(1) = nexttile(tlo(fid));
            plot(xi(ii,:),xid(ii,:),'Color',Color{2+ii},'Marker','o','MarkerEdgeColor','k','MarkerIndices',1:nT:length(time(idx)))
            hold(ax(1),"on")
            grid(ax(1),"on")
            xlabel(ax(1),['$\xi_' int2str(ii) '$ [rad]'])
            ylabel(ax(1),['$\dot{\xi}_' int2str(ii) '$ [rad/s]'])
        end
        fid = fid+1;


        end
        fig(fid) = figure(fid);
        fig(fid).Name = 'Mxifailed'; fid = fid+1;
        if strcmp(damp(kk),'ib')
            fxi = ncread(fn,'elem.joint.14031.f',[1, 1], [1, Inf]);
        else
            fxi = ncread(fn,'elem.joint.14031.m',[3, 1], [1, Inf]);
        end
        plot(psi_nd,fxi(idx));
        xlabel('revolutions')
        ylabel('Failed damper loads')
        grid on
        % PSmat = [hub_x;v_x;hub_y;v_y;xi;xid];
        % MLCE(index,kk) = calculate_MLCE(PSmat,100,10,15,2,fs);
        index = index+1;
    end
end
