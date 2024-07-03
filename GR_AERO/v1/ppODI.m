% Author: Alejandro Alvaro, 2023-2024
% generic post-processing for the MBDyn output files obtained from the
% execution of main_realGRODICzetavsomega
%
% CAVEAT: This function/script is required to execute the main script DO NOT delete
% it and be careful with any modifications
%% Execution of MBDyn an partial PP
% Allocate damping ratio
zeta_col = zeros(length(damp),length(dof_col));
zeta_cyc = zeros(length(damp),length(dof_cyc));
% omega = 26;
% Allocate vectors for further plotting
zetacol_om = zeros(length(damp),length(dof_col),length(omega));
zetacyc_om = zeros(length(damp),length(dof_cyc),length(omega));

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
    
    for kk=1:length(damp)
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
    
        fignames = [damp{kk} int2str(Nb) num2str(round(om,0))];
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
        psi_nd = time*om/(2*pi);
    
        % Pp for fft and truncation of signals
        nT = 256;                           % Steps per period
        revN = fix(length(psi_nd)/nT) - 1;  % Number of revs -1
    
    
        % Extract lead-lag angles and azimut for each blade for MBC
        % transformation
        psi = ncread(fn, ['elem.joint.' num2str(10000) '.Phi'], [3, 1], [1, Inf]);
        psi_b = zeros(Nb,length(time));
        xi = zeros(Nb,length(time));
        for ii = 1:Nb
            psi_b(ii,:) = psi+ii*dpsi;
            blade = int2str(10000+1000*ii + 30);
            % SIGN CRITERIA: In MBDyn the lead-lag is positive when leading,
            % while normally it is considered the other way around
            xi(ii,:) = -ncread(fn, ['elem.joint.',blade,'.Phi'],[3, 1], [1, Inf]);
            % xid(ii,:) = -ncread(fn, ['elem.joint.',blade,'.Omega'],[3, 1], [1, Inf]);
        end
        % blade = int2str(10000+1000*Nb + 30);
        % Store some indexes to perform the logarithmic decrement estimate
        % Perturbation time is from 40*T_REV to 40*T_REV+dt
        idx = 1:nT*12;
        T = 2*pi/om;       
        xi_NR_MB = MBC(psi_b,xi);
        qNR = NRdof(Nb);
        % Logarithmic decrement
        for jj = 1:length(dof_col)
            zeta_col(kk,jj) = logarithmic_decrement(psi_nd(idx),xi_NR_MB(dof_col(jj),idx),'PlotFlag',false,'MinPeakDistance',0.5);
        end
        for jj = 1:length(dof_cyc)
            zeta_cyc(kk,jj) = logarithmic_decrement(psi_nd(idx),xi_NR_MB(dof_cyc(jj),idx),'PlotFlag',false,'MinPeakDistance',0.5);
        end
        
    end

    zetacol_om(:,:,index) = zeta_col;
    zetacyc_om(:,:,index) = zeta_cyc;
    index = index+1;
end

%% Global figure declaration
% Evolution of the damping ratios with omega
fig(fid) = figure(fid);
fig(fid).Name = ['ZetavsOmega' suffix]; fid = fid +1;
% Define variables
numRows = 2; % Number of rows (1, 2, or 3)
numCols = size(zetacol_om, 1); % Maximum columns between the two arrays
% Loop through each row
for ii = 1:numCols
    subplot(numRows, numCols, ii);
    hold on
    for jj = 1:size(zetacol_om, 2)
        zetacoljj(:,jj) = squeeze(zetacol_om(ii,jj,:));
    end
    plot(omega(1:(end-1))*60/(2*pi), zetacoljj(1:(end-1),:));
    title(damp{ii});
    xlabel('$\Omega$ [rpm]')
    ylabel('$\zeta_{\mathrm{col}}$')
    grid on
    legend(qNR{dof_col});
end
for ii = 1: numCols
    subplot(numRows, numCols, numCols + ii);
    hold on
    for jj = 1:size(zetacyc_om, 2)
        zetacycjj(:,jj) = squeeze(zetacyc_om(ii,jj,:));
    end
    plot(omega(1:(end-1))*60/(2*pi), zetacycjj(1:(end-1),:));
    xlabel('$\Omega$ [rpm]')
    ylabel('$\zeta_{\mathrm{cyc}}$')
    grid on
    legend( qNR{dof_cyc});
end

fig(fid) = figure(fid);
fig(fid).Name = ['ZetavsOmegaBAR' suffix]; fid = fid+1;
% Define variables
numRows = 2; % Number of rows (1, 2, or 3)
numCols = size(zetacol_om, 1); % Maximum columns between the two arrays
% Loop through each row
for ii = 1:numCols
    subplot(numRows, numCols, ii);
    hold on
    for jj = 1:size(zetacol_om, 2)
        zetacoljj(:,jj) = squeeze(zetacol_om(ii,jj,:));
    end
    bar(omega(1:(end-1))*60/(2*pi), zetacoljj(1:(end-1),:),...
            'DisplayName', [damp{ii} ' ' qNR{dof_col(jj)}]);
    title(damp{ii});
    xlabel('$\Omega$ [rpm]')
    ylabel('$\zeta_{\mathrm{col}}$')
    grid on
    legend(qNR{dof_col});
end
for ii = 1: numCols
    subplot(numRows, numCols, numCols + ii);
    hold on
    for jj = 1:size(zetacyc_om, 2)
        zetacycjj(:,jj) = squeeze(zetacyc_om(ii,jj,:));
    end
    bar(omega*60/(2*pi), zetacycjj);
    xlabel('$\Omega$ [rpm]')
    ylabel('$\zeta_{\mathrm{cyc}}$')
    grid on
    legend( qNR{dof_cyc});
end

if any(strcmp(damp, 'i2b')) && any(strcmp(damp, 'ib'))

    fig(fid) = figure(fid);
    fig(fid).Name = ['ibvsi2b' suffix];
    idib = find(strcmp(damp, 'ib'));
    idi2b = find(strcmp(damp, 'i2b'));
    for jj = 1:size(zetacyc_om, 2)
        zetacycrat(:,jj) = squeeze(zetacyc_om(idi2b,jj,:))./squeeze(zetacyc_om(idib,jj,:));
    end
    bar(omega*60/(2*pi),zetacycrat)
    xlabel('$\Omega$ [rpm]')
    ylabel('$\zeta_{\mathrm{cyc}}^{i2b}/\zeta_{\mathrm{cyc}}^{ib}$')
    grid on
    legend( qNR{dof_cyc});
    fid = fid+1;
end