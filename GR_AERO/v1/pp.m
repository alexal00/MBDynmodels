% Author: Alejandro Alvaro, 2023-2024
% generic post-processing for the main script main_realGRsensitivityanalysis
%
% CAVEAT: This function/script is required to execute the main script DO NOT delete
% it and be careful with any modifications
%% Execution of MBDyn and PP
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
        % disp('executing MBDyn...');
    
    
        [rc, errmsg] = system(['wsl ' 'MBDYNVARS=' variables ' ' shfile  f_in ' -o ' f_out]);
    
        % [rc, errmsg] = system(['./mbdyn.sh ' f_in ' -o ' fn_base ' > ' fn_base '.txt']);
        if (rc ~= 0)
            error(errmsg);
        end
        % disp('   ... done');
    end
    
    % filename
    fn = [f_out '.nc'];
    % disp(sprintf('reading output from file ''%s''', fn));
    
    
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
end