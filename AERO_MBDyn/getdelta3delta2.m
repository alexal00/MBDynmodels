function sol=getdelta3delta2(params)
    % Starting values for the optimization
    X0          = [0 1.];
    options     = optimset(...
              'MaxFunEvals',300*length(X0),...
              'TolX',1e-4,...
              'TolFun',1e-4);
    myfun = @(X)getErrorfun(X,params);
    % myfun2 = @(X)getErrorfun2(X,params);
    % lb = [-inf -inf];
    % delta2max = params.nuxi^2*(params.nubeta^2+params.lock/8*tan(pi/4))/...
    %     (params.nubeta^2*params.lock*params.lamdda0/6);
    % up = [inf delta2max];
    [X,fval,FLAG]       = fminsearch(myfun,X0,options);
    % [X2,fval2,FLAG]       = fminsearch(myfun2,X0,options);
    % if fval<=fval2
    sol = X;
    err = fval;
    % else
    %     sol = X2;
    %     err = fval2;
    % 
    % end
    format long
    % [X,fval,~]       = fmincon(myfun,X0,[],[],[],[],lb,up,[],options);
    disp(strcat('delta3=',num2str(atan(sol(1)),'%.5f')))
    disp(strcat('delta2=',num2str(sol(2),'%.5f')))
    disp(strcat('Rel. error(%)=',num2str(err)))
    format short
end
function e = getErrorfun(x,params)
beta0 = params.beta0;
xi0 = params.xi0;
nubeta = params.nubeta;
nuxi = params.nuxi;
lock = params.lock;
lambda0 = params.lamdda0;
fb = params.fb;
fx = params.fx;
% delta3 = params.delta3;
Kpt = x(1);

B = [fb;fx];
A_r = [nubeta^2+lock*Kpt/8 -lock/8*x(2);...
           lock*lambda0*Kpt/6 nuxi^2-lock/6*lambda0*x(2)];

e = sum((B-A_r*[beta0;xi0]).^2);
end

function e = getErrorfun2(x,params)
input_t0 = params.INPUT_THETA_0;
beta0 = params.beta0;
xi0 = params.xi0;
theta0 = params.theta0;
delta3 = params.delta3;
e = abs(theta0-(input_t0-tan(delta3)*beta0+x*xi0))/abs(theta0);
end
