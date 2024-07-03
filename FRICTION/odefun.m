function dydt = odefun(t,y,params)
% odefun Function that contains the ODE systme of equations
%        used as comparison with the corresponding results
%        from the simple MBDYn model
%
% Input:
%       * t : time vector
%       * y : angular displacement
%       * params: structure that contains the required parameters
%                 to define the 1 dof system
% Output:
%       * dydt : system of ODE

    Jb = params.Jb;
    Kx = params.K_xi;
    f = params.f;
    alpha = params.alpha;
    Mz = params.M_z;
    omega = params.omega_f;


    dydt(1,1) = y(2);
    dydt(2,1) = -Kx/Jb*y(1)-f/Jb*tanh(alpha*y(2))+Mz*sin(omega*t);


end