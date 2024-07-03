% Author: Alejandro Alvaro, 2023-2024
%
% Sensitivity analysis on the different parameters included in the formulation
% of the non-linear ``elastomeric'' damper constitutive law.
%
% The figures correspond to those containe in Chapter 5.
close all; clear all;

%% Equivalent viscous damping
omega = 2*pi;                   %[rad/s], Angular frequency
t = linspace(0,2*pi/omega,101); % Time vector, [0, T]
A = 1*pi/180;                  % [rad], Oscilation amplitude
% NL viscous law =: f*tanh(alpha*eps_prime)
f_0 = 1000;                 % Amplitude intial value
alpha = 10;                 % Multiplier
%%
% Define a function that calculates C_eq for a given value of f
calculate_C_eq = @(f) calculate_C_eq_given_f(f, omega, A, alpha, t);

% Target value of C_eq
target_C_eq = 4067.5; % Set your desired value here

% Use fminsearch to find the value of f that minimizes the difference
f_optimal = fminsearch(@(f) abs(calculate_C_eq(f) - target_C_eq), f_0);

% Display the optimal value of f
disp(['Optimal value of f: ', num2str(f_optimal)]);

% Calculate C_eq for the optimal value of f
C_eq_optimal = calculate_C_eq(f_optimal);

% Display the calculated C_eq for the optimal value of f
disp(['Calculated C_eq for optimal f: ', num2str(C_eq_optimal)]);

%% Compare with sign function
theta = A*sin(omega*t);         % [rad], orientation
theta_p = omega*A*cos(omega*t); % [rad/s], angular velocity

% NL viscous law =: f*tanh(alpha*eps_prime)
F =@(nu) f_optimal*tanh(alpha*nu);  % Non-linear damping

Ed = trapz(t,F(theta_p).*theta_p);   % Dissipated energy in one cycle
C_eq = Ed/(pi*omega*A^2);   % Equivalent viscous damping


% NL viscous law =: f*sign(eps_prime)
F_s = f_optimal*sign(theta_p);

figure
plot(t,F(theta_p),'DisplayName','N.L $\tanh$'); hold on
plot(t,F_s,'DisplayName','N.L $\mathrm{sgn}$')
legend

%%
% Define the function to calculate C_eq given a value of f
function C_eq = calculate_C_eq_given_f(f, omega, A, alpha, t)
    thetap = omega * A * cos(omega * t);
    F = @(nu) f * tanh(alpha * nu);
    Ed = trapz(t, F(thetap) .* thetap);
    C_eq = Ed / (pi * omega * A^2);
end
