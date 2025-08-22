clearvars, clc, close all 

% -------------------------------------------------------------------------
% Spacecraft Guidance & Navigation (2024/2025)
% Assignment #1, Exercise 3 - Continuous Guidance
% Author: Daniele Paternoster
% -------------------------------------------------------------------------

%% Data of the problem

k1 = 1e-5;
k2 = 1e-4;
Re = 6378.1366;
rho0 = Re + 750;
DU = 7178.1366;
h_i = 800;
h_f = 1e3;
dinc = 0.75;
i2 = deg2rad(0.75);
mu = 398600.435;
Tmax = 3;
m0 = 1e3;
g0 = 9.81;
Isp = 3120;

par.k1 = k1;
par.k2 = k2;
par.rho0 = rho0/DU; %adimensionalized for the debris density calculation

%% 1 - Debris function plot

% Data for the initial and final orbits
r_i = Re + h_i;
r_f = Re + h_f;

% Evaluate debris density function
rho = linspace(r_i-100, r_f+100, 1e3);
q = debris_density(rho/DU, par);

% Plot the debris density function
fontOptions(12)
figure;
plot(rho, q, 'LineWidth',2);
xlabel("$\rho$ [km]"); ylabel("$q$ [1/$DU^3$]");
grid on;
ylim([0 0.104])

% --- Calculate initial and final ORBITAL state - Dimensional units.

% Initialize the two ORBITAL states
x_i = zeros(6,1);
x_f = x_i;

% Initial orbital state
x_i(1) = r_i;
x_i(5) = sqrt(mu/r_i);

% Final orbital state
x_f(1) = r_f;
v_f = sqrt(mu/r_f);
x_f(5) = v_f*cos(i2);
x_f(6) = v_f*sin(i2);

%% Adimensionalization of parameters

m0 = 1e3;
MU = m0;
TU = sqrt(r_i^3/mu);
VU = sqrt(mu/r_i);

% Admiensionalize Thrust, Specific Impulse and Acceleration
par.Tmax = Tmax / (MU*1e3*DU/TU^2);
par.Isp = Isp / TU;
par.g0 = g0 / (1e3*DU/TU^2);

%% 3 - Solve the zero-finding problem of the shooting function

% Adimensionalize full S/C initial state
x_i(1) = x_i(1)/DU;
x_i(5) = x_i(5)/VU;
x_i(7) = m0 / MU;

% Adimensionalize full S/C final state
x_f(1) = x_f(1)/DU;
x_f(5:6) = x_f(5:6)/VU;

% --------------------------CREATE INITIAL GUESSES-------------------------
% Create initial guesses for lambda0 (random numbers between -250 and +250)
% and tf (20pi --> 10revolutions of initial orbit).

X0(1:7) = [-250 + 500*rand; -250 + 500*rand; -250 + 500*rand; ...
    -250 + 500*rand; -250 + 500*rand;-250 + 500*rand;250*rand];
X0(8) = 20*pi + (-pi/2 + pi*rand);

% Best solution found (found by trial&error)
sol_best = 1.0e+02*[-2.149812206879365;
                    -0.103658727446721;
                     0.008855677607339;
                    -0.103929207435009;
                    -2.146104511538979;
                    -1.129453556481870;
                     0.025964491623902;
                     0.644801062443903;
                    ];

% Options used for the found solution
optionsFsolve = optimoptions("fsolve", "Display","iter-detailed", ...
    'FunctionTolerance',1e-7,  'Algorithm','levenberg-marquardt',...
    'MaxFunctionEvaluations',2e4,'MaxIterations', ...
    1e3, 'FiniteDifferenceType','central');

% Zero-finding problem definition and solving procedure.
[sol, fval, exitflag] = fsolve(@(X) objective_fun(X, x_i, par, x_f), sol_best, optionsFsolve);

%% Post-processing on the solution for  T=3N

% Propagate the OCP solution
options = odeset("RelTol", 1e-12, "AbsTol", 1e-12);
[tt, XX] = ode78(@(t,x) dynamics_rhs(t,x,par), [0 sol(end)], [x_i; sol(1:end-1)], options);

fprintf('----------------Solution for T_max = 3.000N-----------------\n')

% Print time of flight
fprintf('Final time: %.4f [mins] \n',sol(end)*TU/60);

% Print final mass
fprintf('Final mass: %.4f [kg] \n', XX(end,7)*MU);

% Print the IC on lambda0
fprintf('Initial co-state [-]: [')
fprintf('%.4f ', sol(1:end-1));
fprintf(']\n');

% Calculate the error on the position and velocity from the solver solution
% using fval
fprintf('Norm of position error: %e [km] \n',norm(fval(1:3)) * DU)
fprintf('Norm of velocity error: %e [m/s] \n',norm(fval(4:6)) * VU * 1e3)

% Calculate the Hamiltonian on the OCP solution - 
HH = hamiltonian(XX,par);
figure;
plot(tt,HH)
grid on;
xlabel('Time [TU]');
ylabel('$H(t) [-]$')

% --- Plot of primer vector in NTW frame ---

%Retrieve the velocity-associated co-state
lambda_v = XX(:,11:13)';

% Initialize the alfa vector in the NTW frame
alfa_NTW = zeros(size(lambda_v));

% For each time instant --> calculate alfa
for k = 1:length(lambda_v)

    % calculate the ECI2NTW position vector transformation matrix
    r = XX(k,1:3)'; 
    v = XX(k,4:6)';
    h = cross(r,v);

    N = h/norm(h); %angular momentum (x)
    T = v/norm(v); %tangential to velocity (y)
    W = cross(N,T)/vecnorm(cross(N,T)); %inside direction (z)
    R_ECI2NTW = [N,T,W]';

    % Convert the primer vector from ECI to NTW --> to get alfa_NTW
    alfa_NTW(:,k) = R_ECI2NTW*(-lambda_v(:,k)/norm(lambda_v(:,k)));
end

% Convert time of thrust in hours
tt_hours = tt*TU/3600; 

% Build the figure of the graphs for alfa NTW components
figure('Name','Primer vector components in NTW frame')

% Set figure position
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.3 0.60 0.60]);

% Set font size and interpreter
fontOptions(14);

t = tiledlayout(3,1);

title(t, '$\mathbf{\alpha}$: NTW frame - $T_{max}$ = 3.000N ','Interpreter','latex','Fontsize', 14)

nexttile
plot(tt_hours,alfa_NTW(1,:),'LineWidth',2)
ylabel("$\alpha_N$ [-]");
grid on;
ylim([-1.05 1.05]);

nexttile
plot(tt_hours,alfa_NTW(2,:),'LineWidth',2);
ylabel("$\alpha_T$ [-]");
grid on;
ylim([-0.05 1.1])

nexttile
plot(tt_hours,alfa_NTW(3,:),'LineWidth',2);
xlabel("$\Delta t_{thrust}$ [h]"); 
ylabel("$\alpha_W$ [-]");
grid on;
ylim([-0.11 0.01])


kep = NaN(6,length(XX));
for k = 1:length(lambda_v)
    kep(:,k) = car2kep(DU*XX(k,1:3), VU*XX(k,4:6), mu);
end


%% Figure - transfer for T = 3.000N

figure;
sgtitle('Orbital transfer for T = 3.000N')
plot3(XX(:,1),XX(:,2),XX(:,3), 'LineWidth', 1.3);
xlabel("$x$ [DU]"); ylabel("$y$ [DU]"); zlabel("$z$ [DU]");
hold on;
grid on;
plot3(XX(1,1), XX(1,2), XX(1,3), 'Marker','x','LineWidth',1.3,'Color','g');
plot3(XX(end,1), XX(end,2), XX(end,3), 'Marker','x','LineWidth',1.3,'Color','k');
plot3(0,0,0,'Marker','o','MarkerSize',5,'LineWidth',5,'Color',[0.3010 0.7450 0.9330])
% x-axis

% Radius increase plot

figure;
sgtitle('Radius of transfer orbit - T = 3.000N')
plot(tt_hours, vecnorm(XX(:,1:3),2,2)*DU ,'LineWidth',2)
xlabel('$\Delta t$ [h]'); ylabel('$||\mathbf{r}||$ [km]');
grid on;
xlim([ 0 17.3])

%% Point 5: new solution for lower thrust

% Compute the Thrust span from 3N to 2.860N (10 elements)
TTadim_span = linspace(3, 2.860, 10)./(MU*1e3*DU/TU^2);

% Initialize fsolve output for the n cases of thrust
sol_NC = [];
fval_NC = [];
exitflag_NC = [];

% Define the first guess of lambda and tf 
guess_NC = sol; 

% Loop over the T_max to perform the numerical continuation of solutions
% k starts from 2 since K=1 is the solution (already found) with T = 3N
for k = 2:length(TTadim_span)

    % Update the T_max parameter
    par.Tmax = TTadim_span(k);

    % Solve the zero-finding problem with new T and IC as solution of
    % previous T.
    optionsFsolve =  optimoptions("fsolve", "Display","iter-detailed",... 
    'FunctionTolerance',1e-7, 'OptimalityTolerance', 1e-6, ... 
    'MaxFunctionEvaluations',2e4,'MaxIterations', 1e3, ...
    'FiniteDifferenceType','central');
    
    [sol_temp, fval_temp, exitflag_temp] = fsolve(@(X) objective_fun(X, x_i, par, x_f), guess_NC, optionsFsolve);

    % Update the guess for the next thrust value with the solution of this
    % step
    guess_NC = sol;

    % Memorization of the solution values (state, value of 0, exit flag)
    sol_NC  = [sol_NC, sol_temp];
    fval_NC = [fval_NC, fval_temp];
    exitflag_NC = [exitflag_NC, exitflag_temp];
end

clear sol_temp fval_temp exitflag_temp;

%% Post-processing on the solution for  T=2.86N

% Propagate the OCP solution
options = odeset("RelTol", 1e-12, "AbsTol", 1e-12);
[tt286, XX286] = ode78(@(t,x) dynamics_rhs(t,x,par), [0 sol_NC(end,end)], [x_i; sol_NC(1:end-1,end)], options);

fprintf('----------------Solution for T_max = 2.860N-----------------\n')

% Print time of flight
fprintf('Final time: %.4f [mins] \n',sol_NC(end,end)*TU/60);

% Print final mass
fprintf('Final mass: %.4f [kg] \n', XX286(end,7)*MU);

% Print the IC on lambda0
fprintf('Initial co-state [-]: [')
fprintf('%.4f ', sol_NC(1:end-1,end));
fprintf(']\n');

% Calculate the error on the position and velocity from the solver solution
% using fval
fprintf('Norm of position error: %e [km] \n',norm(fval_NC(1:3,end)) * DU)
fprintf('Norm of velocity error: %e [m/s] \n',norm(fval_NC(4:6,end)) * VU * 1e3)

% Calculate the Hamiltonian on the OCP solution - 
HH286 = hamiltonian(XX286,par);
figure;
plot(tt,HH)
grid on;
xlabel('Time [TU]');
ylabel('$H(t) [-]$')

% --- Plot of primer vector in NTW frame ---

%Retrieve the velocity-associated co-state
lambda_v286 = XX286(:,11:13)';

% Initialize the alfa vector in the NTW frame
alfa_NTW286 = zeros(size(lambda_v286));

% For each time instant --> calculate alfa
for k = 1:length(lambda_v286)

    % calculate the ECI2NTW position vector transformation matrix
    r = XX286(k,1:3)'; 
    v = XX286(k,4:6)';
    h = cross(r,v);

    N = h/norm(h); %angular momentum (x)
    T = v/norm(v); %tangential to velocity (y)
    W = cross(N,T)/vecnorm(cross(N,T)); %inside direction (z)
    R_ECI2NTW = [N,T,W]';

    % Convert the primer vector from ECI to NTW --> to get alfa_NTW
    alfa_NTW286(:,k) = R_ECI2NTW*(-lambda_v286(:,k)/norm(lambda_v286(:,k)));
end

% Convert time of thrust in hours
tt_hours286 = tt286*TU/3600; 

% Build the figure of the graphs for alfa NTW components
figure('Name','Primer vector components in NTW frame')

% Set figure position
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.3 0.60 0.60]);

% Set font size and interpreter
fontOptions(14);

t = tiledlayout(3,1);

title(t, '$\mathbf{\alpha}$: NTW frame - $T_{max}$ = 2.860N','Interpreter','latex','Fontsize', 14)

nexttile
plot(tt_hours286,alfa_NTW286(1,:),'LineWidth',2)
ylabel("$\alpha_N$ [-]");
grid on;
ylim([-1.05 1.05]);

nexttile
plot(tt_hours286,alfa_NTW286(2,:),'LineWidth',2);
ylabel("$\alpha_T$ [-]");
grid on;
ylim([-1.1 1.1])

nexttile
plot(tt_hours286,alfa_NTW286(3,:),'LineWidth',2);
xlabel("$\Delta t_{thrust}$ [h]"); 
ylabel("$\alpha_W$ [-]");
grid on;
ylim([-0.02 0.01])

kep286 = NaN(6,length(XX286));
for k = 1:length(lambda_v286)
    kep286(:,k) = car2kep(DU*XX286(k,1:3), VU*XX286(k,4:6), mu);
end

radii286 = vecnorm(DU*XX286(:,1:3),2,2);

%% Figure - transfer for T = 2.860N

figure;
sgtitle('Orbital transfer for T = 2.860N')
plot3(XX286(:,1),XX286(:,2),XX286(:,3), 'LineWidth', 1.3);
xlabel("$x$ [DU]"); ylabel("$y$ [DU]"); zlabel("$z$ [DU]");
hold on;
grid on;
plot3(XX286(1,1), XX286(1,2), XX286(1,3), 'Marker','x','LineWidth',1.3,'Color','g');
plot3(XX286(end,1), XX286(end,2), XX286(end,3), 'Marker','x','LineWidth',1.3,'Color','k');
plot3(0,0,0,'Marker','o','MarkerSize',5,'LineWidth',5,'Color',[0.3010 0.7450 0.9330])

% Radius increase plot

figure;
sgtitle('Radius of transfer orbit - T = 2.860N')
plot(tt_hours286, vecnorm(XX286(:,1:3),2,2)*DU ,'LineWidth',2)
xlabel('$\Delta t$ [h]'); ylabel('$||\mathbf{r}||$ [km]');
grid on;
xlim([ 0 17.3])

%% Functions

function q = debris_density(rho, par)

% -------------------------------------------------------------------------
% This function allows to calculate the the debris density in adimensional
% units. 
% INPUT:
% - rho          [n]: array of adimensional radii [DU]
% - par          [-]: struct of parmeters for the function 
%
% OUTPUT:
% - q            [n]: debris denisty value [1/DU^3]
% 
% -------------------------------------------------------------------------
   
k1 = par.k1; k2 = par.k2; rho0 = par.rho0;
q = k1 ./ (k2 + (rho - rho0).^2);

end

function F = dynamics_rhs(~, x, par)
    
% -------------------------------------------------------------------------
% This function allows to calculate the rhs of state/costate dynamics
% coming from the E/L necessary condition of the OCP with respect to
% minimization of the debris impacts. All the units are a-dimensional (the
% characteristic units are in the par struct).
% INPUT:
% - x            [14x1]: full state and co-state vector
% - par          [-]: struct of parmeters for the (constants of the SC)
%
% OUTPUT:
% - F            [14x1]: evaluated RHS of the state/co-state dynamics
% 
% -------------------------------------------------------------------------
% DEPENDENCIES:
% - none
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% state 
r = x(1:3);
rnorm = norm(r);
v = x(4:6);
m = x(7);

% co-state
lambda_r = x(8:10);
lambda_v = x(11:13);
lambda_v_norm = norm(lambda_v);

% Adimensional parameters for the problem of interest
Tmax = par.Tmax;
g0 = par.g0;
Isp = par.Isp;
r0 = par.rho0;
k1 = par.k1;
k2 = par.k2;

% Term of the first component of co-state dynamics relative to debris
% density function
g1 = 2*k1*(rnorm - r0)/(k2 + (rnorm-r0)^2)^2;

% rhs for both dynamics of state (f) and co-state (g)
f = [v; 
     -r/rnorm^3 - (Tmax/m)*lambda_v/lambda_v_norm;
     -Tmax/(Isp*g0)];

g = [ g1*r/rnorm - (3/(rnorm^5))*dot(r,lambda_v)*r + lambda_v/(rnorm^3);
    - lambda_r;
    - lambda_v_norm*Tmax/m^2];

% Build full RHS (state/co-state)
F = [f;g];

end

function H = hamiltonian(X,par)

% -------------------------------------------------------------------------
% This function allows to calculate the Hamiltonian of a given solution
% x(t) for the OCP which minimizes the debris impacts.
% 
% INPUT:
% - X            [nx14]: matrix expressing x(t) for both state and costate
% - par          [-]: struct of parmeters for the function 
%
% OUTPUT:
% - H            [nx1]: hamiltonian H(t)
% 
% -------------------------------------------------------------------------
% DEPENDENCIES:
% - debris_density(r_normalized, par); dynamics_rhs(t,x,par)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Retrieve orbital state from X --> r and v
r = X(:,1:3);
normr = vecnorm(r,2,2);
v = X(:,4:6);

% Calculate the debris density
l = debris_density(normr,par);

% Retrieve the full co-state lambda(t)
lambda = X(:,8:end);

% Calculate the second term of the Hamiltonian definition for each time
% instant
for k = 1:length(X)
     F(:,k) = dynamics_rhs(0, X(k,:)',par);
     f(:,k) = F(1:7,k);
     h2(k) = dot(lambda(k,:), f(:,k));
end

% Calculate the Hamiltonian with its definition
H = l + h2';

end

function Hf = hamiltonian_final(Xf, par)
% -------------------------------------------------------------------------
% This function allows to calculate the Hamiltonian of a given solution
% at final time t_f for agiven solution x of the OCP which minimizes the 
% debris impacts.
% 
% INPUT:
% - Xf           [1x14]: state and costate at t_f
% - par          [-]: struct of parmeters for the function 
%
% OUTPUT:
% - Hf            [1]: hamiltonian H at t_f
% 
% -------------------------------------------------------------------------
% DEPENDENCIES:
% - debris_density(r_normalized, par); dynamics_rhs(t,x,par)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------
   
% Caclulate the RHS of the state/costate dynamics at final state Xf
F = dynamics_rhs(0, Xf, par);

% Retrieve the RHS only for the state components
f = F(1:7);

% Calculate the Hamiltonian by its definition.
Hf = debris_density(norm(Xf(1:3)), par) + Xf(8:14)'*f;

end

function f = objective_fun(X, X0p, par, xf)
% -------------------------------------------------------------------------
% This function allows to calculate the shooting function (objective
% function) for the solution of the TPBVP defined by the OCP which 
% minimizes the debris impacts.
% 
% INPUT:
% - X            [8x1]: Variables to find: initial costate and t_f (independent variables)
% - X0p          [7x1]: Physical state of the S/C
% - par          [-]: struct of parmeters for the function 
% - xf           [6x1]: Final desired orbital state
% OUTPUT:
% - f            [7x1]: value of the shooting function
% 
% -------------------------------------------------------------------------
% DEPENDENCIES:
% -hamiltonian_final(x,par) ; dynamics_rhs(t,x,par)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Optimization variable is X: it contains lambda0 and final time tf
lambda0 = X(1:7);
tfin = X(8);

% Initial condition build up (physical initial state and variables from
% co-state at t0)
XX0 = [X0p; lambda0];

% Propagation of the full state&co-state coupled dynamics from 0 to X(8) = t_f
options = odeset("RelTol", 1e-12, "AbsTol", 1e-12);
[~, XX] = ode78(@(t,x) dynamics_rhs(t,x,par), [0 tfin], XX0, options);

% Desired final orbital state
rf = xf(1:3); 
vf = xf(4:6);

% Transpose output of ode78 and calculate obtained Hamiltonian at final
% time
XXf = XX(end,:)';
Hf = hamiltonian_final(XXf, par);

% Define the objective function as: obtained quantity - desired quantity
% For mass-associated co-state and Hamiltonian at final time the
% desired quantities are zero
f = [XX(end,1:3)' - rf; XX(end,4:6)' - vf; XX(end,14); Hf];

end

function fontOptions(DefaultFontSize)
% Function that sets the fontsize for the axes of plot. and latex
% interpreter for axes labels and legend.
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaulttextinterpreter', 'latex');
    set(0, 'defaultAxesFontSize', DefaultFontSize)
end

