
clc; clearvars; close all; cspice_kclear

% -------------------------------------------------------------------------
% Spacecraft Guidance & Navigation (2024/2025)
% Assignment #1, Exercise 2 - Impulsive Guidance
% Author: Daniele Paternoster
% -------------------------------------------------------------------------
% NOTES: 
% - put the script in the same path of the "kernels" folder with
%   the spice kernels inside. 
% - Multiple Shooting solution can take some minutes to be found at each
%   run if the "stringentOptions" are used. The already found solution is
%   saved in a .mat file (go to line 150 for more details).
% -------------------------------------------------------------------------

%% Upload Spice kernels in kernel pool

cspice_furnsh('kernels\gm_de432.tpc')
cspice_furnsh('kernels\pck00010.tpc');
cspice_furnsh('kernels\de432s.bsp');
cspice_furnsh('kernels\naif0012.tls');

% Extract Moon & Earth Kernel data

% Gravitational parameters for Earth and Moon
muE = cspice_bodvrd('399', 'GM', 1);
muM = cspice_bodvrd('301', 'GM', 1);

% Radii of Earth
radiiE = cspice_bodvrd('399', 'RADII', 3);

% Equatorial radius considered
rE = radiiE(1);

% Radii of Moon
radiiM = cspice_bodvrd('301', 'RADII', 3);
% All radii are equal (flatness = 0)
rM = radiiM(1);

% compute mu parameter for the CRTBP
mu = muM / (muE + muM);

% --- Problem Data --- 

% Dimensional altitudes
hi_d = 167;
hi_a = 100;

% Earth-Moon distance
l_em = 3.84405000e5;

% Omega Earth-Moon [rad/s]
om_em = 2.66186135e-6;

% Adimensional angular velocity of Sun in synodic frame
om_s  = -9.25195985e-1;

% Adimensional Sun - (Earth/Moon) distance
rho = 3.88811143e2;

% Adimensional mass of Sun
m_s = 3.28900541e5;

% --- Scaling Units ---

% Time unit [s]
TU = 4.34256461 * 24 * 3600;

% Distance unit [km]
DU = l_em;

% Velocity unit [km/s]
VU = 1.02454018;

% --- Initial guess parameters --- 
alfa = 0.2*pi; %1.5*pi;
beta = 1.41;
delta = 4; % 7;
ti = 2; %0;

% Struct which contains parameters of the problem
par.mu = mu;
par.rho = rho;
par.om_s = om_s;
par.m_s = m_s;
par.rE = rE;
par.rM = rM;
par.DU = DU;

%% 1) Initial guess & propagation

% Build-up the initial guess in cartesian state (adimensional) from the parameters
ri = (rE + hi_d)/DU;
rf = (rM + hi_a)/DU;
v0 = beta*sqrt((1-mu)/ri);
x0 = [ri*cos(alfa) - mu
           ri*sin(alfa)
     -(v0-ri)*sin(alfa)
     (v0-ri)*cos(alfa)];

% Propagation of the initial guess
options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[tt, xx] = ode78(@(t,x) PBRFBP_rhs(t,x,par), [ti ti + delta], x0, options);

% Plot the propagation in the EMB centered rotating (synodic) frame
twoDplot_rotframe(xx,mu,ri,rf, 'First Guess propagated - @EM rot. frame');

% Plot the propagation in the Earth centered inertial frame
twoDplot_ECinertial(xx, tt, mu, 'First Guess propagated - @Inertial frame')

%% 2) Simple shooting 

% Build the initial guess for the shooting problem (NLP variabile) 
y0 = [x0; ti; ti+delta];

% --- Simple shooting without gradients ---
 
SIMPLE_NOGRAD = simple_shoot(y0, par, ri, rf, false);

% --- Simple shooting with the gradients ---

SIMPLE_WGRAD = simple_shoot(y0, par, ri, rf, true);

%% Post-processing of simple shooting solution

fprintf('---------------SIMPLE SHOOTING noGradients -----------------\n');
fprintf('The found deltaV is: %f [km/s]\n', SIMPLE_NOGRAD.minvalue * VU);
fprintf('The found deltaT is: %f [days]\n', (SIMPLE_NOGRAD.state(end) - SIMPLE_NOGRAD.state(end-1))*TU/(24 * 3600));
fprintf('The max constraint violation is: %e [-]\n', max(SIMPLE_NOGRAD.c_f_at_sol))

fprintf('---------------SIMPLE SHOOTING withGradients -----------------\n');
fprintf('The found deltaV is: %f [km/s]\n', SIMPLE_WGRAD.minvalue * VU);
fprintf('The found deltaT is: %f [days]\n', (SIMPLE_WGRAD.state(end) - SIMPLE_WGRAD.state(end-1))*TU/(24 * 3600));
fprintf('The max constraint violation is: %e [-]\n', max(SIMPLE_WGRAD.c_f_at_sol))

% Simple shooting plots 

% w/out gradient
[tt_SSnograd, xx_SSnograd] = ode78(@(t,x) PBRFBP_rhs(t,x,par), [SIMPLE_NOGRAD.state(5) SIMPLE_NOGRAD.state(6)], SIMPLE_NOGRAD.state(1:4), options);
twoDplot_rotframe(xx_SSnograd,mu,ri, rf, 'Simple Shooting no/grad - @EMB Rot. frame');
twoDplot_ECinertial(xx_SSnograd, tt_SSnograd, mu,'Simple Shooting no/grad - @Earth Inertial frame')

% w/ gradient
[tt_SSwgrad, xx_SSwgrad] = ode78(@(t,x) PBRFBP_rhs(t,x,par), [SIMPLE_WGRAD.state(5) SIMPLE_WGRAD.state(6)], SIMPLE_WGRAD.state(1:4), options);
twoDplot_rotframe(xx_SSwgrad,mu,ri, rf, 'Simple Shooting w/grad - @EMB Rot. frame');
twoDplot_ECinertial(xx_SSwgrad, tt_SSwgrad, mu,'Simple Shooting w/grad - @Earth Inertial frame')


%% 3) multiple shooting with N = 4;

y01 = x0;
tspan = linspace(ti, ti+delta, 4);

% !!! Uncomment stringentOptions to propagate the refined solution (it takes 10 mins)
% and put them in the propagation of line 159 !!!

% stringentOptions = odeset('RelTol',1e-12, 'AbsTol', 1e-12);
[t_0, y00] = ode78(@(t,x) PBRFBP_rhs(t,x,par), tspan, x0); %, stringentOptions); 

y02 = y00(2,:);
y03 = y00(3,:);
y04 = y00(4,:);

yM0 = [y01; y02'; y03'; y04'; ti; ti+delta];

[MULTIPLE_WGRAD] = multiple_shoot(yM0, par, ri, rf);

%% Post processing on Multiple Shoot solution

% Piece-wise propagation of the found solution
tt_MS_span = linspace(MULTIPLE_WGRAD.state(end-1), MULTIPLE_WGRAD.state(end),4);
options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[tt_MS_1, xx_MS_1] = ode78(@(t,x) PBRFBP_rhs(t,x,par), linspace(tt_MS_span(1),tt_MS_span(2),500), MULTIPLE_WGRAD.state(1:4), options);
[tt_MS_2, xx_MS_2] = ode78(@(t,x) PBRFBP_rhs(t,x,par), linspace(tt_MS_span(2),tt_MS_span(3),500), MULTIPLE_WGRAD.state(5:8), options);
[tt_MS_3, xx_MS_3] = ode78(@(t,x) PBRFBP_rhs(t,x,par), linspace(tt_MS_span(3),tt_MS_span(4),500), MULTIPLE_WGRAD.state(9:12), options);

% Reconstruct the solution (three arcs)
xx_Msol1 = [xx_MS_1; xx_MS_2; xx_MS_3];
tt_Msol1 = [tt_MS_1; tt_MS_2; tt_MS_3];

% Indexes where the solution is patched (for future plots)
indxMS = [1 500 1000 1500];

% Plot the multiple shooting solution
twoDplot_rotframe(xx_Msol1,mu, ri, rf, 'Multiple Shooting solution - @EMB rot. frame',indxMS);
twoDplot_ECinertial(xx_Msol1, tt_Msol1, mu, 'Multiple Shooting solution - @EMB rot. frame',indxMS)

%% Refined solution (it is the found solution without the need to run multipleshoot)

% Load the already found solution (it takes around 10 mins to find it using
% stringentTolerances). The solution inside the .mat file is called the:
% MULTIPLE_WGRAD_refIG --> [struct] with fields
% .state --> 18x1 NLP variable as presented in the report
% .minvalue --> adimesional cost [V]
% .extiflag of the solver
% .out --> solver details (from fmincon).

load MSHOOT_refined.mat

% Define the times at which solution is discretized
tt_MSREF_span = linspace(MULTIPLE_WGRAD_refIG.state(end-1), MULTIPLE_WGRAD_refIG.state(end), 4);

% Propagate the three arcs 
[ttREF_Msol1, xxREF_Msol1] = ode78(@(t,x) PBRFBP_rhs(t,x,par), linspace(tt_MSREF_span(1),tt_MSREF_span(2),500), MULTIPLE_WGRAD_refIG.state(1:4), options);
[ttREF_Msol2, xxREF_Msol2] = ode78(@(t,x) PBRFBP_rhs(t,x,par), linspace(tt_MSREF_span(2),tt_MSREF_span(3),500), MULTIPLE_WGRAD_refIG.state(5:8), options);
[ttREF_Msol3, xxREF_Msol3] = ode78(@(t,x) PBRFBP_rhs(t,x,par), linspace(tt_MSREF_span(3),tt_MSREF_span(4),500), MULTIPLE_WGRAD_refIG.state(9:12),options);

% Reconstruct the arcs (patch them)
xxREF_Msol = [xxREF_Msol1; xxREF_Msol2; xxREF_Msol3];
ttREF_Msol = [ttREF_Msol1; ttREF_Msol2; ttREF_Msol3];

% Indexes for plots
indxMS = [1 500 1000 1500];

% Plots of MS solution loaded
twoDplot_rotframe(xxREF_Msol,mu, ri, rf, 'Multiple Shooting refined solution - @EMB rot. frame',indxMS);
twoDplot_ECinertial(xxREF_Msol, ttREF_Msol, mu,'Multiple Shooting refined solution - @EMB rot. frame',indxMS)

% --- MS solutions ---

fprintf('---------------MULTIPLE SHOOTING-----------------\n');
fprintf('The found deltaV is: %f [km/s]\n', MULTIPLE_WGRAD.minvalue * VU);
fprintf('The found deltaT is: %f [days]\n', (MULTIPLE_WGRAD.state(end) - MULTIPLE_WGRAD.state(end-1))*TU/(24 * 3600));
fprintf('The max constraint violation is: %e [-]\n', max(MULTIPLE_WGRAD.psi_eq))

fprintf('---------------REFINED MULTIPLE SHOOTING-----------------\n');
fprintf('The found deltaV is: %f [km/s]\n', MULTIPLE_WGRAD_refIG.minvalue * VU);
fprintf('The found deltaT is: %f [days]\n', (MULTIPLE_WGRAD_refIG.state(end) - MULTIPLE_WGRAD_refIG.state(end-1))*TU/(24 * 3600));
[~, ceq] = psi(MULTIPLE_WGRAD_refIG.state,par,ri,rf);
fprintf('The max constraint violation is: %e [-]\n', max(ceq));


%% 4) - Match with real world condition

fprintf('---------------------MATCH WITH REAL WORLD---------------------\n')

% Find the theta angle of the Sun at the solution
theta_i = wrapTo2Pi(par.om_s* SIMPLE_WGRAD.state(end-1));

% Define the initial guess date
initial_date = '2024-09-28-00:00:00.000 TDB';

% Convert the date to ephemeris time
et0 = cspice_str2et(initial_date);

% Define the approximate periodicity for theta
deltaet = (2*pi)/(-om_s * om_em );

% --- Plot the function to define the search interval ---

% Define the timespan in ET
ett = linspace(et0, et0 + deltaet, 1e3);

% Allocate space for the theta difference 
diff = NaN(size(ett));

% Calculate the difference (removing vertical line at discontinuity)
for k = 1:length(ett)
    diff(k) = thetaECLIP(ett(k), theta_i);
    if k > 1 && diff(k-1)*diff(k) < 0
        diff = [diff(1:k-1), NaN, diff(k:end)];
        ett  = [ett(1:k-1), NaN, ett(k:end)];
    end
end

% Plot dTheta as function of time
figure;

plot(ett, rad2deg(diff), 'LineWidth', 2)

% Create labels with UTC date
et_dates = linspace(ett(1), ett(end), 10);
for k = 1:length(et_dates)
    label_temp = cspice_et2utc(et_dates(k), 'C', 1);
    labels{k} = label_temp(1:11);
end

% Put char labels instead of numbers
xticks(et_dates)
xticklabels(labels);

grid minor; 
xlabel('UTC date'); ylabel('$\Delta \theta \, [^{\circ}]$')

% Use the graphical information to select the interested domain of
% zero-finding problem
% Reasonable interval in ET --> [7.815e8 7.83e8]
fzeroOPT = optimset('Display','Iter');
[et, fval, flag_et, output_dthzero] = fzero(@(et) thetaECLIP(et, theta_i), [7.815e8 7.83e8],fzeroOPT);


% -- Conversion of initial condition from EMB to ECI and physical units--

% Extract adimensional initial condition 
x_SS_grad = (SIMPLE_WGRAD.state)';
% Transform the initial condition from EMB2ECI
[R, V] = EMB2ECI(x_SS_grad(1:2), x_SS_grad(3:4), x_SS_grad(5), par.mu);
% Define the time of propagation
det = (x_SS_grad(end) - x_SS_grad(end-1))*TU;

% Print the epoch of interest (for the next propagation)
fprintf('Starting epoch of propagation: [UTC] %s\n',cspice_et2utc(et, 'ISOC', 14))
fprintf('Ending epoch of propagation  : [UTC] %s\n',cspice_et2utc(et+det, 'ISOC', 14))

% Collect the new state in vector X
X = [R*DU,0,V*VU,0,et,et+det];

%% N-body propagation

% Cell of labels of main attractors (orderer in magnitude)
labels = {'Moon';
          'Sun';
          'Venus';
          'Mars Barycenter';
          'Mercury';
          'Jupiter Barycenter';
          'Saturn Barycenter';
          'Uranus Barycenter';
          'Neptune Barycenter';
          'Pluto Barycenter'};

% Define a new cell sized as labels 
bodies = cell(size(labels));

% Populate the bodies cell with names and GM of each attractor
for i = 1:length(labels)
    bodies{i}.name = labels{i};
    bodies{i}.GM   = cspice_bodvrd(labels{i}, 'GM', 1);
end

% Set options for dimensional n-body propagation
optionsDIM = odeset('RelTol',1e-12, 'AbsTol', [1e-8*ones(3,1); 1e-11*ones(3,1)]);

% N-body propagator wrt ECI frame
[TT, XX] = ode78(@(t,x) nbody_ECI_rhs(t,x,bodies), linspace(X(end-1), X(end), 2e3), X(1:6), optionsDIM);

% Propagation of SS solution with PBRFBP 
[tt_sol1, xx_sol1] = ode78(@(t,x) PBRFBP_rhs(t,x,par), linspace(SIMPLE_WGRAD.state(end-1), SIMPLE_WGRAD.state(end), 2e3) , SIMPLE_WGRAD.state(1:4), options);
% Convert the solution of the PBRFBP from EMB rotating to ECI frame
[R_sol1 ,V_sol1] = EMB2ECI(xx_sol1(:,1:2), xx_sol1(:,3:4), tt_sol1, mu);
% Create the dimensional trajectory by adding units and z components
XX_sol1 = [DU*R_sol1, zeros(size(tt_sol1)), VU*V_sol1, zeros(size(tt_sol1))];

%% 3d qualitative plots - comparison of N-body and PBRFBP propagations

% Retrieve Moon orbit in ECLIPJ2000 frame Earth centered
r_MOON_ECLIP = cspice_spkpos('MOON', TT','ECLIPJ2000','NONE','EARTH');

% Create plots
figure;
sgtitle('N-Body propagation - @Earth:ECLIPJ2000')

% N-body orbit
plot3(XX(:,1), XX(:,2),XX(:,3), 'LineWidth',2, 'Color',[0 0.4470 0.7410]);
hold on; grid on; axis square;
plot3(XX(1,1), XX(1,2),XX(1,3), 'o','MarkerSize',3,'LineWidth',3, 'Color',[0 0.4470 0.7410]);
plot3(XX(end,1), XX(end,2),XX(end,3), 'diamond','MarkerSize',3,'LineWidth',3, 'Color',[0 0.4470 0.7410]);

% PBRFBP orbit
plot3(XX_sol1(:,1), XX_sol1(:,2), XX_sol1(:,3), 'LineWidth',2,'Color',[0.8500 0.3250 0.0980])
plot3(XX_sol1(1,1), XX_sol1(1,2), XX_sol1(1,3),'o', 'MarkerSize', 3, 'LineWidth',3,'Color',[0.8500 0.3250 0.0980])
plot3(XX_sol1(end,1), XX_sol1(end,2), XX_sol1(end,3), 'diamond', 'MarkerSize', 3,'LineWidth',3,'Color',[0.8500 0.3250 0.0980])

% Moon orbit
plot3(r_MOON_ECLIP(1,:), r_MOON_ECLIP(2,:), r_MOON_ECLIP(3,:),'LineStyle','-.','Color', [0.7 0.7 0.7]);
plot3(r_MOON_ECLIP(1,80:200), r_MOON_ECLIP(2,80:200), r_MOON_ECLIP(3,80:200),'Color', 'k','LineWidth',3, 'HandleVisibility','off');
plot3(r_MOON_ECLIP(1,1), r_MOON_ECLIP(2,1), r_MOON_ECLIP(3,1), 'o', 'MarkerSize', 3 ,'LineWidth',3,'Color',[0.4940 0.1840 0.5560])
plot3(r_MOON_ECLIP(1,end), r_MOON_ECLIP(2,end), r_MOON_ECLIP(3,end), 'diamond', 'MarkerSize', 3 ,'LineWidth',3,'Color',[0.4940 0.1840 0.5560])

xlabel('$x$ [km]'); ylabel('$y$ [km]'); zlabel('$z$ [km]');
legend('nBD','nBD-start','nBD-end','PBRFBP','PBRFBP-start','PBRFBP-end','Moon Orbit','Moon Start','Moon End','FontSize', 8)


% --- Quantitative plots ---

% Calculate the difference on the position
err_R_x = XX_sol1(:,1) - XX(:,1);
err_R_y = XX_sol1(:,2) - XX(:,2);
err_R_z = XX_sol1(:,3) - XX(:,3);

% Calculate the difference on the velocity.
err_V_x = XX_sol1(:,4) - XX(:,4);
err_V_y = XX_sol1(:,5) - XX(:,5);
err_V_z = XX_sol1(:,6) - XX(:,6);

% Show quantitative plots
figure
t = tiledlayout(2,1); %, "TileSpacing","compact",'Padding','compact');
TT_d = (TT - TT(1))/(3600*24);

% Position error
nexttile
plot(TT_d, err_R_x, 'LineWidth',1.6);
hold on;
grid on;
plot(TT_d, err_R_y, 'LineWidth',1.6);
plot(TT_d, err_R_z, 'LineWidth',1.6);
legend("$x_{FBP} - x_{nBD}$", "$y_{FBP} - y_{nBD}$", "$z_{FBP} - z_{nBD}$");
ylim([-8 1.4]*1e4);
xlabel('$\Delta t$ [days]');
ylabel('$\epsilon_r$[km]')


% Velocity error
nexttile
plot(TT_d, err_V_x, 'LineWidth',1.6);
hold on
grid on;
plot(TT_d, err_V_y, 'LineWidth',1.6);
plot(TT_d, err_V_z, 'LineWidth',1.6);
legend("$\dot{x}_{FBP} - \dot{x}_{nBD}$", "$\dot{y}_{FBP} - \dot{y}_{nBD}$", "$\dot{z}_{FBP} - \dot{z}_{nBD}$");
xlabel('$\Delta t$ [days]');
ylabel('$\epsilon_v$[km/s]')
ylim([-2.5 1.4])


%% FUNCTIONS

function [dxdt] = nbody_ECI_rhs(t,x,bodies)
% -------------------------------------------------------------------------
% This function allows to calculate the rhs of n-body dynamics centerd 
% in the Earth Centered Inertial frame with directions defined by the
% ECLIPJ2000.
% All the units are dimensional.
% INPUT:
% - t            [1]: time 
% - xx           [6x1]: state vector
% - bodies       [-]: cell containing bodies and GMs for the propagation
% OUTPUT:
% - dxdt         [6x1]: evaluated RHS nbody dynamics
% DEPENDENCIES:
% cspice_bodvrd, cspice_spkpos
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Allocate memory for dxdt
dxdt = NaN(6,1);

% First three components of dxdt are the velocity components
dxdt(1:3) = x(4:6);

% Define r as S/C position and its norm
r = x(1:3);
rnorm = norm(r);

% Retrieve GM for Earth
GM0 = cspice_bodvrd('Earth', 'GM', 1);

% Calculate acceleration of central body (Earth)
dxdt(4:6) = -GM0 * r / rnorm^3;

% Loop over the bodies
for k = 1:size(bodies)

    % Define the position vector of k-th body wrt Earth in ECLIPJ2000
    rho = cspice_spkpos(bodies{k}.name, t, 'ECLIPJ2000', 'NONE', 'EARTH');
    
    % Calculate distance vector d between S/C and k-th body
    d = r - rho;
    dnorm = norm(d);
    
    % Calculate the Battin parameters
    q = dot(r, r-2*rho) / dot(rho,rho);
    f = q*(3+3*q+q^2)/(1 + (1+q)^1.5);

    % Calculate the perturbation terms with Battin parameters
    dxdt(4:6) = dxdt(4:6) - bodies{k}.GM *(r + rho*f) / dnorm^3;
end

end

function [R,V] = EMB2ECI(r,v,t,mu)
% -------------------------------------------------------------------------
% This function allows to convert the full state vector from EMB rotating 
% frame to ECI frame.
% INPUT:
% - r            [nx2]: position vector in rotating EMB frame
% - v            [nx2]: velocity vector in rotating EMB frame
% - t              [n]: time span 
%
% OUTPUT:
% - R            [nx2]: position vector in inertial EC frame
% - V            [nx2]: velocity vector in inertial EC frame
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Define the position and velocity components
x = r(:,1); 
y = r(:,2); 
vx = v(:,1);
vy = v(:,2);

% Apply the transformation (see Topputo,2013)
R(:,1) = (x+mu).*cos(t) - y.*sin(t);
R(:,2) = (x+mu).*sin(t) + y.*cos(t);
V(:,1) = (vx-y).*cos(t) - (vy+x+mu).*sin(t);
V(:,2) = (vx-y).*sin(t) + (vy+x+mu).*cos(t);

end

function diff_th = thetaECLIP(et_val, th_i)
% -------------------------------------------------------------------------
% This function allows to calculate the difference in angles between the
% value th_i and the angle between the Sun-EMB vector projected onto the 
% EM motion plane and Moon-EMB vector (x axis of the EMBrotating frame)
% INPUT:
% - et_val         [1]: ephemeris time at which the angle must be computed
% - th_i           [1]: angle to which subtract th
% OUTPUT:
% - diff_th        [1]: difference between th and th_i
% DEPENDENCIES:
% cspice_spkezr, cspice_spkpos
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Moon position wrt EMB in ECLIPJ200 @et_val
xxM = cspice_spkezr('MOON', et_val, 'ECLIPJ2000', 'NONE', '3');

% Define Moon position and velocity
rM = xxM(1:3);
vM = xxM(4:6);

% Sun position wrt EMB in ECLIPJ2000 @et_val
rsun = cspice_spkpos('SUN', et_val, 'ECLIPJ2000', 'NONE', '3');

% Construct x,y,z axis of the EMB centered rotating frame expressed in
% ECLIPJ2000
xrot = rM / norm(rM);
zrot = cross(rM,vM);
zrot = zrot/norm(zrot);

yrot = cross(zrot, xrot);
yrot = yrot/norm(yrot);

% rotation matrix (EMB2ECLIP)
R = [xrot, yrot, zrot];

% Sun direction expressed in the EMB centered rotating frame (3d)
rsun_EMB = R' * rsun;

% Calculate theta angle between sun direction and x axis (EMB direction)
th = wrapTo2Pi(atan2(rsun_EMB(2), rsun_EMB(1)));

% Calculate the difference between theta values
diff_th = th - th_i;

end

function PBRFBP_rhs = PBRFBP_rhs(t, xx, par)
% -------------------------------------------------------------------------
% This function allows to calculate the rhs of the PBRFBP dynamics. All the 
% units are a-dimensional (the characteristic units are in the par struct).
% INPUT:
% - t            [1]: time 
% - xx           [4x1]: state vector
% - par          [-]: struct of parmeters 
%
% OUTPUT:
% - PBRFBP_rhs   [4x1]: evaluated RHS of the state dynamics
%
% DEPENDENCIES:
% - OM4(t,x(1,2), par)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Calculate gradient of the PBRFBP potential (OM4)    
[~, dOM4] = OM4(t,xx(1:2), par);

% Re-name state at position 3 and 4 for clarity
dotx = xx(3); doty = xx(4);

% Calculate the rhs 
PBRFBP_rhs =  [dotx;
               doty;
               2*doty + dOM4(1);
              -2*dotx + dOM4(2)];  
end

function rhs_var = rhs_variational(t, yy, par)
% -------------------------------------------------------------------------
% This function allows to calculate the rhs of the PBRFBP dynamics coupled 
% with the variational propagation of the STM (coupled dynamics). All the 
% units are a-dimensional (the characteristic units are in the par struct).
% INPUT:
% - t            [1]: time 
% - yy           [20x1]: state and vectorized STM vector
% - par          [-]: struct of parmeters 
% OUTPUT:
% - rhs_var      [20x1]: evaluated RHS 
% DEPENDENCIES:
% PBRFBP_rhs(t,x,par)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Evaluate the rhs of the PBRFBP dynamics
f_rhs =  PBRFBP_rhs(t, yy(1:4), par);

% Define x and y variables from yy
x = yy(1); y = yy(2);

% Define the parameters from the par structure
om_s = par.om_s;  mu = par.mu; rho = par.rho; m_s = par.m_s;

% Calculate the Hessian of the potential 
d2OM = [(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - m_s/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) + (3*m_s*(2*x - 2*rho*cos(om_s*t))^2)/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2)^(5/2)) + 1, (3*m_s*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)); 
        (3*m_s*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - m_s/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) + (3*m_s*(2*y - 2*rho*sin(om_s*t))^2)/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1
       ];

% Build up the A(t) = df/dx matrix which defines the STM dynamics
AA = [zeros(2,2),     eye(2);
      d2OM      , [0,2;-2,0]];

% Calculate the rhs of the STM dynamics ( STM_dot = A*STM)
STM_dot = AA * reshape(yy(5:end), [4 4]);

% Vectorize the STM_dot matrix and concantenate with PBRFBP_rhs
rhs_var= [f_rhs; reshape(STM_dot, [16 1])]; 

end

function [OM4, dOM4] = OM4(t, xx, par)
% -------------------------------------------------------------------------
% This function allows to calculate the potential and its gradient for the
% PBRFBP written in the EMB-centered and rotating frame. All in
% adimensional units
% INPUT:
% - t            [1]: time 
% - yy           [4x1]: state vector
% - par          [-]: struct of parmeters 
%
% OUTPUT:
% - OM4          [1x1]: potential value
% - dOM4         [2x1]: gradient of the potential
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Extract problem parameters from par struct
om_s = par.om_s;  mu = par.mu; rho = par.rho; m_s = par.m_s;

% Extract position components from state vector
x = xx(1); y = xx(2); 

% Calculate distance of S/C from Earth
r1 = sqrt((x+mu)^2 + y^2); 

% Calculate distance of S/C from Earth
r2 = sqrt((x+mu-1)^2 + y^2);

% Calculate distance of S/C from Sun
r3 = sqrt((x-rho*cos(om_s*t))^2 + (y-rho*sin(om_s*t))^2);
    
% Calculate the potential
OM4 = 0.5*(x^2 + y^2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu) + m_s/r3 - (m_s/rho^2) * (x*cos(om_s*t) + y*sin(om_s*t));
  
% Calculate the gradient of the potential
dOM4 = [ x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - ...
        (m_s*cos(om_s*t))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - ...
        (m_s*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
         y - (m_s*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - ...
        (m_s*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) ...
         + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2)
       ];
end

function [xf, STM_f] = propagate_PBRFBP(x0, ti, ta, par)
% -------------------------------------------------------------------------
% This function allows to calculate both state and STM at final time, given
% the initial condition x_i, the initial time and final time. All in
% adimensional units.
% INPUT:
% - ti            [1]: initial time of propagation 
% - ta            [1]: final time of propagation
% - x0            [4x1]: initial condition (state cartesian vector)
% - par           [-]: struct of parmeters 
% OUTPUT:
% - x_f           [4x1]: state at final time (ta)
% - STM_f         [4x4]: STM at final time (ta)
% DEPENDENCIES:
% rhs_variational(t,y,par)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------    
% Construct the initial state of propagation (vectorize the STM at t_i)
y0 = [x0; reshape(eye(4), [16 1])];

% Set propagation tolerances
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Propagate the coupled state/STM dynamics
[~, yy] = ode78(@(t,y) rhs_variational(t,y,par), [ti ta], y0, options);

% Extract final state propagated
xf = yy(end,1:4)';

% Extract final STM propagated and make in matrix form
STM_f = reshape(yy(end, 5:20), [4 4]);

end

function twoDplot_rotframe(xx, mu, ri, rf, strName, indx_MS)
% -------------------------------------------------------------------------
% This function plots the trajectory xx in the EMB rotating frame.
% INPUT:
% - xx            [nx4]: trajectory
% - mu            [1]: EM gravitational parameter adimensionalized
% - ri            [1]: initial radius
% - rf            [1]: final radius
% - strName       [string]: title for the plotted trajectory 
% - indx_MS(OPTIONAL) [Nx1]: indexes for dots of the MS solution
%
% OUTPUT:
% - figure
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------    

% Set a default fontsize of 12
fontOptions(12);

% Create figure and give title
figure('Name', 'Orbit in the Rotating Frame');
sgtitle(strName);

% Add plot of trajectory
plot(xx(:,1), xx(:,2), 'LineWidth',2);
grid on
axis equal
hold on;
xlabel('$x$ [DU]'); ylabel('$y$ [DU]');

% Attractors 
plot(-mu, 0, 'o', 'MarkerSize',3, 'LineWidth',3);
plot(1-mu, 0, 'o', 'MarkerSize',3, 'LineWidth',3);

% Initial and final orbits
fimplicit(@(x,y) (x+mu).^2 + y.^2 - ri^2, 'LineStyle', '--', 'Color', "#808080" ,'LineWidth', 1.2);
plot(1-mu, 0); hold on;
fimplicit(@(x,y) (x+mu-1).^2 + y.^2 - rf^2, 'LineStyle', '--', 'Color', "#808080" ,'LineWidth', 1.2);

legend('Orbit','Earth','Moon');

% If in input --> plot the MS patch points
if nargin > 5
    for k = 1:length(indx_MS)
        plot(xx(indx_MS(k),1), xx(indx_MS(k),2), 'o','Color','r', 'HandleVisibility','off');
    end
end

% Put limits on x and y based on trajectory size
xlim([min(xx(:,1))-0.15, max(xx(:,1)) + 0.15]);
ylim([min(xx(:,2))-0.15, max(xx(:,2)) + 0.15]);

% Create a box for zoom around the Earth (initial orbit zoomed)
axes('position',[.4 .28 .15 .15])
box on % put box around new pair of axes
plot(xx(1:200,1), xx(1:200,2), 'LineWidth',2); % plot on new axes
hold on;
plot(-mu, 0, 'o', 'MarkerSize',3, 'LineWidth',3);
fimplicit(@(x,y) (x+mu).^2 + y.^2 - ri^2, 'LineStyle', '--', 'Color', "#808080" ,'LineWidth', 1.2);
axis tight
grid minor
axis equal
xlim([-0.04 0.01])
ylim([-0.02 0.03]);

% Create a box for zoom around the Moon (final orbit zoomed)
axes('position',[.65 .28 .15 .15])
box on % put box around new pair of axes
plot(xx(end-200:end,1), xx(end-200:end,2), 'LineWidth',2); % plot on new axes
hold on;
plot(1-mu, 0, 'o', 'MarkerSize',3, 'LineWidth',3, 'Color',[0.9290 0.6940 0.1250]);
fimplicit(@(x,y) (x+mu-1).^2 + y.^2 - rf^2, 'LineStyle', '--', 'Color', "#808080" ,'LineWidth', 1.2);
axis tight
grid minor
axis equal
xlim([(1-mu)-rf-5*1e-4 (1-mu)+rf+5*1e-4])
ylim([-5.5 5.5].*1e-3);

end

function twoDplot_ECinertial(xx, tt, mu, strName, indx_MS)
% -------------------------------------------------------------------------
% This function plots the trajectory xx in the Earth centered inertial frame.
% INPUT:
% - xx            [nx4]: trajectory in EMB rotating frame
% - tt            [nx1]: time span of trajectory
% - mu            [1]  : gravitational par. of CRTBP for EM
% - strName       [string]: title for the plotted trajectory 
% - indx_MS(OPTIONAL) [Nx1]: indexes for dots of the MS solution
%
% OUTPUT:
% - figure
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------   

% Set default font size
fontOptions(12);

% Convert the trajectory from EMB rotating frame to ECI frame
XX(:, 1) = (xx(:,1) + mu).*cos(tt) - xx(:,2).*sin(tt);
XX(:, 2) = (xx(:,1) + mu).*sin(tt) + xx(:,2).*cos(tt);

% Convert the Moon trajectory in EMB rotating frame to ECI frame
XX_moon(:,1) = (1 + mu)*cos(tt);
XX_moon(:,2) = (1 + mu)*sin(tt);

% Create figure and give the title
figure('Name','Undefined Orbit in ECI frame');
sgtitle(strName);

% Plot the trajectories
plot(XX(:,1), XX(:,2), 'LineWidth',2);
grid on
axis equal
hold on;
xlabel('$x$ [DU]'); ylabel('$y$ [DU]');
plot(0, 0, 'o', 'MarkerSize',3, 'LineWidth',3);
plot(XX_moon(:,1), XX_moon(:,2), '--', 'Color', [.7 .7 .7]);
plot(XX_moon(1,1), XX_moon(1,2), 'diamond', 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560])
plot(XX_moon(end,1), XX_moon(end,2), 'x', 'LineWidth', 2,'Color', [0.4940 0.1840 0.5560])
legend('Orbit','Earth', 'Moon Orbit','Moon @$t_i$', 'Moon @$t_f$');

% Limits on x axis and y axis based on size of S/C and Moon trajectory
xlim([min([XX(:,1); XX_moon(:,1)])-0.15, max([XX(:,1); XX_moon(:,1)]) + 0.15]);
ylim([min([XX(:,2); XX_moon(:,2)])-0.15, max([XX(:,2); XX_moon(:,2)]) + 0.15]);

% If requested plot MS patch points in red
if nargin > 4
    for k = 1:length(indx_MS)
        plot(XX(indx_MS(k),1), XX(indx_MS(k),2), 'o','Color','r', 'HandleVisibility','off');
    end
end

end

function [solution] = simple_shoot(y0, par, ri, rf, gradOptions)
% -------------------------------------------------------------------------
% This parent function allows to find an initial state for a
% trajectory between Earth and the Moon which minimizes the cost function 
% and satisfies the constraint starting from an initial guess. The NLP 
% problem is solved used using the Simple Shooting method 
% (both nested functions defined inside this parent function). 
% INPUT:
% - y0            [6x1]: NLP-variable initial guess 
% - par           [-]: struct of parmeters 
% - ri            [1]: constraint on initial circular orbit radius
% - rf            [1]: constraint on final circular orbit radius
% - gradOptions   [bool]: boolean value for the gradient (used in fmincon)
% OUTPUT:
% - solution      [-]: struct which contains information on found solution
% DEPENDENCIES:
% propagatePBRFBP(x,ti,tf,par), PBRFBP_rhs(t,x,par)
%
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------
% --- Define common variables used by both the nested functions (cost and
% constrain functions)

% Variable which contains the last NLP variable which has been used
yLast = [];
% Variable which contains the last propagated state (at final time)
xf = [];
% Variable which contains the last propagated STM (at final time)
STM = [];

    function [f, df] = cost_function(y)
        % -----------------------------------------------------------------
        % Nested function for the cost calculation
        % INPUT:
        % - y     [6x1]: NLP variable
        % OUTPUT:
        % - f     [1]: cost value
        % - df    [6x1]: gradient of cost
        % -----------------------------------------------------------------
       
        % Extract mu
        mu = par.mu;

        % Extract inital state from NLP variable
        xi = y(1:4);
        % Define repeated quantites based on xi
        dx1 = xi(3) - xi(2);
        dx2 = xi(4) + xi(1) + mu;

        % Extract inital time from NLP variable
        ti = y(5);
        % Extract final time from NLP variable
        ta = y(6);

        % Check if the same values of (xi, ti, tf) have been already used
        % for the same PBRFBP propagation
        if ~isequal(y, yLast)
            
            % Propagation of PBRFBP
            [xf, STM] = propagate_PBRFBP(xi, ti, ta, par);
            
            % Save the actual y value to the last y value used for the
            % propagation
            yLast = y;
        end

        % Define repeated quantities based on xf
        dx3 = xf(3) - xf(2);
        dx4 = xf(4) + xf(1) + mu - 1; 

        % Calculate the cost as sum of dv1 and dv2
        dv1 = sqrt( dx1^2 + dx2^2 ) - sqrt((1-mu)/ri);
        dv2 = sqrt( dx3^2 + dx4^2 ) - sqrt(mu/rf) ;

        %dv1 = sqrt( (xi(3) - xi(2))^2 + (xi(4) + xi(1) + par.mu)^2 ) - sqrt((1-par.mu)/ri);
        %dv2 = sqrt( (xf(3) - xf(2))^2 + (xf(4) + xf(1) + par.mu - 1)^2 ) - sqrt(par.mu/rf) ;

        f = dv1 + dv2;

        % If the requested output is > 1, then calculate the gradient of
        % cost
        if nargout > 1

            % Calculate the rhs of PBRFBP at ti
            rhsi = PBRFBP_rhs(ti, xi, par);
            % Calculate the rhs of PBRFBP at tf
            rhsf = PBRFBP_rhs(ta, xf, par);

            % Calculate the denominator of d(dvi)/d(xi)
            dvi = sqrt(dx1^2 + dx2^2);

            % Calculate the denominator of d(dvf)/d(xf)
            dvf = sqrt(dx3^2 +dx4^2);

            % Calculate the total d(dvi)/d(xi)
            ddvidxi = [dx2; -dx1; dx1; dx2]/dvi;

            % Calculate the total d(dvf)/d(xf)
            ddvfdxf = [dx4; -dx3; dx3; dx4]/dvf;

            % Calculate the full grad(f) (incorporating the STM)
            df = [ddvidxi + STM' * ddvfdxf; ddvfdxf' * (-STM*rhsi);  ddvfdxf' * rhsf];
        end

    end

    function [c, ceq, gradc, gradceq] = psi(y)
    % -----------------------------------------------------------------
    % Nested function for the non-linear constraints calculation (psi)
    % INPUT:
    % - y        [6x1]: NLP variable
    % OUTPUT:
    % - c          [-]: inequality constraints
    % - ceq      [4x1]: equality constraints
    % - gradc      [-]: Jacobian of inequality constraints
    % - gradceq  [6x4]: Jacobian of equality constraints
    % -----------------------------------------------------------------

        % Extract mu
        mu = par.mu;
        
        % Extract inital state from NLP variable
        xi = y(1:4);
        % Define repeated quantites based on xi
        dx1 = xi(3) - xi(2);
        dx2 = xi(4) + xi(1) + mu;
        % Extract inital time from NLP variable
        ti = y(5);
        % Extract final time from NLP variable
        ta = y(6);


        % Check if the same values of (xi, ti, tf) have been already used
        % for the same PBRFBP propagation
        if ~isequal(y, yLast)

            % Propagation of PBRFBP
            [xf, STM] = propagate_PBRFBP(xi, ti, ta, par);
            
            % Save the actual y value to the last y value used for the
            % propagation
            yLast = y;
        end

        % Define repeated quantities based on xf
        dx3 = xf(3) - xf(2);
        dx4 = xf(4) + xf(1) + mu - 1;  
        
       % Calculate the non-linear equality constraints
       ceq = [
            (xi(1) + mu)^2 + xi(2)^2 - ri^2;
            (xi(1) + mu)*dx1 + xi(2)*dx2;
            (xf(1) + mu - 1)^2 + xf(2)^2 - rf^2;
            (xf(1) + mu - 1)*dx3 + xf(2)*dx4;
            ];

        % Non-linear inequality constraints are not present
        c   = [];

        % If reqeusted, calculate Jacobian of non linear equality
        % contraints
        if nargout > 2
    
            % Calculate rhs at initial and final time
            rhs  = PBRFBP_rhs(ti, xi, par);
            rhsf = PBRFBP_rhs(ta, xf, par);

            % Allocate space for jacobian
            gradc = [];
            gradceq = zeros(6,4);

            % Calculate first two columns of gradient of non-lin equality
            gradceq(:,1) = [2*(xi(1) + mu); 2*xi(2); 0; 0; 0; 0];
            gradceq(:,2) = [xi(3); xi(4); xi(1) + mu; xi(2); 0; 0];

            % Construct the third column of non-lin equality
            dg3dxf = [2*(xf(1) + mu - 1); 2*xf(2); 0; 0];
            gradceq(:,3) = [STM' * dg3dxf; dg3dxf'*(-STM*rhs); dg3dxf'*rhsf];

            % Construct the fourth colum of non-lin equality
            dg4dxf = [xf(3); xf(4); xf(1) + mu - 1; xf(2)];
            gradceq(:,4) = [STM' * dg4dxf; dg4dxf'*(-STM*rhs); dg4dxf'*rhsf];
        end
    end

% Set options for the fmincon solver
optS = optimoptions(@fmincon,'SpecifyObjectiveGradient', gradOptions, 'SpecifyConstraintGradient',gradOptions, ...
    'Display', 'iter', 'Algorithm', 'active-set', 'ConstraintTolerance', 1e-10, 'StepTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, 'FunctionTolerance', 1e-6, ...
    'MaxFunctionEvaluations', 3e4);

% Set linear inequality constraints (ti < tf)s
A = [zeros(1,4) 1 -1];
b = 0;

% Set lb and ub
v_esc = sqrt(2*(1-par.mu)/ri);
ti_max = 2*pi/abs(par.om_s);
dT_max = 23;
lb = [-par.mu - ri; -ri; -v_esc; -v_esc; 0; 0];
ub = [-par.mu + ri;  ri;  v_esc;  v_esc; ti_max; ti_max + dT_max];

% Call the solver - fmincon
tic
[y_opt_as, dV_opt_as, exitFlag, output_as] = fmincon(@(y) cost_function(y), y0, A, b, [], [], lb, ub,  @(y) psi(y), optS);
toc

% --- Construct the solution struct ---

% Save found solution in terms of y (NLP variable)
solution.state    = y_opt_as;

% Save cost function evaluated at solution
solution.minvalue = dV_opt_as;

% Save the exitflag of the optimizer (for check of solution)
solution.exitflag = exitFlag;

% Save outpu of the optimizer (for check of solution)
solution.out      = output_as;

% Evaluate he constraint at the found solution and save them in solution
% struct
[~, PSI_eq, ~, ~] = psi(y_opt_as);
solution.c_i_at_sol = PSI_eq(1:2);
solution.c_f_at_sol = PSI_eq(3:4);

end

function [solution] = multiple_shoot(y0, par, ri, rf)
% -------------------------------------------------------------------------
% This function allows to find an initial state for a
% trajectory between Earth and the Moon which minimizes the cost function 
% and satisfies the constraint starting from an initial guess. The NLP 
% problem is solved used using the Multiple Shooting method 
%
% INPUT:
% - y0            [18x1]: NLP-variable initial guess 
% - par           [-]: struct of parmeters 
% - ri            [1]: constraint on initial circular orbit radius
% - rf            [1]: constraint on final circular orbit radius
% - gradOptions   [bool]: boolean value for the gradient (used in fmincon)
% OUTPUT:
% - solution      [-]: struct which contains information on found solution
% DEPENDENCIES:
% cost_function(y,mu,ri,rf), PBRFBP_rhs(y,par,ri,rf)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------   

% Set options for fmincon
optMS = optimoptions(@fmincon,'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient',true,'CheckGradients', false, ...
    'Display', 'iter', 'Algorithm', 'active-set', ...
    'ConstraintTolerance', 1e-10,  'StepTolerance', 1e-6, 'FunctionTolerance', 1e-6, ...
    'OptimalityTolerance', 1e-6, 'MaxIterations',3e4,'MaxFunctionEvaluations',3e4);

% solve the optimization problem with fmincon and specified options
tic
[y_opt_as, dV_opt_as, exitFlag, output_as] = fmincon(@(y)cost_function(y,par.mu,ri,rf), y0, [], [], [], [], [], [],  @(y)psi(y,par,ri,rf), optMS);
toc

% --- Build the solution struct

% Save the found optimal state
solution.state    = y_opt_as;
% Save the found objective function value
solution.minvalue = dV_opt_as;
% Save the exitflag
solution.exitflag = exitFlag;
% Save the output of the solver
solution.out      = output_as;
% Evaluate equality constraints at solution
[~, ceq] = psi(y_opt_as,par,ri,rf);
solution.psi_eq = ceq;

end

function [f, df] = cost_function(y,mu,ri,rf)
% -----------------------------------------------------------------
% MS cost function to be minimized.
% INPUT:
% - y     [18x1]: NLP variable
% OUTPUT:
% - f     [1]: cost value
% - df    [18x1]: gradient of cost
% DEPENDENCIES:
% NONE
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Extract the state relative to t1
x1 = y(1:4);

% Extract the state relative to t4
x4 = y(13:16);   

% Precompute differences (for efficiency)
dx1 = x1(3) - x1(2);
dy1 = x1(4) + x1(1) + mu;
dx4 = x4(3) - x4(2);
dy4 = x4(4) + x4(1) + mu - 1;

% Calculate cost function
dv1 = sqrt(dx1^2 + dy1^2) - sqrt((1-mu)/ri);
dv2 = sqrt(dx4^2 + dy4^2) - sqrt(mu/rf) ;

f = dv1 + dv2;

    % If requested, calculate the gradient of the cost
    if nargout > 1
    
        % Calculate denominator of d(dvi)/dxi and d(dvf)/dxf
        dvi = sqrt(dx1^2 + dy1^2);
        dvf = sqrt(dx4^2 + dy4^2);
    
        % Derivative of dvi wrt initial state xi (ddvidxi)
        ddvidxi = [dy1; -dx1; dx1; dy1]/dvi;
        
        % Derivative of dvf wrt initial state xf (ddvfdxf)
        ddvfdxf = [dy4; -dx4; dx4; dy4]/dvf;
    
        % Build-up the gradient of cost f
        df = [ddvidxi; zeros(8,1); ddvfdxf; zeros(2,1)];
    end

end
      
function [c, ceq, gradc, gradceq] = psi(y,par,ri,rf)
% -----------------------------------------------------------------
% MS constraint function (psi)
% INPUT:
% - y        [18x1]: NLP variable
% OUTPUT:
% - c         [9x1]: inequality constraints
% - ceq      [16x1]: equality constraints
% - gradc    [18x9]: Jacobian of inequality constraints
% - gradceq  [18x16]: Jacobian of equality constraints
% DEPENDENCIES:
% - PBRFBP_rhs(t,x,par), propagate_PBRFBP(t,x,par)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Extract mu
mu = par.mu;

%Extract DU, Re and Rm
DU = par.DU;
rE = par.rE;
rM = par.rM;

% Number of discrete points of multiple shooting
N = 4;
% Construct x as the ensemble of state at discretized times
x = reshape(y(1:16), 4, []);
% Extract initial time t1 and final time t4
t1 = y(17); t4 = y(18);
% Calculate time step
h = (t4-t1)/(N-1);

% Calculate uniform time span from t1 to t4 with dt
tt = t1:h:t4;
% Allocate space for the flow evaluated at t_i propagate from
% x_(i-1) with i = 2,3,4
phi = zeros(4,N-1);
% Allocate space for the rhs evaluated at phi_i, i=2,3,4
f_phi_eval = zeros(4,N-1);
% Allocate space for the STM evaluated at t_i, i=2,3,4
STM = zeros(4,4,N-1);

% RHS evaluated at x-variables: x1, x2, x3
f_var_eval = [PBRFBP_rhs(tt(1),x(:,1),par), PBRFBP_rhs(tt(2),x(:,2),par), PBRFBP_rhs(tt(3),x(:,3),par)];

% Loop over the N-1 propagation to perform
for i = 1:N-1
    
    % Propagate the state x_i from t_i to t_i+1
    [phi_i, STM_i] = propagate_PBRFBP(x(:,i), tt(i), tt(i+1), par);
    
    % Save the obtained flow
    phi(:,i) = phi_i;

    % Save the RHS evaluated with ODE flow 
    f_phi_eval(:,i) = PBRFBP_rhs(tt(i+1), phi(:,i), par);

    % Save the STM
    STM(:,:,i) = STM_i;
end


% Define the non-linear inequality constraints
c   = [
    (rE/DU)^2 - (x(1,1) + mu)^2 - x(2,1)^2;
    (rM/DU)^2 - (x(1,1) + mu - 1)^2 - x(2,1)^2;
    (rE/DU)^2 - (x(1,2) + mu)^2 - x(2,2)^2;
    (rM/DU)^2 - (x(1,2) + mu - 1)^2 - x(2,2)^2;
    (rE/DU)^2 - (x(1,3) + mu)^2 - x(2,3)^2;
    (rM/DU)^2 - (x(1,3) + mu - 1)^2 - x(2,3)^2;
    (rE/DU)^2 - (x(1,4) + mu)^2 - x(2,4)^2;
    (rM/DU)^2 - (x(1,4) + mu - 1)^2 - x(2,4)^2;
    t1 - t4;
    ];

% Define the non-linear equality constraints
ceq = [
    phi(:,1) - x(:,2);
    phi(:,2) - x(:,3);
    phi(:,3) - x(:,4);
    (x(1,1) + mu)^2 + x(2,1)^2 - ri^2;
    (x(1,1) + mu)*(x(3,1) - x(2,1)) + x(2,1)*(x(4,1) + x(1,1) + mu);
    (x(1,4) + mu - 1)^2 + x(2,4)^2 - rf^2;
    (x(1,4) + mu - 1)*(x(3,4) - x(2,4)) + x(2,4)*(x(4,4) + x(1,4) + mu - 1);
    ];

% If requested, calculate the jacobian of equality and ineqauliy
% constraints
if nargout > 2

    % Define 2 handle functions used for jacobian(ceq): derivatives of
    % j-th defect (j=2,3,4) with respect to t1 and t4
    dzjdt1 = @(j) ( (-STM(:,:,j)*f_var_eval(:,j))*(N-j)/(N-1) + f_phi_eval(:,j)*(N-1-j)/(N-1) )';
    dzjdtN = @(j) ( (-STM(:,:,j)*f_var_eval(:,j))*(j-1)/(N-1) + f_phi_eval(:,j)*j/(N-1) )';
    
    % Define handle function for the inequality constraints derivatives
    ineq = @(j) [-2*(x(1,j)+mu)   -2*x(2,j) 0 0;
                 -2*(x(1,j)+mu-1) -2*x(2,j) 0 0 ];

    % Define the gradient of the ineqaulity constraints by calls
    % to the handle function ineq
    gradc            = zeros(18,9);
    gradceq          = zeros(18,16);

    gradc(1:4,1:2)   = ineq(1)';
    gradc(5:8,3:4)   = ineq(2)';
    gradc(9:12,5:6)  = ineq(3)';
    gradc(13:16,7:8) = ineq(4)';
    % Gradient wrt to t1 and t4
    gradc(17:18, 9)  = [1; -1];

    % Define derivatives of the constraints psi_i and psi_f with
    % respect to x1 and x4
    dpsiidx1 = [2*(x(1,1)+mu) 2*x(2,1) 0 0; x(3,1) x(4,1) (x(1,1)+mu) x(2,1)]';
    dpsifdxN = [2*(x(1,4)+mu-1) 2*x(2,4) 0 0; x(3,4) x(4,4) (x(1,4)+mu-1) x(2,4)]';

    % Build the Jacobian of the non-linear equality constraints
    gradceq = [
               STM(:,:,1)'  zeros(4,4)    zeros(4,4)     dpsiidx1    zeros(4,2);
               -eye(4)      STM(:,:,2)'   zeros(4,4)     zeros(4,2)  zeros(4,2);
               zeros(4,4)     -eye(4)      STM(:,:,3)'   zeros(4,2)  zeros(4,2);
               zeros(4,4)     zeros(4,4)   -eye(4)        zeros(4,2)  dpsifdxN;
               dzjdt1(1)       dzjdt1(2)   dzjdt1(3)      zeros(1,2)  zeros(1,2);
               dzjdtN(1)       dzjdtN(2)   dzjdtN(3)      zeros(1,2)  zeros(1,2);
              ];

end

end

function fontOptions(DefaultaontSize)
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaulttextinterpreter', 'latex');
    set(0, 'defaultAxesFontSize', DefaultaontSize)
end
