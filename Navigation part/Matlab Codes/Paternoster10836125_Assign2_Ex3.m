clc; clearvars; cspice_kclear; close all

% -------------------------------------------------------------------------
% Spacecraft Guidance & Navigation (2024/2025)
% Assignment #2, Exercise 3
% Author: Daniele Paternoster
% -------------------------------------------------------------------------
% NOTE:
% - IMPORTANT: requires to have two kernel files (made by me for the GUESS
% coordinates of the lander): "lander_guess.bsp", "lander_guess.tf"
% - requires the assignment02.tm meta kernel in the same position
% - assignment02.tm requires also kernels\ folder in the same position
% -------------------------------------------------------------------------
% Paths of sgp4 utilities
addpath('sgp4');

% Upload meta kernel to the pool (contains all needed kernels in \kernels)
cspice_furnsh('assignment02.tm')

% Upload the created kernel for the lander position (lander_guess.bsp lander_guess.tf)
cspice_furnsh('lander_guess.bsp');
cspice_furnsh('lander_guess.tf');

%% Point 1: Visibility window

% Initial epoch definition t0 (UTC)
t0_utc = '2024-11-18T16:30:00.000';

% Initial epoch definition in ET et0
et0 = cspice_str2et(t0_utc);

% Final epoch definition tf (UTC)
tf_utc = '2024-11-18T20:30:00.000';

% Final epoch definition tf in ET etf
etf = cspice_str2et(tf_utc);

% Initial S/C state definition in Moon Centered Inertial RF (MCI) @et0
r0 = [4307.844185282820; -1317.980749248651; 2109.210101634011];
v0 = [-0.110997301537882; -0.509392750828585; 0.815198807994189];

% Complete SC state in MCI @et0
x0 = [r0; v0];

% --- Compute the SC trajectory in the MCI fram between et0 and etf ---

% Number of elements in time span: 30.0s span
n_el = round((etf-et0)/30.0) + 1;

% time span
et_span = linspace(et0, etf, n_el);

% Set propagation tolerances
optionsODE = odeset('RelTol', 1e-13, 'AbsTol', 1e-20);

% Propagation of  the orbiter state in MCI from et0 to etf
[~, xxORB_MCI] = ode113(@(t,x) kepMOON_rhs(0,x), et_span, x0, optionsODE);

% Lander initial guess position in MCMF-LATLON coordinates
LAT_land = 78; %[deg]
LONG_land = 15; %[deg]
ALT_land = 0; %[km]

% --- Find Lander initial guess position in MCI frame (from t0 to tf) ---

% Extract Moon radius (all three are equal, flatness computed just for completeness)
radius_M = cspice_bodvrd('MOON','RADII',3);
flat = (radius_M(1) - radius_M(3))/radius_M(1);

% Convert LATLONG coordinates into cartesian coordinates in MCMF frame
% Since flattness is 0, I could use cspice_reclat instead of pgrrec
rLAND_MCMF = cspice_pgrrec('MOON', cspice_rpd*LONG_land, cspice_rpd*LAT_land, 0, radius_M(1), flat);

% The lander is not moving in the MCMF frame
vLAND_MCMF = [0;0;0];

% State (r,v) transformation matrices for each time instant --> (MCMF2MCI) 
T_MCMF2MCI = cspice_sxform('IAU_MOON','J2000',et_span);

% Find the lander cartesian state in the MCI frame for each time instant
xxLAND_MCI = zeros(6,n_el);
for k = 1:n_el
    xxLAND_MCI(:,k) = T_MCMF2MCI(:,:,k)*[rLAND_MCMF; vLAND_MCMF];
end

% --- Check relative visibility between the lander and the orbiter ---

% Call the antenna_data function to get ideal observation of the orbiter
% from the lander
% NOTE: the Lander kernel was created with pinpoint and located in the
% kernels folder. The name of the created file are: 'lander_guess.tf' and
% 'lander_guess.bsp', with definition file saved as 'lander_guess.def'. The
% the SITE name is 'LANDER'.

[AZ_ORB, EL_ORB, RNG_ORB] = antenna_data('LANDER', et_span, xxORB_MCI);

% --- Plots of AZ, EL, RNG profiles --- 

% Create figure and set fontOptions (size and style)
figure('Name', 'AZ/EL/RANGE in visibility');
fontOptions(12);

% Compute hours past initial epoch in ET
delta_ET_hrs = (et_span - et0)/3600;

% Position the figure in center and resize it
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.35 0.6 0.45]);  

% Create a tiled layout for AZ EL and RANGE and set title of the figure
t1 = tiledlayout(1,3, 'TileSpacing','compact','Padding','compact');
title(t1, 'Orbiter Azimuth, Elevation and Range from Lander', ...
    'Reference epoch $t_0$ [UTC] : 2024-11-18T16:30:00.000', ...
    'interpreter', 'latex')

% Position the tiles containing AZ, EL and RANGE plots
nexttile;
plot(delta_ET_hrs, cspice_dpr*AZ_ORB, 'LineWidth',2)
grid minor
xlabel('$(t - t_0)$ [h]'); ylabel('AZ [deg]')
axis square

nexttile;
plot(delta_ET_hrs, cspice_dpr*EL_ORB,'LineWidth',2)
xlabel('$(t - t_0)$ [h]'); ylabel('EL [deg]');
yline(0, 'r--', 'LineWidth',2)
grid minor
legend('', '$EL_{min} \, [deg]$','Location','southeast')
axis square
ylim([-2, 100]);

nexttile;
plot(delta_ET_hrs,RNG_ORB,'LineWidth',2)
xlabel('$(t - t_0)$ [h]'); ylabel('RANGE [km]')
grid minor
axis square

%% Orbit plot


[~, xxfull_MCI] = ode113(@(t,x) kepMOON_rhs(0,x), [et0 et0 + 8.2*3600], x0, optionsODE);

fig_traj = figure;

% Plot Moon
[X, Y, Z] = sphere;
MoonRadii = cspice_bodvrd('MOON','RADII',3);
MoonSurface = surf(X*MoonRadii(1), Y*MoonRadii(2), Z*MoonRadii(3));
set(MoonSurface,'FaceColor',[0.5 0.5 0.5])
hold on; grid on; axis equal;

% Full trajectory of orbiter
plot3(xxfull_MCI(:, 1), xxfull_MCI(:, 2),xxfull_MCI(:, 3), 'r-.','LineWidth',1.3);

% In-visibility Trajectory of orbiter from lander
plot3(xxORB_MCI(:, 1), xxORB_MCI(:, 2), xxORB_MCI(:, 3), 'r','LineWidth',2);
plot3(xxORB_MCI(1, 1), xxORB_MCI(1, 2), xxORB_MCI(1, 3),'o', 'Color','r','LineWidth',2);
plot3(xxORB_MCI(end, 1), xxORB_MCI(end, 2), xxORB_MCI(end, 3),'+','Color','r','LineWidth',2);

% In-visibility Trajectory of orbiter from lander
%plot3(xxLAND_MCI(1,:), xxLAND_MCI(2,:), xxLAND_MCI(3,:), 'b','LineWidth',2);
plot3(xxLAND_MCI(1,1), xxLAND_MCI(2,1), xxLAND_MCI(3,1),'o', 'Color','b','LineWidth',2);
%plot3(xxLAND_MCI(1,end), xxLAND_MCI(2,end), xxLAND_MCI(3,end),'+','Color','b','LineWidth',2);

title('Orbiter Trajectory in MCI frame')
xlabel('x [km]'); ylabel('y [km]'); zlabel('z [km]');

legend('Moon','Full orbit', '$t_0$/$t_f$', 'S/C @$t_0$', 'S/C @$t_f$', 'Lander Pos', ...
    'NumColumns',1,'FontSize',10,'Location','best')


%% Point 2: Measurements simulation

% Noise for orbiter position in MCI and relative range (ORBITER - LANDER)
sigma_pos = 0.1; %[km]

% Add noise to position of the orbiter in MCI
pos_COV = diag(sigma_pos^2*ones(3,1));
rr_ORB_real = mvnrnd(xxORB_MCI(:,1:3), pos_COV);

% Retrieve lander position from the ground truth kernel - moon_lander
[rr_LAND_real, ~] = cspice_spkpos('MOONLANDER', et_span,'J2000','NONE','MOON');

%% Compute Range as norm difference between orbiter and lander pos. in J2000

RNG_ORB_mean = vecnorm(xxORB_MCI(:,1:3) - rr_LAND_real',2,2);

% Add noise to relative range between orbiter and lander
range_COV = sigma_pos^2;
RNG_ORB_real = mvnrnd(RNG_ORB_mean, range_COV);

% Graphical verification of 3SIGMA boundaries

% Create figure for the simulated measured orbiter position in MCI
figure('Name','Measured orbiter position in MCI');

% Set figure position and size
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.05 0.35 0.6 0.45]);  

% Create tiled layout object and title the figure
t2 = tiledlayout(1,3, 'TileSpacing','compact','Padding','compact');
title(t2, ['Residuals of simulated measurements of Orbiter position in ' ...
    'MCI frame'], 'Reference epoch $t_0$ [UTC] : 2024-11-18T16:30:00.000', ...
    'interpreter', 'latex')

% Create tiles and respective graphs
nexttile;
plot(delta_ET_hrs, rr_ORB_real(:,1) - xxORB_MCI(:,1), '*')
hold on
plot(delta_ET_hrs, +3*sigma_pos*ones(size(delta_ET_hrs)), 'r--')
plot(delta_ET_hrs, -3*sigma_pos*ones(size(delta_ET_hrs)), 'r--')
grid minor
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_{x,MCI}$ [km]');
axis square

nexttile;
plot(delta_ET_hrs, rr_ORB_real(:,2) - xxORB_MCI(:,2) , '*')
hold on
plot(delta_ET_hrs,  +3*sigma_pos*ones(size(delta_ET_hrs)), 'r--','LineWidth',1)
plot(delta_ET_hrs,  -3*sigma_pos*ones(size(delta_ET_hrs)), 'r--','LineWidth',1)
grid minor
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_{y,MCI}$ [km]');
axis square

nexttile;
plot(delta_ET_hrs, rr_ORB_real(:,3) - xxORB_MCI(:,3), '*')
hold on
plot(delta_ET_hrs, +3*sigma_pos*ones(size(delta_ET_hrs)), 'r--','LineWidth',1)
plot(delta_ET_hrs, -3*sigma_pos*ones(size(delta_ET_hrs)), 'r--','LineWidth',1)
grid minor
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_{z,MCI}$ [km]');
axis square

% Measured range plot
figure('Name','Relative range ORB - LANDER measurement')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.65 0.35 0.3 0.55]);  
titleRANGE = sgtitle('Residualts of simulated range measurement: ORB-LANDER','FontSize', 12);
plot(delta_ET_hrs, RNG_ORB_real - RNG_ORB_mean, '*');
hold on
plot(delta_ET_hrs, 3*sigma_pos*ones(size(delta_ET_hrs)), 'r--','LineWidth',1)
plot(delta_ET_hrs, -3*sigma_pos*ones(size(delta_ET_hrs)), 'r--','LineWidth',1)
grid minor
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_{\rho}$ [km]');
ylim([-0.4 0.4]);
xlim([0 4])
axis square

%% Point 3: estimate the lunar orbiter state in MCI with UKF

% Initialization variables for the UKF filter
P0 = diag([10,1,1,0.001,0.001,0.001]);
x0_pert = mvnrnd(x0, P0)';

% Measurements noise covariance matrix (relative to the noisy position
% vector of the orbiter in MCI frame)
R = diag(sigma_pos^2 * ones(3,1));

% Sequence of measurements in chronological order (orbiter position in MCI)
Ym = rr_ORB_real';

% Call UKF: estimate the orbiter full state and covariance matrix (in MCI)
[x, P] = UnscentedKalmanFilter(x0_pert, P0, et_span, Ym, R);

%% Post-processing on UKF state estimation

% Error estimate considering xxORB as the truth
err_r = vecnorm(xxORB_MCI(:,1:3) - x(1:3,:)',2,2);
err_v = vecnorm(xxORB_MCI(:,4:6) - x(4:6,:)',2,2);

% Retrieve the sigma values from the filter-estimated covariance
for k = 1:length(x)
    sigma_r(k) = squeeze(sqrt(max(eig(P(1:3,1:3,k)))));
    sigma_v(k) = squeeze(sqrt(max(eig(P(4:6,4:6,k)))));
end

% --- Error and 3-sigma plots --- 

fontOptions(14);

figure('Name','Estimated error - UKF');
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.25 0.6 0.65]);  

t3 = tiledlayout(2,5, 'TileSpacing','compact','Padding','compact');
title(t3, 'Estimated error for $r_{i,MCI}$ and $v_{i,MCI}$ with $3\sigma$ bounds', ...
     'Reference epoch $t_0$ [UTC] : 2024-11-18T16:30:00.000','Interpreter','latex')

nexttile([1 4]);
semilogy(delta_ET_hrs, 3*sigma_r,'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs, err_r,'LineWidth',1.3)
% semilogy(delta_ET_hrs, err_x,'LineWidth',1.3)
% semilogy(delta_ET_hrs, err_y,'LineWidth',1.3)
% semilogy(delta_ET_hrs, err_z,'LineWidth',1.3)
grid minor
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_r$ [km]');
%legend('3$\sigma_r$','$\epsilon_x$','$\epsilon_y$','$\epsilon_z$','Location','Northeast','NumColumns',4)
legend('3$\sigma_r$','$||\epsilon_r||$','Location','Northeast','NumColumns',2)
xlim([-0.05 4]);
set(gca, 'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

nexttile([1 1])
semilogy(delta_ET_hrs(1:10), 3*sigma_r(1:10),'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs(1:10), err_r(1:10),'LineWidth',1.3)
% semilogy(delta_ET_hrs(1:10), err_x(1:10),'LineWidth',1.3)
% semilogy(delta_ET_hrs(1:10), err_y(1:10),'LineWidth',1.3)
% semilogy(delta_ET_hrs(1:10), err_z(1:10),'LineWidth',1.3)
grid minor
set(gca,'YAxisLocation','right')
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_r$ [km]');
set(gca, 'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

nexttile([1 4]);
semilogy(delta_ET_hrs, 3*sigma_v,'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs, err_v,'LineWidth',1.3)
% semilogy(delta_ET_hrs, err_vx,'LineWidth',1.3)
% semilogy(delta_ET_hrs, err_vy,'LineWidth',1.3)
% semilogy(delta_ET_hrs, err_vz,'LineWidth',1.3)
grid minor
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_v$ [km/s]');
%legend('3$\sigma_v$','$\epsilon_{vx}$','$\epsilon_{vy}$','$\epsilon_{vz}$','Location','Northeast','NumColumns',4)
legend('3$\sigma_v$','$||\epsilon_v||$','Location','Northeast','NumColumns',2)
xlim([-0.05 4]);
set(gca, 'YTick',[1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1])

nexttile([1 1])
semilogy(delta_ET_hrs(1:10), 3*sigma_v(1:10),'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs(1:10), err_v(1:10),'LineWidth',1.3)
% semilogy(delta_ET_hrs(1:10), err_vx(1:10),'LineWidth',1.3)
% semilogy(delta_ET_hrs(1:10), err_vy(1:10),'LineWidth',1.3)
% semilogy(delta_ET_hrs(1:10), err_vz(1:10),'LineWidth',1.3)
grid minor
set(gca,'YAxisLocation','right')
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_v$ [km/s]');
set(gca, 'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

%% Point 4: estimate the lunar orbiter state in MCI and Lander LATLONG with UKF

% Initialization variables for the UKF filter WITH parameter estimation
P0_aug = diag([10,1,1,0.001,0.001,0.001,0.00001,0.00001]);
x0_aug = mvnrnd([x0; LAT_land*cspice_rpd; LONG_land*cspice_rpd], P0_aug)';

% Sequence of measurements in chronological order - orb-position in MCI and
% range measurement ORB-LANDER
Ym_aug = [Ym; RNG_ORB_real'];

% Measurements noise covariance matrix (relative to the noisy position
% vector of the orbiter in MCI frame and relative range ORB-LANDER)
R_aug  = diag(sigma_pos^2 * ones(4,1));

% Call the UKF function with Lander Position estimation
[x_aug, P_aug, yyy] = UnscentedKalmanFilterLandEst(x0_aug, P0_aug, et_span, Ym_aug, R_aug);

% --- Extract UKF estimation of LATLONG of lander at initial and final
% epoch ---

% Initial LATLONG (initialization variables)
LATLONG0 = x_aug(7:8, 1);
P0_LATLON = P_aug(7:8,7:8,1);

% Final mean and covariance after filtering process
LATLONGf = x_aug(7:8, end);
Pf_LATLON = P_aug(7:8,7:8,end);

%% Lat Long errors - plots

% Retrieve the true lat-lon of the lander from the moon-lander kernel
[r_land_truth, ~] = cspice_spkpos('MOONLANDER', et_span,'IAU_MOON','NONE','MOON');
[rad_truth, lon_truth, lat_truth ] = cspice_reclat(r_land_truth);

% Error estimate considering xxORB as the truth
err_r_wpe = vecnorm(xxORB_MCI(:,1:3) - x_aug(1:3,:)',2,2);
err_v_wpe = vecnorm(xxORB_MCI(:,4:6) - x_aug(4:6,:)',2,2);

% Retrieve the sigma values from the filter-estimated covariance (WPE: with
% parameters estimation)

for k = 1:length(x_aug)
    sigma_r_wpe(k) = squeeze(sqrt(max(eig(P_aug(1:3,1:3,k)))));
    sigma_v_wpe(k) = squeeze(sqrt(max(eig(P_aug(4:6,4:6,k)))));
    sigma_p(k)     = squeeze(sqrt(max(eig(P_aug(7:8,7:8,k)))));
end

% Plot 3sigma bounds and norm errors for r and v

fontOptions(14);

figure('Name','Estimated error - UKF with parameter estimation');

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2 0.10 0.6 0.85]);  

t4 = tiledlayout(3,5, 'TileSpacing','compact','Padding','compact');
title(t4, 'Estimated error for $r_{i,MCI}$, $v_{i,MCI}$, LAT and LONG with $3\sigma$ bounds', ...
     'Reference epoch $t_0$ [UTC] : 2024-11-18T16:30:00.000','Interpreter','latex')

nexttile([1 4]);
semilogy(delta_ET_hrs, 3*sigma_r_wpe,'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs, err_r_wpe,'LineWidth',1.3)
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_r$ [km]');
grid minor
legend('$3\sigma_r$','$||\epsilon_r||$','Location','Northeast')
xlim([0 4]);
set(gca, 'YTick',[1e-2 1e-1 1e0])

nexttile([1 1])
semilogy(delta_ET_hrs(1:10), 3*sigma_r_wpe(1:10),'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs(1:10), err_r_wpe(1:10),'LineWidth',1.3)
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_r$ [km]');
grid minor
set(gca,'YAxisLocation','right')

nexttile([1 4]);
semilogy(delta_ET_hrs, 3*sigma_v_wpe,'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs, err_v_wpe,'LineWidth',1.3)
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_v$ [km/s]');
grid minor
legend('$3\sigma_v$','$||\epsilon_v||$','Location','Northeast')
xlim([0 4]);
set(gca, 'YTick',[1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

nexttile([1 1])
semilogy(delta_ET_hrs(1:10), 3*sigma_v_wpe(1:10),'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs(1:10), err_v_wpe(1:10),'LineWidth',1.3)
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_r$ [km/s]');
grid minor
set(gca,'YAxisLocation','right')
set(gca, 'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

% Plot 3sigma bounds and errors lat and lon

nexttile([1 4]);
semilogy(delta_ET_hrs, 3*sigma_p*cspice_dpr,'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs, 3*sqrt(squeeze(P_aug(end-1, end-1,:)))*cspice_dpr)
semilogy(delta_ET_hrs, 3*sqrt(squeeze(P_aug(end, end,:)))*cspice_dpr)
semilogy(delta_ET_hrs,abs(x_aug(7,:)-lat_truth)*cspice_dpr,'LineWidth',1.3);
semilogy(delta_ET_hrs,abs(x_aug(8,:)-lon_truth)*cspice_dpr,'LineWidth',1.3);
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_{lat,long}$ [deg]');
grid minor
legend('$3\sigma_{p}$','','','LAT','LONG','Location','Northeast','NumColumns',3);

xlim([0 4]);
set(gca, 'YTick',[1e-5 1e-4 1e-3 1e-2 1e-1 1e0])

nexttile([1 1])
semilogy(delta_ET_hrs(1:10), 3*sigma_p(1:10)*cspice_dpr,'LineWidth',1.5,'LineStyle','--')
hold on;
semilogy(delta_ET_hrs(1:10),abs(x_aug(7,1:10)-lat_truth(1:10))*cspice_dpr,'LineWidth',1.3);
semilogy(delta_ET_hrs(1:10),abs(x_aug(8,1:10)-lon_truth(1:10))*cspice_dpr,'LineWidth',1.3);
xlabel('$(t - t_0)$ [h]'); ylabel('$\epsilon_{lat,long}$ [deg]');
grid minor
set(gca,'YAxisLocation','right')


%% Plots with LATLONG ground truth

fontOptions(12)

figure;

t5 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(delta_ET_hrs, lat_truth*cspice_dpr,'LineWidth',1.5,'LineStyle','--');
hold on
plot(delta_ET_hrs, x_aug(7,:)*cspice_dpr,'LineWidth',1.3);
xlabel('$(t - t_0)$ [h]'); ylabel('LAT [deg]');
legend('LAT truth', 'LAT est.');
xlim([-0.05 4]);
grid on;

nexttile
plot(delta_ET_hrs, lon_truth*cspice_dpr,'LineWidth',1.5,'LineStyle','--');
hold on
plot(delta_ET_hrs, x_aug(8,:)*cspice_dpr,'LineWidth',1.3);
xlabel('$(t - t_0)$ [h]'); ylabel('LONG [deg]');
legend('LONG truth', 'LONG est.');
xlim([-0.05 4]);
grid on;

%% Correlation on the parameters covariance submatrix 

for k = 1:length(P_aug)
    P_latlon_corr(:,:,k) = corrcov(P_aug(end-1:end,end-1:end,k));
end

figure;
plot(delta_ET_hrs, squeeze(P_latlon_corr(1,2,:)),'LineWidth',2);
xlabel('$(t - t_0)$ [h]'); ylabel('$\rho_{P}$ [-]');
legend('Correlation parameter - LATLON covariance submatrix');
xlim([-0.05 4.2]);
ylim([-0.05 1.2]);
grid on;
%% Functions

function f = kepMOON_rhs(t, x)
%--------------------------------------------------------------------------
% This function allows to calculate the RHS for the keplerian dynamics
% around the Moon written in MCI frame.
% -------------------------------------------------------------------------
% INPUT:
% - t              [1]: time
% - x            [6x1]: orbital state of the S/C in MCI frame
% OUTPUT:
% - f            [6x1]: rhs 
% DEPENDECIES:
% - it requires to have the gm_de432.tpc kernel uploaded
%--------------------------------------------------------------------------

% Allocate rhs 
f = zeros(6,1);

% Assign position and velocity from the full state
r = x(1:3);
v = x(4:6);

% Calculate norm of the position vector
rnorm = norm(r);

% Get the gravitational parameter for Moon using spice routines
muM = cspice_bodvrd('MOON', 'GM', 1);

% Calculate the two sub-vectors of the rhs 
f(1:3) = v;
f(4:6) = -muM*r/rnorm^3;

end

function [AZ, EL, RANGE] = antenna_data(station, et_span, x_MCI)
% -------------------------------------------------------------------------
% This function allows to obtain the ideal values of AZ, EL, RANGE 
% of a space object with respect to a fixed point on the moon, defined by 
% a spice kernel.
% NOTE: visibility constraint are not applied here.
%
% INPUT:
% - station name [string]
% - et_span      [n]: ET span of interest
% - x_MCI        [nx6]: orbital state of the space object in ECI frame
%
% OUTPUT:
% - AZ           [1xn] AZ vector
% - EL           [1xn] EL vector
% - RANGE        [1xn] RANGE vector
% 
% -------------------------------------------------------------------------


% Initialize output quantities
AZ = zeros(size(et_span));
EL = AZ;
RANGE = AZ;

% Define topocentric frame
topocentric_frame = [station, '_TOPO'];

% Loop over epochs in et_span
for k = 1:length(et_span)

    % Define the k-th epoch of interest
    et_k = et_span(k);

    % Get ECI position of the station (from input)
    r_station_MCI = cspice_spkpos(station, et_k, 'J2000', 'NONE', 'MOON');

    % Get ECI position of SC wrt to station
    r_sat_station_MCI = x_MCI(k,1:3)' - r_station_MCI;
    
    % Transformation from ECI to topocentric frame
    ROT_MCI2TOPO = cspice_pxform('J2000', topocentric_frame, et_k);
    
    % Convert state into topocentric frame
    r_sat_station_TOPO = ROT_MCI2TOPO*r_sat_station_MCI;
    
    % Obtain Azimuth and Elevation coordinates with spice reclat routine
    [RANGE(k), AZ(k),  EL(k)] = cspice_reclat(r_sat_station_TOPO);

end

end

function [RANGE] = antenna_data_UKF(LAT, LON, et, x_MCI)
% -------------------------------------------------------------------------
% This function allows to obtain the ideal value of RANGE seen from a
% station (LATLON defined) on the Moon.
% NOTE: visibility constraint are not applied here.
%
% INPUT:
% - LATLON       [2x1]: Latitude and Longitude of the observation site
% - et           [1]: epoch of the measurement
% - x_MCI        [6x1]
% OUTPUT:
% - RANGE        [1] RANGE value
% 
% -------------------------------------------------------------------------

% Get Moon information for cspice_pgrrec 
radiiMoon = cspice_bodvrd('MOON', 'RADII', 3);
re = radiiMoon(1); rp = radiiMoon(3); 
flat = (re - rp)/re;

% LAT LONG --> in radians 

% Get lander site position
pos_lander_MCMF = cspice_pgrrec('MOON', LON, LAT, 0, re, flat);

% Convert Orbiter Full state from MCI2MCMF
T_MCI2MCMF = cspice_pxform('J2000', 'IAU_MOON', et);
r_MCMF = T_MCI2MCMF*x_MCI(1:3);

% Get relative position of orbiter wrt lander in MCMF
pos_ORBwrtLAND_MCMF = r_MCMF - pos_lander_MCMF;

RANGE = norm(pos_ORBwrtLAND_MCMF);

end

function fontOptions(DefaultFontSize)
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaulttextinterpreter', 'latex');
    set(0, 'defaultAxesFontSize', DefaultFontSize)
end

function [x, P] = UnscentedKalmanFilter(x0, P0, et_span, Ym, R)
% -------------------------------------------------------------------------
% This function is the implementation of a UKF applied in the context of
% sequential orbital state estimation for a lunar orbiter given the 
% inertial position measurements.
% INPUT:
% - x0           [6x1]: initial state at t0 epoch (
% - P0           [6x6]: covariance of the state at t0
% - et_span      [nx1]: time span over which measurements are taken
% - Ym           [3xn]: matrix of measurements in chronological order
% - R            [3x3]: measurements noise covariance matrix 
% OUTPUT:
% - x            [6xn]: estimated orbiter state vector in MCI 
% - P            [6x6xn]: estimated covariance of x in MCI
% -------------------------------------------------------------------------

% Initialize estimated state vector size
n = length(x0);
n_el = length(et_span)-1;

% --- Allocate memory for the variables ---
% Nomenclature (holds for state x, measurements y and state x covariance):
% km1: updated(+) k-1 quantity (a-posteriori from previous step)
% kp: predicted(-) k quantity (a-priori actual step)
% k: updated(+) k quantity (a-posteriori actual step)
% Covariances nomenclature:
% P_k_ee: innovation covariance of step k
% P_k_xy: cross-correlation covariance at step k

x_km1   = x0;
x_kp    = zeros(n,1);
y_kp    = zeros(3,1);
x_k     = [];

% First mean element (at t0)
x       = x0;

P_km1    = P0;
P_kp     = zeros(n,n);
P_k_ee   = R;
P_k_xy   = zeros(n,3);
P_k      = [];
P        = zeros(6,6,n_el + 1);

% First covariance element (at t0)
P(:,:,1) = P0;

% Allocate memory for weight vectors of the UT.
% Mean values weights
W_mean        = zeros(2*n+1,1);
% Covariances weights
W_cov         = W_mean;

%.--- Calculate sigma points --- 
% UT parameters
alfa = 0.01;
beta = 2;
lambda = (alfa^2 - 1)*n;

% Calculate UT weights
W_mean(1) = lambda/(n+lambda);
W_cov(1)  = lambda/(n+lambda) + (1 - alfa^2 + beta);
W_mean(2:end) = 1/(2*(n+lambda));
W_cov(2:end) = W_mean(2:end);

% Allocation of variables (sigma points, evaluation of sigma points, weights)
sigma_pts_km1 = zeros(n,2*n+1);
sigma_pts_k   = zeros(n,2*n+1);
gamma_pts_k   = zeros(3,2*n+1); 

% Options for propagation of ODE
options = odeset('RelTol',1e-13, 'AbsTol', 1e-20);

% Loop over the measurements (index k of loop is one unit behind the k used
% for nomenclature)
for k = 1:n_el

% --- Sigma points ---
% Calculate square root of matrix
sqrtMAT   = sqrtm((n+lambda)*P_km1);

% Definition of first sigma point - sigma0 (1)
sigma_pts_km1(:,1) = x_km1;
[~,sigma_temp] = ode113(@(t,x) kepMOON_rhs(0,x), [et_span(k) et_span(k+1)], sigma_pts_km1(:,1),options);
sigma_pts_k(:,1) = sigma_temp(end, :)';

% Definition of first gamma points
gamma_pts_k(:,1) = sigma_pts_k(1:3,1); %[AZ_temp; EL_temp; RNG_temp];

    % Definition of sigma points - sigma_i (2N) 
    for s = 2:n+1
        % Define the sigma points at k-1
        sigma_pts_km1(:,s)   = x_km1 + sqrtMAT(:,s-1);
        sigma_pts_km1(:,n+s) = x_km1 - sqrtMAT(:,s-1);
        
        % Propagate the other 1-->n sigma points at k
        [~,sigma_temp] = ode113(@(t,x) kepMOON_rhs(0,x), [et_span(k) et_span(k+1)], sigma_pts_km1(:,s),options);
        sigma_pts_k(:,s) = sigma_temp(end, :)';

        % Apply measurement equation: basis gamma point
        gamma_pts_k(:,s) = sigma_pts_k(1:3,s); 
    
        % Propagate the other 1-->n sigma points at k
        [~,sigma_temp] = ode113(@(t,x) kepMOON_rhs(0,x), [et_span(k) et_span(k+1)], sigma_pts_km1(:,n+s),options);
        sigma_pts_k(:,n+s) = sigma_temp(end, :)';
        
        % Apply measurement equation: basis gamma point
        gamma_pts_k(:,n+s) = sigma_pts_k(1:3,n+s); 
    end
   
    % --- Apply UT weigths to sigma and gamma points ---
    
    % Calculate the predicted quantites with weighted mean
    for l = 1:(2*n+1)
        % Calculate the predicted mean state (sample mean)
        x_kp = x_kp + W_mean(l) * sigma_pts_k(:,l);
        % Calculate the predicted measurements (sample mean)
        y_kp = y_kp + W_mean(l) * gamma_pts_k(:,l);
    end
    
    % Calculate predicted state covariance, innovation covariance and
    % cross-covariance as weighted sums
    for l = 1:(2*n+1)
        % Calculate the predicted covariance (sample covariance)
        P_kp   = P_kp    + W_cov(l) * (sigma_pts_k(:,l) - x_kp)*(sigma_pts_k(:,l) - x_kp)';
        % Calculate the predicted covariance for measurements
        P_k_ee = P_k_ee  + W_cov(l) * (gamma_pts_k(:,l) - y_kp)*(gamma_pts_k(:,l) - y_kp)';
        % Calculate the cross covariance
        P_k_xy = P_k_xy  + W_cov(l) * (sigma_pts_k(:,l) - x_kp)*(gamma_pts_k(:,l) - y_kp)';
    end
    
    % Compute Kalman Gain
    K_k = P_k_xy * (P_k_ee\eye(size(R)));
    
    % Update step
    x_k = x_kp + K_k*(Ym(:,k+1) - y_kp);
    P_k = P_kp - K_k*P_k_ee*K_k';
    
    % Store sequential updates of the UKF
    x          = [x, x_k];
    P(:,:,k+1) = P_k;
    
    % Update new k-1 variables to k
    P_km1 = P_k;
    x_km1 = x_k;
    
    % Re-allocate the variables to be calculated in the next step
    P_kp    = zeros(n,n);
    P_k_ee  = R;
    P_k_xy  = zeros(n,3);
    x_kp    = zeros(n,1);
    y_kp    = zeros(3,1);

end
end

function [x_aug, P_aug, y_pred] = UnscentedKalmanFilterLandEst(x0_aug, P0_aug, et_span, Ym, R)
% -------------------------------------------------------------------------
% This function is the implementation of a UKF applied in the context of
% sequential orbital state estimation for a lunar orbiter and parameter 
% estimation (LATLON of the lander position) given the  inertial
% position measurements of the SC and the relative range bw orbiter and
% lander.
% INPUT:
% - x0_aug       [8x1]: augmented state at t0 (x_orbiter; lat; lon)
% - P0           [8x8]: Covariance matrix of augmented state at t0
% - et_span      [nx1]: time span over which measurements are taken
% - Ym           [4xn]: matrix of measurements in chronological order
% - R            [4x4]: measurements noise covariance matrix 
% OUTPUT:
% - x_aug        [8xn]: estimated augmented state
% - P_aug        [8x8xn]: estimated covariance of the augmented state
% -------------------------------------------------------------------------

% Define the size of the augmented state and the number of measurements
n = length(x0_aug);
n_el = length(et_span)-1;

% --- Allocate memory for the variables ---
% Nomenclature (holds for state x, measurements y and state x covariance):
% km1: updated(+) k-1 quantity (a-posteriori from previous step)
% kp: predicted(-) k quantity (a-priori actual step)
% k: updated(+) k quantity (a-posteriori actual step)
% Covariances nomenclature:
% P_k_ee: innovation covariance of step k
% P_k_xy: cross-correlation covariance at step k

x_km1   = x0_aug;
x_kp    = zeros(n,1);
y_kp    = zeros(4,1);
x_k     = [];
x_aug   = x0_aug;

P_km1    = P0_aug;
P_kp     = zeros(n,n);
P_k_ee   = R;
P_k_xy   = zeros(n,4);
P_k      = [];
P_aug       = zeros(n,n,n_el + 1);
P_aug(:,:,1) = P0_aug;

% Allocate memory for weight vectors of the UT.
% Mean values weights
W_mean        = zeros(2*n+1,1);
% Covariances weights
W_cov         = W_mean;

%.--- Calculate sigma points --- 
% UT parameters
alfa = 0.01;
beta = 2;
lambda = (alfa^2 - 1)*n;

% Calculate UT weights
W_mean(1) = lambda/(n+lambda);
W_cov(1)  = lambda/(n+lambda) + (1 - alfa^2 + beta);
W_mean(2:end) = 1/(2*(n+lambda));
W_cov(2:end) = W_mean(2:end);

% Allocation of variables (sigma points, evaluation of sigma points)
sigma_pts_km1 = zeros(n,2*n+1);
sigma_pts_k   = zeros(n,2*n+1);
gamma_pts_k   = zeros(4,2*n+1); 

% Options for propagation of ODE
options = odeset('RelTol',1e-13, 'AbsTol', 1e-20);

% Loop over the measurements (index k of loop is one unit behind the k used
% for nomenclature)
for k = 1:n_el

% --- Sigma points ---
% Calculate square root of aposteriori covariance matrix at previous step
sqrtMAT   = sqrtm((n+lambda)*P_km1);

% --- Definition of first sigma point - sigma0 (1)
% Sigma point at k-1
sigma_pts_km1(:,1) = x_km1;

% --- Obtain sigma at k
[~,sigma_temp] = ode113(@(t,x) kepMOON_rhs(0,x), [et_span(k) et_span(k+1)], sigma_pts_km1(1:6,1),options);
% first six elements are for the orbiter state
sigma_pts_k(1:6,1) = sigma_temp(end,:)';
% last two elements are for latitude and longitude
sigma_pts_k(7:8,1) = sigma_pts_km1(7:8,1);

% Find gamma points (range and orbiter position) given sigma points at t_k
RNG_temp = antenna_data_UKF(sigma_pts_k(7,1), sigma_pts_k(8,1), et_span(k+1), sigma_pts_k(1:6,1));
% Definition of first gamma points
gamma_pts_k(1:3,1) = sigma_pts_k(1:3,1); 
gamma_pts_k(4,1)   = RNG_temp; 

% Definition of the other sigma points - sigma_i (2N) 
for s = 2:n+1
    % Define the sigma points at k-1
    sigma_pts_km1(:,s)   = x_km1 + sqrtMAT(:,s-1);
    sigma_pts_km1(:,n+s) = x_km1 - sqrtMAT(:,s-1);
    
    % Propagate the first six sigma points at t_k (et(k+1))
    [~,sigma_temp] = ode113(@(t,x) kepMOON_rhs(0,x), [et_span(k) et_span(k+1)], sigma_pts_km1(1:6,s),options);
    sigma_pts_k(1:6,s) = sigma_temp(end, :)';
    [~,sigma_temp] = ode113(@(t,x) kepMOON_rhs(0,x), [et_span(k) et_span(k+1)], sigma_pts_km1(1:6,n+s),options);
    sigma_pts_k(1:6,n+s) = sigma_temp(end, :)';
    % Sigma points relative to LATLON - no dynamics
    sigma_pts_k(7:8,s)   = sigma_pts_km1(7:8,s);
    sigma_pts_k(7:8,n+s) = sigma_pts_km1(7:8,n+s);

    % Project propagated sigma points onto the measurements sub-space through
    % the measurement model
    RNG_temp = antenna_data_UKF(sigma_pts_k(7,s), sigma_pts_k(8,s), et_span(k+1), sigma_pts_k(1:6,s));
    gamma_pts_k(4,s)   = RNG_temp;
    gamma_pts_k(1:3,s) = sigma_pts_k(1:3,s); 
    
    RNG_temp = antenna_data_UKF(sigma_pts_k(7,n+s), sigma_pts_k(8,n+s), et_span(k+1), sigma_pts_k(1:6,n+s));
    gamma_pts_k(4,n+s)   = RNG_temp;
    gamma_pts_k(1:3,n+s) = sigma_pts_k(1:3,n+s); 
    
end

% --- Apply UT weigths to sigma and gamma points ---

% Calculate mean
for l = 1:(2*n+1)
    % Calculate the predicted mean state (sample mean)
    x_kp = x_kp + W_mean(l) * sigma_pts_k(:,l);
    % Calculate the predicted measurements (sample mean)
    y_kp = y_kp + W_mean(l) * gamma_pts_k(:,l);
end

% Save predicted measurements
y_pred(k) = y_kp(4);

% Calculate covariance
for l = 1:(2*n+1)
    % Calculate the predicted covariance (sample covariance)
    P_kp   = P_kp    + W_cov(l) * (sigma_pts_k(:,l) - x_kp)*(sigma_pts_k(:,l) - x_kp)';
    % Calculate the predicted covariance for measurements
    P_k_ee = P_k_ee  + W_cov(l) * (gamma_pts_k(:,l) - y_kp)*(gamma_pts_k(:,l) - y_kp)';
    % Calculate the cross covariance
    P_k_xy = P_k_xy  + W_cov(l) * (sigma_pts_k(:,l) - x_kp)*(gamma_pts_k(:,l) - y_kp)';
end

% Compute Kalman Gain
K_k = P_k_xy * (P_k_ee\eye(size(R)));

% Update step
x_k = x_kp + K_k*(Ym(:,k+1) - y_kp);
P_k = P_kp - K_k*P_k_ee*K_k';

% Store sequential updates of the UKF
x_aug          = [x_aug, x_k];
P_aug(:,:,k+1) = P_k;

% Update new k-1 variables to k
P_km1 = P_k;
x_km1 = x_k;

% Re-allocate the variables to be calculated in the next step
P_kp    = zeros(n,n);
P_k_ee  = R;
P_k_xy  = zeros(n,4);
x_kp    = zeros(n,1);
y_kp    = zeros(4,1);

end



end

