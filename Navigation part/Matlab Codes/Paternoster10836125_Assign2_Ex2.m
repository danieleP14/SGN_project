clc; clearvars; cspice_kclear; close all

% -------------------------------------------------------------------------
% Spacecraft Guidance & Navigation (2024/2025)
% Assignment #2, Exercise 2
% Author: Daniele Paternoster
% -------------------------------------------------------------------------
% NOTE 
% - requires to have "sgp4" folder in same position.
% - requires the "assignment02.tm" meta kernel in the same position
% - assignment02.tm requires also "kernels\" folder in the same position
% - requires "sgp4" folder with utilities in the same position
% - parts which used antenna_pointing.p to validate my personal function 
%   results are commented (they were there just as a check)
% -------------------------------------------------------------------------

%% Data
% Paths of sgp4 utilities
addpath('sgp4');

% Upload meta kernel to the pool (contains all needed kernels in \kernels)
cspice_furnsh('assignment02.tm')

% parameters definition for SGP4
typerun    = 'u';  % user-provided inputs to SGP4 Matlab function
opsmode    = 'a';  % afspc approach ('air force space command')
whichconst =  72;  % WGS72 constants (radius, gravitational parameter)

% TLE of SMOS ( tle\36036.3le file)
longstr1 = '1 36036U 09059A   24323.76060260  .00000600  00000-0  20543-3 0  9995';
longstr2 = '2 36036  98.4396 148.4689 0001262  95.1025 265.0307 14.39727995790658';

% Conversion from arcseconds to radians (for ddpsi and ddeps)
arcsec2rad = pi / (180*3600);

%% Point 1: Compute visibility windows of SMOS from KOUROU, TROLL, SVALBARD 

% --- Find reference epoch t_ref in ET ---

% Extract SMOS record from the TLE 
SMOSrec = twoline2rv(longstr1, longstr2, typerun, 'e', opsmode, whichconst);

% Find the date from the Julian date using invjday routine (sgp4\invjday.m)
[year,mon,day,hr,Min,sec] = invjday(SMOSrec.jdsatepoch, SMOSrec.jdsatepochf);

% Print the reference epoch of the TLE in calendar format
SMOS_refepoch_str = sprintf('%d-%02d-%02dT%02d:%02d:%02.6f', [year,mon,day,hr,Min,sec]);

% Find the reference epoch ot the TLE in ephemeris time (using spice str2et
% routine)
SMOS_refepoch_et = cspice_str2et(SMOS_refepoch_str);

% --- Find the osc. orbital elements of SMOS at ref. epoch of the TLE ---

% Use the sgp4 to obtain the state (r,v) @ref.epoch expressed in TEME
[SMOS_rec, r_teme, v_teme] = sgp4(SMOSrec, 0.0);

% Use the cspice routine oscelt to find the oscl. orb. els (input in TEME) 
[SMOS_oscelt_ref] = cspice_oscelt([r_teme;v_teme], SMOS_refepoch_et, SMOSrec.mu);

% Print the found oscl. orbital elements for SMOS at the et_reference of
% TLE
fprintf('---------------------------------------------------------------\n')
fprintf('*** Osculating orbital elements for SMOS at TLE @t_ref \n***')
fprintf('SMA   [km]:  %.5f \n', SMOS_oscelt_ref(1)/(1-SMOS_oscelt_ref(2)));
fprintf('ECC   [-]:  %.8f \n', SMOS_oscelt_ref(2));
fprintf('INC  [deg]: %.5f \n', SMOS_oscelt_ref(3)*cspice_dpr());
fprintf('RAAN [deg]: %.5f \n', SMOS_oscelt_ref(4)*cspice_dpr());
fprintf('ARGP [deg]: %.5f \n', SMOS_oscelt_ref(5)*cspice_dpr());
fprintf('M.AN [deg]: %.5f \n', SMOS_oscelt_ref(6)*cspice_dpr());
fprintf('---------------------------------------------------------------\n')

% --- Compute the SMOS state in ECI frame at the TLE ref. epoch ---

% Define "fake" acceleration (needed by the teme2eci routine)
a_teme = [0;0;0];

% Nutation parameters obtained from:
% (https://celestrak.org/SpaceData/EOP-Last5Years.txt)
% EPOCH used: 2024/11/18 (reference epoch of the TLE)
ddpsi_refepoch = -0.114752*arcsec2rad;
ddeps_refepoch = -0.007529*arcsec2rad;

% Centuries from 2000-01-01T00:00:00.00 TDT (for precession)
ttt_refepoch = cspice_unitim(SMOS_refepoch_et, 'ET', 'TDT')/cspice_jyear()/100;

% Use the teme2eci routine (sgp4\teme2eci.m) to compute the state of SMOS
% in the ECI frame at the reference epoch of TLE 
[r_eci, v_eci, ~] = teme2eci(r_teme, v_teme, a_teme, ...
    ttt_refepoch, ddpsi_refepoch, ddeps_refepoch);

% Define the initial and final epochs for the visibility window as strings (UTC)
init_epoch_str = '2024-11-18T20:30:00.0000';
fin_epoch_str  = '2024-11-18T22:15:00.0000';

% Convert UTC formats to ET for initial and final epoch.
et0 = cspice_str2et(init_epoch_str);
etf = cspice_str2et(fin_epoch_str);

% --- Compute the state of SMOS in ECI @ the initial epoch of the
% visibility window --- 

% Set the propagation tolerances for the propagation of the dynamics
options = odeset('RelTol',1e-13, 'AbsTol',1e-20);

% Propagate the state of SMOS from reference epoch of TLE (et_ref) to the initial
% epoch of the visibility window (et_0).
[~, xx_eci] = ode113(@(t,x) kepEARTH_rhs(t,x,false), [SMOS_refepoch_et et0], [r_eci; v_eci], options);

% Get the cartesian state in ECI coordinates for SMOS @et0 (last element)
xSMOS_eci_et0 = xx_eci(end,1:6);

% Calculate number of elements which are equispaced and between et0 and
% etf
n_el = round((etf-et0)/30.0)+1;

% Define the visibility window span from et0 to etf
et_window_span = linspace(et0, etf, n_el);

% Propagate the state from et0 to etf
[~, xxSMOS_eci] = ode113(@(t,x) kepEARTH_rhs(t,x,false), et_window_span, xSMOS_eci_et0, options);

% Use correct frequencies of acquisition (60s for KOU and SVAL, 30s for TROLL)
et_KOU = et_window_span(1:2:end);
et_TROLL = et_window_span;
et_SVAL = et_window_span(1:2:end);

% Use correct frequencies of acquisition (60s for KOU and SVAL, 30s for TROLL)
xxSMOS_eciKOU = xxSMOS_eci(1:2:end,:);
xxSMOS_eciTROLL = xxSMOS_eci;
xxSMOS_eciSVAL = xxSMOS_eci(1:2:end,:);

% Retrieve ideal antenna measurements from Kourou, Troll and Svalbard GSs
% Commented lines are relative to the lab4 function to validate results

[AZ_KOU, EL_KOU, RANGE_KOU] = antenna_data('KOUROU', et_KOU, xxSMOS_eciKOU);
%[AZ_KOUvv, EL_KOUvv, RANGE_KOUvv] = antenna_pointing('KOUROU', et_KOU, xxSMOS_eciKOU');

[AZ_TROLL, EL_TROLL, RANGE_TROLL] = antenna_data('TROLL', et_TROLL, xxSMOS_eciTROLL);
%[AZ_TROLLvv, EL_TROLLvv, RANGE_TROLLvv] = antenna_pointing('TROLL', et_TROLL, xxSMOS_eciTROLL');

[AZ_SVAL, EL_SVAL, RANGE_SVAL] = antenna_data('SVALBARD', et_SVAL, xxSMOS_eciSVAL);
%[AZ_SVALvv, EL_SVALvv, RANGE_SVALvv] = antenna_pointing('SVALBARD', et_SVAL, xxSMOS_eciSVAL');

% Filter-out the unfeasible elements for the GSs - basing on minimum
% detectable elevation
index_KOU = EL_KOU > cspice_rpd*6;
index_TROLL = EL_TROLL > 0;
index_SVAL = EL_SVAL > cspice_rpd*8;

% Compute the epochs for the visibility windows
et_KOU_visible   = et_KOU(index_KOU);
et_TROLL_visible = et_TROLL(index_TROLL);
et_SVAL_visible  = et_SVAL(index_SVAL);

% Convert the ET visibility (initial & final) epochs in calendar format UTC
UTC_KOU_visibile_i = cspice_et2utc(et_KOU_visible(1), 'C', 10);
UTC_KOU_visibile_f = cspice_et2utc(et_KOU_visible(end), 'C', 10);

UTC_TROLL_visibile_i = cspice_et2utc(et_TROLL_visible(1), 'C', 10);
UTC_TROLL_visibile_f = cspice_et2utc(et_TROLL_visible(end), 'C', 10);

UTC_SVAL_visibile_i = cspice_et2utc(et_SVAL_visible(1), 'C', 10);
UTC_SVAL_visibile_f = cspice_et2utc(et_SVAL_visible(end), 'C', 10);

% Overall visibilities windows in calendar format UTC - KOU; TROLL; SVAL
fprintf('Kourou visibility starts at %s and ends at %s \n', UTC_KOU_visibile_i, UTC_KOU_visibile_f)
fprintf('Troll visibility starts at %s and ends at %s \n', UTC_TROLL_visibile_i, UTC_TROLL_visibile_f)
fprintf('Svalbard visibility starts at %s and ends at %s \n', UTC_SVAL_visibile_i, UTC_SVAL_visibile_f)

% --- Plot the Az and El (deg) ideal profiles for all the stations --- 

figure

t = tiledlayout(2,3, 'TileSpacing','compact','Padding','compact');
nexttile;
plot(cspice_dpr*AZ_KOU(index_KOU),cspice_dpr*EL_KOU(index_KOU),'+','MarkerSize',7,'LineWidth',1.5)
xlabel('AZ [deg]'); ylabel('EL [deg]')
title('KOUROU - AZEL')
grid on;
ylim([5 25]);

nexttile;
plot(cspice_dpr*AZ_TROLL(index_TROLL), cspice_dpr*EL_TROLL(index_TROLL),'+','MarkerSize',7,'LineWidth',1.5)
xlabel('AZ [deg]'); ylabel('EL [deg]')
grid on;
title('TROLL - AZEL')
ylim([0 7]);

nexttile;
plot(cspice_dpr*AZ_SVAL(index_SVAL), cspice_dpr*EL_SVAL(index_SVAL),'+','MarkerSize',7,'LineWidth',1.5)
xlabel('AZ [deg]'); ylabel('EL [deg]')
grid on;
title('SVAL - AZEL')
ylim([0 70]);

nexttile;
plot(RANGE_KOU(index_KOU),'o','MarkerSize',2,'Color','r','LineWidth',1.2)
xlabel('$n^{th}$acq'); ylabel('$\rho$ [km]')
title('KOUROU - RANGE')
grid on;

nexttile;
plot(RANGE_TROLL(index_TROLL),'o','MarkerSize',2,'Color','r','LineWidth',1.2)
xlabel('$n^{th}$acq');ylabel('$\rho$ [km]')
grid on;
title('TROLL - RANGE')

nexttile;
plot(RANGE_SVAL(index_SVAL),'o','MarkerSize',2,'Color','r','LineWidth',1.2)
xlabel('$n^{th}$acq'); ylabel('$\rho$ [km]')
grid on;
title('SVAL - RANGE')

%% Point 2: Simulate antenna measurements (KOU, TROLL, SVAL)

% --- KOUROU measurements (Az & El in radians) ---

% Define Kourou struct
KOUstr.name = 'KOUROU';
KOUstr.acqSTime = 60;
KOUstr.minEL = cspice_rpd*6;
KOUstr.noiseMatrix = diag([cspice_rpd*125e-3, cspice_rpd*125e-3, 0.01].^2);

% Generate measurements for Kourou
[measAZ_KOU, measEL_KOU, measRNG_KOU, mean_meas_KOU, index_vis_KOU] = ...
    antenna_measurements(KOUstr, SMOS_rec, SMOS_refepoch_et, ...
    et_KOU_visible, [ddpsi_refepoch ddeps_refepoch]);

% After the measurements simulation --> the visibility window changes
et_KOU_visible_meas = et_KOU_visible(index_vis_KOU);

% --- TROLL measurements (minimum EL and Az,El noises in radians) ---

% Define Troll struct
TROLLstr.name = 'TROLL';
TROLLstr.acqSTime = 30;
TROLLstr.minEL = 0;
TROLLstr.noiseMatrix = diag([cspice_rpd*125e-3, cspice_rpd*125e-3, 0.01].^2);

% Define Troll measurements
[measAZ_TROLL, measEL_TROLL, measRNG_TROLL, mean_meas_TROLL, index_vis_TROLL] = ...
    antenna_measurements(TROLLstr, SMOS_rec, SMOS_refepoch_et, ...
    et_TROLL_visible, [ddpsi_refepoch ddeps_refepoch]);

% After the measurements simulation --> the visibility window changes
et_TROLL_visible_meas = et_TROLL_visible(index_vis_TROLL);

% --- SVAL measurements (Az & El in radians) --- 

% Define Sval struct
SVALstr.name = 'SVALBARD';
SVALstr.acqSTime = 60;
SVALstr.minEL = cspice_rpd*8;
SVALstr.noiseMatrix = diag([cspice_rpd*125e-3, cspice_rpd*125e-3, 0.01].^2);

% Define Svalbard measurements
[measAZ_SVAL, measEL_SVAL, measRNG_SVAL, mean_meas_SVAL, index_vis_SVAL] = ...
    antenna_measurements(SVALstr, SMOS_rec, SMOS_refepoch_et, ...
    et_SVAL_visible, [ddpsi_refepoch ddeps_refepoch]);

% After the measurements simulation --> the visibility window changes
et_SVAL_visible_meas = et_SVAL_visible(index_vis_SVAL);

% AZ-EL plots: the off-set between meas and ideal is not due to random
% Gaussian noise but also use of sgp4 and kep propagation in the two cases.

% ---------------------------PLOT MEASUREMENTS-----------------------------
% --- KOU measurements ---

fontOptions(12);

figure('Position', [100, 100, 1100, 700]);  % Resize figure
tK = tiledlayout(3,3, 'TileSpacing','compact','Padding','compact');
sgtitle('Simulated measurements')

nexttile(tK)
plot(cspice_dpr*measEL_KOU, '+')
hold on
grid minor
plot(cspice_dpr*EL_KOU(index_KOU), '*','MarkerSize',3);
plot(cspice_dpr*mean_meas_KOU(2,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('EL [deg]')
yline(KOUstr.minEL*cspice_dpr, 'r--')
lgd = legend('Measured', 'Real','Mean','$EL_{min}$','Orientation','horizontal');
lgd.FontSize = 10;

nexttile(tK)
plot(cspice_dpr*measAZ_KOU, '+')
hold on
grid minor
plot(cspice_dpr*AZ_KOU(index_KOU), '*','MarkerSize',3);
plot(cspice_dpr*mean_meas_KOU(1,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('AZ [deg]')

nexttile(tK)
plot(measRNG_KOU, '+')
hold on
grid minor
plot(RANGE_KOU(index_KOU), '*','MarkerSize',3);
plot(mean_meas_KOU(3,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('$\rho$ [km]')

annotation('textbox', [0.03, 0.7, 0.1, 0.1], 'String', 'KOUROU', 'EdgeColor', 'none','Rotation',90,'FontSize',12,'Color','k');

% --- TROLL measurements ---

nexttile(tK)
plot(cspice_dpr*measEL_TROLL, '+')
hold on
grid minor
plot(cspice_dpr*EL_TROLL(index_TROLL), '*','MarkerSize',3);
plot(cspice_dpr*mean_meas_TROLL(2,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('EL [deg]')
yline(TROLLstr.minEL*cspice_dpr, 'r--')
ylim([-1 7]);

nexttile(tK)
plot(cspice_dpr*measAZ_TROLL, '+')
hold on
grid minor
plot(cspice_dpr*AZ_TROLL(index_TROLL), '*','MarkerSize',3);
plot(cspice_dpr*mean_meas_TROLL(1,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('AZ [deg]')

nexttile(tK)
plot(measRNG_TROLL, '+')
hold on
grid minor
plot(RANGE_TROLL(index_TROLL), '*','MarkerSize',3);
plot(mean_meas_TROLL(3,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('$\rho$ [km]')

annotation('textbox', [0.03, 0.4, 0.1, 0.1], 'String', 'TROLL', 'EdgeColor', 'none','Rotation',90,'FontSize',12,'Color','k');

% --- SVAL plots ---

nexttile(tK)
plot(cspice_dpr*measEL_SVAL, '+')
hold on
grid minor
plot(cspice_dpr*EL_SVAL(index_SVAL), '*','MarkerSize',3);
plot(cspice_dpr*mean_meas_SVAL(2,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('EL [deg]')
yline(SVALstr.minEL*cspice_dpr, 'r--')
ylim([1 65]);

nexttile(tK)
plot(cspice_dpr*measAZ_SVAL, '+')
hold on
grid minor
plot(cspice_dpr*AZ_SVAL(index_SVAL), '*','MarkerSize',3);
plot(cspice_dpr*mean_meas_SVAL(1,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('AZ [deg]')

nexttile(tK)
plot(measRNG_SVAL, '+')
hold on
grid minor
plot(RANGE_SVAL(index_SVAL), '*','MarkerSize',3);
plot(mean_meas_SVAL(3,:), 'o');
xlabel('$n^{th} \; acq.$'); ylabel('$\rho$ [km]')

annotation('textbox', [0.03, 0.06, 0.1, 0.1], 'String', 'SVALBARD', 'EdgeColor', 'none','Rotation',90,'FontSize',12,'Color','k');

%% Point 3: mimimum variance Navigation solution 

% define sqrt of weight matrix (for lsqnonlin)
W_mKOU = inv(sqrtm(KOUstr.noiseMatrix));
W_mSVAL = inv(sqrtm(SVALstr.noiseMatrix));
W_mTROLL = inv(sqrtm(TROLLstr.noiseMatrix));

% -- Minimum-variance estimates ---

% Lsqnonlin options
options_lsqnonlin = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt','Display', 'iter', ...
   'FiniteDifferenceType','central');

% Flag for J2 (run filter without J2 model on Kourou measurements and all
% GSs)
J2flag = false;

% --- Batch filter with KOUROU measurements - pure keplerian dynamics;

% Kourou only, pure Keplerian - COST FUNCTION
fun_lsqKOU = @(x) cost(x,et0,et_KOU_visible_meas,[measAZ_KOU measEL_KOU measRNG_KOU],  W_mKOU, KOUstr,J2flag);

% Minimum variance estimate
[est_x_KOU,resnorm_KOU,residual_KOU, exitflag_KOU,~,~,jacobian_KOU] = lsqnonlin(fun_lsqKOU,xSMOS_eci_et0, [], [], options_lsqnonlin);

% ---- Post process the obtained quantites ---
Jacobian_KOU = full(jacobian_KOU);
P_KOU_ls = resnorm_KOU/(length(residual_KOU)-length(xSMOS_eci_et0)) .*inv(Jacobian_KOU'*Jacobian_KOU);

% Square root of the trace - KOU
SQRT_TRACE_POS_KOU = sqrt(trace(P_KOU_ls(1:3,1:3)));
SQRT_TRACE_VEL_KOU = sqrt(trace(P_KOU_ls(4:6,4:6)));

% --- Batch filter with ALL (KOU, TROLL, SVAL) measurements - pure keplerian dynamics

% All groundstations, pure Keplerian - COST FUNCTION
fun_lsqALL = @(x)  [cost(x,et0,et_KOU_visible_meas,[measAZ_KOU measEL_KOU measRNG_KOU], W_mKOU, KOUstr,J2flag);
                       cost(x,et0,et_TROLL_visible_meas,[measAZ_TROLL measEL_TROLL measRNG_TROLL], W_mTROLL, TROLLstr,J2flag);
                       cost(x,et0,et_SVAL_visible_meas,[measAZ_SVAL measEL_SVAL measRNG_SVAL], W_mSVAL, SVALstr,J2flag);
                  ];

% Minimum variance estimate
[est_x_ALL,resnorm_ALL,residual_ALL,exitflag_ALL,~,~,jacobian_ALL] = lsqnonlin(fun_lsqALL,xSMOS_eci_et0, [], [], options_lsqnonlin);

% ---- Post process the obtained quantites ---
Jacobian_ALL = full(jacobian_ALL);
P_ALL_ls = resnorm_ALL/(length(residual_ALL)-length(xSMOS_eci_et0)).*inv(Jacobian_ALL'*Jacobian_ALL);

% Square root of the trace - ALL
SQRT_TRACE_POS_ALL = sqrt(trace(P_ALL_ls(1:3,1:3)));
SQRT_TRACE_VEL_ALL = sqrt(trace(P_ALL_ls(4:6,4:6)));

% --- Batch filter with ALL (KOU, TROLL, SVAL) measurements - with
% keplerian and J2 ---

% Set J2 flag to true to account for the J2 perturbation in the dynamical
% model of the filter
J2flag = true;

% All groundstations, Keplerian with J2 - COST FUNCTION
fun_lsqALLwJ2 = @(x)  [cost(x,et0,et_KOU_visible_meas,[measAZ_KOU measEL_KOU measRNG_KOU], W_mKOU, KOUstr,J2flag);
                       cost(x,et0,et_TROLL_visible_meas,[measAZ_TROLL measEL_TROLL measRNG_TROLL], W_mTROLL, TROLLstr,J2flag);
                       cost(x,et0,et_SVAL_visible_meas,[measAZ_SVAL measEL_SVAL measRNG_SVAL], W_mSVAL, SVALstr,J2flag);
                  ];

% Minimum variance filter
[est_x_ALLwJ2,resnorm_ALLwJ2,residual_ALLwJ2,exitflag_ALLwJ2,~,~,jacobian_ALLwJ2] = lsqnonlin(fun_lsqALLwJ2,xSMOS_eci_et0, [], [], options_lsqnonlin);

% ---- Post process the obtained quantites ---

Jacobian_ALLwJ2 = full(jacobian_ALLwJ2);
P_ALLwJ2_ls = resnorm_ALLwJ2/(length(residual_ALLwJ2)-length(xSMOS_eci_et0)).*inv(Jacobian_ALLwJ2'*Jacobian_ALLwJ2);

% Square root of the trace - ALLwJ2
SQRT_TRACE_POS_ALLwJ2 = sqrt(trace(P_ALLwJ2_ls(1:3,1:3)));
SQRT_TRACE_VEL_ALLwJ2 = sqrt(trace(P_ALLwJ2_ls(4:6,4:6)));

%% --- Estimated CO-Variances for Keplerian elements elements --- 

% Perform conversion from cartesian uncertainty to a,i uncertainty using
% the car2kepEarth_v2 function 

% KOU
[kepKOU, dkep_KOU] = car2kepEarth_v2(est_x_KOU(1:3), est_x_KOU(4:6));
P_kepKOU = dkep_KOU * P_KOU_ls * dkep_KOU';
sigma_sma_KOU = sqrt(P_kepKOU(1,1));
sigma_inc_KOU = cspice_dpr*sqrt(P_kepKOU(2,2));

% ALL 
[kepALL, dkep_ALL] = car2kepEarth_v2(est_x_ALL(1:3), est_x_ALL(4:6));
P_kepALL = dkep_ALL * P_ALL_ls * dkep_ALL';
sigma_sma_ALL = sqrt(P_kepALL(1,1));
sigma_inc_ALL = cspice_dpr*sqrt(P_kepALL(2,2));

% ALL wJ2
[kepALLwJ2, dkep_ALLwJ2] = car2kepEarth_v2(est_x_ALLwJ2(1:3), est_x_ALLwJ2(4:6));
P_kepALLwJ2 = dkep_ALLwJ2 * P_ALLwJ2_ls * dkep_ALLwJ2';
sigma_SMA_ALLwJ2 = sqrt(P_kepALLwJ2(1,1));
sigma_inc_ALLwJ2 = cspice_dpr*sqrt(P_kepALLwJ2(2,2));

%% ---- PrintOut: Results for the 3 cases: KOU, ALL, ALLwJ2 ----

fprintf('--------KOUROU measurements & pure Keplerian dynamics--------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_KOU);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('sqrt(Tr(r)) %e [km]\n', SQRT_TRACE_POS_KOU)
fprintf('sqrt(Tr(v)) %e [km/s]\n', SQRT_TRACE_VEL_KOU)
fprintf('SMA std dev %e [km]\n', sigma_sma_KOU)
fprintf('INC std dev %e [deg]\n', sigma_inc_KOU)
fprintf('---------------------------------------------------------------\n');

fprintf('--------ALL measurements & pure Keplerian dynamics--------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_ALL);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('sqrt(Tr(r)) %e [km]\n', SQRT_TRACE_POS_ALL)
fprintf('sqrt(Tr(v)) %e [km/s]\n', SQRT_TRACE_VEL_ALL)
fprintf('SMA std dev %e  [km]\n', sigma_sma_ALL)
fprintf('INC std dev %e [deg]\n', sigma_inc_ALL)
fprintf('---------------------------------------------------------------\n');

fprintf('--------ALL measurements & J2 perturbation effect- -------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_ALLwJ2);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('sqrt(Tr(r)) %e [km]\n', SQRT_TRACE_POS_ALLwJ2)
fprintf('sqrt(Tr(v)) %e [km/s]\n', SQRT_TRACE_VEL_ALLwJ2)
fprintf('SMA std dev %e  [km]\n', sigma_SMA_ALLwJ2)
fprintf('INC std dev %e  [deg]\n', sigma_inc_ALLwJ2)
fprintf('---------------------------------------------------------------\n');

%% --------------- Point 4: Trade-off analysis ----------------------------
% Run lsqnonlin on all possible combinations:
% - all with J2
% - also single stations are considered 
% - order: KOUwJ2, TROLL, SVAL, KOU-SVAL, SVAL-TROLL, KOU-TROLL.
% - all runs organized in the just defined order in separate sections to
% allow re-run on single batch
%--------------------------------------------------------------------------
%% ---- KOU (with J2) ----

% Define cost function for Kou with J2
fun_lsqKOUwJ2 = @(x) cost(x,et0,et_KOU_visible_meas,[measAZ_KOU measEL_KOU measRNG_KOU], W_mKOU, KOUstr,true);

% Run lsqnonlin on defined cost function
[est_x_KOUwJ2,resnorm_KOUwJ2,residual_KOUwJ2,exitflag_KOUwJ2,~,~,jacobian_KOUwJ2] = lsqnonlin(fun_lsqKOUwJ2,xSMOS_eci_et0, [], [], options_lsqnonlin);

% --- Post processing of solution ---
Jacobian_KOUwJ2 = full(jacobian_KOUwJ2);
P_KOUwJ2_ls = resnorm_KOUwJ2/(length(residual_KOUwJ2)-length(xSMOS_eci_et0)).*inv(Jacobian_KOUwJ2.'*Jacobian_KOUwJ2);

[kepKOUwJ2, dkep_KOUwJ2] = car2kepEarth_v2(est_x_KOUwJ2(1:3), est_x_KOUwJ2(4:6));
P_kepKOUwJ2 = dkep_KOUwJ2 * P_KOUwJ2_ls * dkep_KOUwJ2';

%% ---- TROLL (with J2) ----

% Define cost function for Troll with J2
fun_lsqTROLL = @(x) cost(x,et0, et_TROLL_visible_meas,[measAZ_TROLL measEL_TROLL measRNG_TROLL], W_mTROLL, TROLLstr,true);

% Run lsqnonlin on defined cost function
[est_x_TROLL,resnorm_TROLL,residual_TROLL,exitflag_TROLL,~,~,jacobian_TROLL] = lsqnonlin(fun_lsqTROLL,xSMOS_eci_et0, [], [], options_lsqnonlin);

% --- Post processing of solution ---
Jacobian_TROLL= full(jacobian_TROLL);
P_TROLL_ls = resnorm_TROLL/(length(residual_TROLL)-length(xSMOS_eci_et0)).*inv(Jacobian_TROLL.'*jacobian_TROLL);

[kepTROLL, dkep_TROLL] = car2kepEarth_v2(est_x_TROLL(1:3), est_x_TROLL(4:6));
P_kepTROLL = dkep_TROLL * P_TROLL_ls * dkep_TROLL';

%% ---- SVAL (with J2) ----

% Define cost function for Sval with J2
fun_lsqSVAL = @(x) cost(x,et0, et_SVAL_visible_meas,[measAZ_SVAL measEL_SVAL measRNG_SVAL], W_mSVAL, SVALstr,true);

% Run lsqnonlin on defined cost function
[est_x_SVAL,resnorm_SVAL,residual_SVAL,exitflag_SVAL,~,~,jacobian_SVAL] = lsqnonlin(fun_lsqSVAL,xSMOS_eci_et0, [], [], options_lsqnonlin);

% --- Post processing of solution ---
Jacobian_SVAL= full(jacobian_SVAL);
P_SVAL_ls = resnorm_SVAL/(length(residual_SVAL)-length(xSMOS_eci_et0)).*inv(Jacobian_SVAL.'*jacobian_SVAL);

[kepSVAL, dkep_SVAL] = car2kepEarth_v2(est_x_SVAL(1:3), est_x_SVAL(4:6));
P_kepSVAL = dkep_SVAL * P_SVAL_ls * dkep_SVAL';

%% ---- KOU-SVAL (with J2) ----

% Define cost function for KouSval with J2
fun_lsqKOUSVAL = @(x) [cost(x,et0, et_KOU_visible_meas,[measAZ_KOU measEL_KOU measRNG_KOU], W_mKOU, KOUstr,true);
                       cost(x,et0, et_SVAL_visible_meas,[measAZ_SVAL measEL_SVAL measRNG_SVAL], W_mSVAL, SVALstr,true)
                       ];
  
% Run lsqnonlin on defined cost function
[est_x_KOUSVAL,resnorm_KOUSVAL,residual_KOUSVAL,exitflag_KOUSVAL,~,~,jacobian_KOUSVAL] = lsqnonlin(fun_lsqKOUSVAL,xSMOS_eci_et0, [], [], options_lsqnonlin);

% --- Post processing of solution ---
Jacobian_KOUSVAL= full(jacobian_KOUSVAL);
P_KOUSVAL_ls = resnorm_KOUSVAL/(length(residual_KOUSVAL)-length(xSMOS_eci_et0)).*inv(Jacobian_KOUSVAL.'*Jacobian_KOUSVAL);

[kepKOUSVAL, dkep_KOUSVAL] = car2kepEarth_v2(est_x_KOUSVAL(1:3), est_x_KOUSVAL(4:6));
P_kepKOUSVAL = dkep_KOUSVAL * P_KOUSVAL_ls * dkep_KOUSVAL';

%% ---- SVAL-TROLL (with J2)

% Define cost function for SvalTroll with J2
fun_lsqSVALTROLL = @(x) [ cost(x,et0, et_SVAL_visible_meas,[measAZ_SVAL measEL_SVAL measRNG_SVAL], W_mSVAL, SVALstr,true);
                          cost(x,et0, et_TROLL_visible_meas,[measAZ_TROLL measEL_TROLL measRNG_TROLL], W_mTROLL, TROLLstr,true)
                       ];

% Run lsqnonlin on defined cost function
[est_x_SVALTROLL,resnorm_SVALTROLL,residual_SVALTROLL,exitflag_SVALTROLL,~,~,jacobian_SVALTROLL] = lsqnonlin(fun_lsqSVALTROLL,xSMOS_eci_et0, [], [], options_lsqnonlin);

% --- Post processing of solution ---
Jacobian_SVALTROLL= full(jacobian_SVALTROLL);
P_SVALTROLL_ls = resnorm_SVALTROLL/(length(residual_SVALTROLL)-length(xSMOS_eci_et0)).*inv(Jacobian_SVALTROLL.'*jacobian_SVALTROLL);

[kepSVALTROLL, dkep_SVALTROLL] = car2kepEarth_v2(est_x_SVALTROLL(1:3), est_x_SVALTROLL(4:6));
P_kepSVALTROLL = dkep_SVALTROLL * P_SVALTROLL_ls * dkep_SVALTROLL';

%% ---- KOU-TROLL (with J2) ----

% Define cost function for KouTroll with J2
fun_lsqKOUTROLL = @(x) [cost(x, et0, et_KOU_visible_meas,[measAZ_KOU measEL_KOU measRNG_KOU], W_mKOU, KOUstr,true);
                        cost(x,et0, et_TROLL_visible_meas,[measAZ_TROLL measEL_TROLL measRNG_TROLL], W_mTROLL, TROLLstr,true)
                       ];

% Run lsqnonlin on defined cost function
[est_x_KOUTROLL,resnorm_KOUTROLL,residual_KOUTROLL,exitflag_KOUTROLL,~,~,jacobian_KOUTROLL] = lsqnonlin(fun_lsqKOUTROLL,xSMOS_eci_et0, [], [], options_lsqnonlin);

% --- Post processing of solution ---
Jacobian_KOUTROLL= full(jacobian_KOUTROLL);
P_KOUTROLL_ls = resnorm_KOUTROLL/(length(residual_KOUTROLL)-length(xSMOS_eci_et0)).*inv(Jacobian_KOUTROLL.'*jacobian_KOUTROLL);

[kepKOUTROLL, dkep_KOUTROLL] = car2kepEarth_v2(est_x_KOUTROLL(1:3), est_x_KOUTROLL(4:6));
P_kepKOUTROLL = dkep_KOUTROLL * P_KOUTROLL_ls * dkep_KOUTROLL';

%% ------------------SGP4 propagated state at t0---------------------------

% Propagation with SGP4 from t-ref to t0
[~,rteme_SMOSt0,vteme_SMOSt0] = sgp4(SMOS_rec, (et0-SMOS_refepoch_et)/60.0);

% Centuries from 2000-01-01T00:00:00.00 TDT for precession
ttt0 = cspice_unitim(et0, 'ET', 'TDT')/cspice_jyear()/100;

% Convert state from TEME to ECI
 [reci_SMOSt0, veci_SMOSt0] = ...
    teme2eci(rteme_SMOSt0, vteme_SMOSt0, [0.0;0.0;0.0],  ttt0, ddpsi_refepoch, ddeps_refepoch);

% Define the state vector of SMOS in ECI coordinates @ et_vw_span(k)
x_SMOS_eci_t0 = [reci_SMOSt0; veci_SMOSt0];

%--------------------------------------------------------------------------

%% Post-processing of trade-off
% Calculation of the estimated error using as truth the value from sgp4
% (not used in the report).

xx_TOFF = [est_x_KOUwJ2; est_x_TROLL; est_x_SVAL; est_x_KOUSVAL; est_x_SVALTROLL; est_x_KOUTROLL; est_x_ALLwJ2];
err_R = vecnorm(xx_TOFF(:,1:3) - reci_SMOSt0',2,2);
err_V = vecnorm(xx_TOFF(:,4:6) - veci_SMOSt0',2,2);


%% Print-out of all the results on the different combinations GSs

% Extract SMA sigma
sigma_SMA = [sqrt(P_kepKOUwJ2(1,1)), sqrt(P_kepTROLL(1,1)), sqrt(P_kepSVAL(1,1)), ...
    sqrt(P_kepKOUSVAL(1,1)), sqrt(P_kepSVALTROLL(1,1)), sqrt(P_kepKOUTROLL(1,1)), sqrt(P_kepALLwJ2(1,1))];

% Extract INC sigma
sigma_INC = [sqrt(P_kepKOUwJ2(2,2)), sqrt(P_kepTROLL(2,2)), sqrt(P_kepSVAL(2,2)), ...
    sqrt(P_kepKOUSVAL(2,2)), sqrt(P_kepSVALTROLL(2,2)), sqrt(P_kepKOUTROLL(2,2)), sqrt(P_kepALLwJ2(2,2))]*cspice_dpr;


fprintf('--------KOUROU measurements & J2 perturbed dynamics--------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_KOUwJ2);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('SMA std dev %e [km]\n', sigma_SMA(1))
fprintf('INC std dev %e [deg]\n', sigma_INC(1))
fprintf('---------------------------------------------------------------\n');

fprintf('--------TROLL measurements & J2 perturbed dynamics--------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_TROLL);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('SMA std dev %e [km]\n', sigma_SMA(2))
fprintf('INC std dev %e [deg]\n', sigma_INC(2))
fprintf('---------------------------------------------------------------\n');

fprintf('--------SVALBARD measurements & J2 perturbed dynamics--------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_SVAL);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('SMA std dev %e [km]\n', sigma_SMA(3))
fprintf('INC std dev %e [deg]\n', sigma_INC(3))
fprintf('---------------------------------------------------------------\n');

fprintf('-----------------COMBINATIONS OF GROUNDSTATIONS----------------\n');

fprintf('--------KOUROU and SVAL & J2 perturbed dynamics--------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_KOUSVAL);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('SMA std dev %e [km]\n', sigma_SMA(4))
fprintf('INC std dev %e [deg]\n', sigma_INC(4))
fprintf('---------------------------------------------------------------\n');

fprintf('----------SVALBARD and TROLL & J2 perturbed dynamics--------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_SVALTROLL);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('SMA std dev %e [km]\n', sigma_SMA(5))
fprintf('INC std dev %e [deg]\n', sigma_INC(5))
fprintf('---------------------------------------------------------------\n');

fprintf('------KOUROU and TROLL measurements & J2 perturbed dynamics------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_KOUTROLL);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('SMA std dev %e [km]\n', sigma_SMA(6))
fprintf('INC std dev %e [deg]\n', sigma_INC(6))
fprintf('---------------------------------------------------------------\n');

fprintf('-----------------ALL GROUNDSTATIONS----------------\n');

fprintf('------ALL measurements & J2 perturbed dynamics------\n');
fprintf('The estimated state x at %s epoch is: \n [', init_epoch_str);
fprintf('%f ', est_x_ALLwJ2);
fprintf('] [km; km/s] in ECI frame\n');
fprintf('SMA std dev %e [km]\n', sigma_SMA(7))
fprintf('INC std dev %e [deg]\n', sigma_INC(7))
fprintf('---------------------------------------------------------------\n');
%% Point 5: Long term analysis 
% Propagation of the state for n_dats --> compute the antenna measurements
% (ideal) from different GS sites

% Days of propagation
n_days = 3;

% Define time span of propagation
et_span_LONG = [et0 et0+3600*24*n_days];

% Propagate 
[etODE,xxLONGTERM] = ode113(@(t,x) kepEARTH_rhs(t,x,true), et_span_LONG, xSMOS_eci_et0, options);

% Compute ideal angles of observations and range from Kou, Sval and Troll
[AZ_KOU_LT, EL_KOU_LT, RNG_KOU_LT] = antenna_data('KOUROU', etODE, xxLONGTERM);
[AZ_TROLL_LT, EL_TROLL_LT, RNG_TROLL_LT]= antenna_data('TROLL', etODE, xxLONGTERM);
[AZ_SVAL_LT, EL_SVAL_LT, RNG_SVAL_LT] = antenna_data('SVALBARD', etODE, xxLONGTERM);

% Plots for long terms analysis
deltaET_hours = (etODE - et0)/(3600*24);

% Plot the long term Elevation 
figure
fontOptions(12);

tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
sgtitle({'Elevation angle - Long term from 2024 NOV 18 20:30:00.000UTC'})

nexttile
plot(deltaET_hours, cspice_dpr*EL_KOU_LT);
xlabel('[days]')
ylabel('EL [deg]')
grid minor
yline(6,'r--','LineWidth',2)
annotation('textbox', [0.03, 0.7, 0.1, 0.1], 'String', 'KOUROU', 'EdgeColor', 'none','Rotation',90,'FontSize',12,'Color','k');

nexttile
plot(deltaET_hours, cspice_dpr*EL_TROLL_LT);
xlabel('[days]')
ylabel('EL [deg]')
grid minor
yline(0,'r--','LineWidth',2)
annotation('textbox', [0.03, 0.4, 0.1, 0.1], 'String', 'TROLL', 'EdgeColor', 'none','Rotation',90,'FontSize',12,'Color','k');

nexttile
plot(deltaET_hours, cspice_dpr*EL_SVAL_LT);
xlabel('[days]')
ylabel('EL [deg]')
grid minor
yline(6,'r--','LineWidth',2)
annotation('textbox', [0.03, 0.06, 0.1, 0.1], 'String', 'SVALBARD', 'EdgeColor', 'none','Rotation',90,'FontSize',12,'Color','k');


%% functions

function f = kepEARTH_rhs(et,x,J2flag)
%--------------------------------------------------------------------------
% This function allows to calculate the RHS for the keplerian dynamics
% written in ECI frame. It can also add the J2 perturbation if specified by
% the flag.
% -------------------------------------------------------------------------
% INPUT:
% - et           [1]: time
% - x            [6x1]: orbital state of the S/C in ECI frame
% - J2flag       [bool]: if set to 'true' the function add the J2 acc.
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

% Get the gravitational parameter for Earth using spice routines
muE = cspice_bodvrd('EARTH', 'GM', 1);

% calculate the two sub-vectors of the rhs 
f(1:3) = v;
f(4:6) = -muE*r/rnorm^3;

% add J2 term if asked
if J2flag==true

    % Define J2 parameters
    J2 = 0.0010826269;
    Radii = cspice_bodvrd('EARTH','RADII',3);
    Re_mean = (Radii(1) + Radii(3))/2;

    % Define transformation from ECI to ECEF frame (in which acc. defined)
    rot_ECI2ECEF = cspice_pxform('J2000','ITRF93',et);
    
    % Transform the position from ECI to ECEF
    r_ecef = rot_ECI2ECEF*r;

    % Write the acceleration of J2 in ECEF
    a_j2_ECEF = 3/2*muE*J2*(r_ecef(1:3)/rnorm^3)*(Re_mean/rnorm)^2.*(5*(r_ecef(3)/rnorm)^2 - [1;1;3]);

    % Define transformation matrix from ECEF to ECI
    rot_ECEF2ECI = cspice_pxform('ITRF93','J2000',et);

    % Transform acceleration in ECI fram
    a_j2_ECI = rot_ECEF2ECI * a_j2_ECEF;
     
    f(4:6) = f(4:6) + a_j2_ECI;
end

end

function [AZ, EL, RANGE] = antenna_data(station, et_span, x_eci)

% -------------------------------------------------------------------------
% This function allows to obtain the ideal values of AZ, EL, RANGE 
% of a space object with respect to a groundstation ('KOROU', 'TROLL', 'SVALBARD'). 
% NOTE: visibility constraint are not applied here.
%
% INPUT:
% - station name [string]: (one between: 'TROLL', 'KOROU', 'SVALBARD')
% - et_span      [n]: ET span of interest
% - x_eci        [nx6]: orbital state of the space object in ECI frame
%
% OUTPUT:
% - AZ           [1xn] AZ vector
% - EL           [1xn] EL vector
% - RANGE        [1xn] RANGE vector
% 
% -------------------------------------------------------------------------

% Check for compatibility of stations
if ~strcmp(station,'TROLL') && ~strcmp(station,'KOUROU') && ~strcmp(station,'SVALBARD')
    error('Invalid station name!')
end

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
    r_station_eci = cspice_spkpos(station, et_k, 'J2000', 'NONE', 'EARTH');
    
    % Get ECI position of SC wrt to station
    r_sat_station_eci = x_eci(k,1:3)' - r_station_eci;
    
    % Transformation from ECI to topocentric frame
    ROT_ECI2TOPO = cspice_pxform('J2000', topocentric_frame, et_k);
    
    % Convert state into topocentric frame
    r_sat_station_topo = ROT_ECI2TOPO*r_sat_station_eci;
    
    % Obtain Azimuth and Elevation coordinates with spice reclat routine
    [RANGE(k), AZ(k),  EL(k)] = cspice_reclat(r_sat_station_topo(:,:));
    
end

end

function [AZ_meas, EL_meas, RANGE_meas, mean_meas, index_vis] = antenna_measurements(station_str, satrec, et_TLE_ref, et_vw, nutation_par)

% -------------------------------------------------------------------------
% This function allows to obtain the simulated measurements of AZ, EL, RANGE 
% of a space object with respect to a groundstation ('KOROU', 'TROLL', 'SVALBARD'). 
% The noise of the measurements is normally-distributed for AZ, EL, RANGE
% with nor correlation between same epoch of measurement nor along time
% correlation
%
% INPUT:
% - station_str  [struct]: contains main features of the station
% - satrec       [struct]: contains the record of the sat. from TLE
% - et_TLE_ref   [1x1]: reference epoch of the TLE
% - et_vw        [1xn]: ideal visibility window (based on ideal measurements)
% - nutation par [ddpsi ddeps]: nutation parameters @ ref epoch of the TLE
%
% OUTPUT:
% - AZ_meas      [mx1] AZ vector
% - EL_meas      [mx1] EL vector
% - RANGE_meas   [mx1] RANGE vector
% 
% NOTE:
% Visibility restriction are applied inside (based on minimum elevation).
% The input number of visibility epoch is n (based on ideal measruements),
% after applying the noise on the quantities the number of measruements
% could be less or equal (m is the number of simulated measurements)
% -------------------------------------------------------------------------

% Check for compatibility of stations
if ~strcmp(station_str.name,'TROLL') && ~strcmp(station_str.name,'KOUROU') && ~strcmp(station_str.name,'SVALBARD')
    error('Invalid station name!')
end

% Retrieve data of the station
stationName  = station_str.name;
min_EL       = station_str.minEL;
noiseCovMat  = station_str.noiseMatrix;

% Define covariance matrix for measurements simulations
R = noiseCovMat;

% Define the topocentric frame
topocentric_frame = [stationName, '_TOPO'];

% Extract nutation parameter from input vector
ddpsi = nutation_par(1);
ddeps = nutation_par(2);

% Loop over the epoch of the ideal visibility window
for k=1:length(et_vw)

    % Time since reference epoch in minutes
    et_since = (et_vw(k)-et_TLE_ref)/60.0;

    % Propagation with SGP4
    [~,rteme_sat,vteme_sat] = sgp4(satrec, et_since);

    % Centuries from 2000-01-01T00:00:00.00 TDT for precession
    ttt = cspice_unitim(et_vw(k), 'ET', 'TDT')/cspice_jyear()/100;

    % Convert state from TEME to ECI
     [reci_sat, veci_sat] = ...
        teme2eci(rteme_sat, vteme_sat, [0.0;0.0;0.0],  ttt, ddpsi, ddeps);
    
    % Define the state vector of SMOS in ECI coordinates @ et_vw_span(k)
    x_sat_eci = [reci_sat; veci_sat];

    % Get station position ECI
    x_station_eci = cspice_spkezr(stationName, et_vw(k), 'J2000', 'NONE', 'EARTH');
    
    % Get ECI position of SC wrt to station
    x_sat_station_eci = x_sat_eci - x_station_eci;
    
    % Transformation from ECI to topocentric frame
    ROT_ECI2TOPO = cspice_sxform('J2000', topocentric_frame, et_vw(k));
    
    % Convert state into topocentric frame
    x_sat_station_topo = ROT_ECI2TOPO*x_sat_station_eci;
    
    % Convert coordinates from rectangular to latitudinal coordinates using
    % spice routines - these are the expected value of the measurements
    [RANGE(k), AZ(k),  EL(k)] = cspice_reclat(x_sat_station_topo(1:3));

end

    % Expected mean measurements
    mean_meas = [AZ; EL; RANGE];

    % Measurement simulation for AZ, EL, RANGE
    measTOT = mvnrnd(mean_meas', R);

    % Extract the measurement components
    AZ_measTOT    = measTOT(:,1);
    EL_measTOT    = measTOT(:,2);
    RANGE_measTOT = measTOT(:,3);
  
    % Exclude not visible measurements - elevation constraint
    index_vis = EL_measTOT > min_EL;

    % Output only compatible measurements with minimum elevation detectable
    AZ_meas = AZ_measTOT(index_vis);
    EL_meas = EL_measTOT(index_vis);
    RANGE_meas = RANGE_measTOT(index_vis);

end

function residuals = cost(xk,et_k,et_meas,Ym,W_m,station_str,J2flag)
 
%--------------------------------------------------------------------------
% This function allows to calculate the objective function needed by
% lsqnonlin in order to solve a minimum-variance estimation problem.
%
% INPUT: 
% - xk              [6x1]: state to estimate (S/C state in ECI)
% - et_k              [1]: reference epoch to evaluate xk
% - et_meas         [nx1]: epochs at which the measurements are acquired
% - Ym              [nx3]: batch of measurements
% - W_m             [3x3]: square root of the inverse of the cov.matrix
% - station_str  [struct]: structure of the Groundstation
% - J2flag         [bool]: flag for the j2 model
% OUTPUT:
% - residuals      [nx3]: matrix which contains the residuals (3x1) at each
%                         time (n times)
% NOTE:
% Notice that the function must not compute the product of the residual,
% for the same reason the SQUARE root of the inverse of the covariance is
% used.
%--------------------------------------------------------------------------

% Extract station information
station_name = station_str.name;

% Number of measurements
n_meas = length(et_meas);

% Initialize variables
residuals = [];

% Propagation options
options = odeset('RelTol',1e-13, 'AbsTol',1e-20);

% Loop over the measurements
for k=1:n_meas
    
    % Propagate from reference epoch to measurement k-th epoch
    [~, xx_eci] = ode113(@(t,x) kepEARTH_rhs(t,x,J2flag), [et_k et_meas(k)], xk, options);
    
    % Generate ideal measurement
    [AZ, EL, RANGE] = antenna_data(station_name, et_meas(k) ,xx_eci(end,:));
   
    % Collect ideal measurements
    y_eci = [AZ, EL, RANGE];

    % Compute residual between angles
    res_AZEL = angdiff(Ym(k,1:2), y_eci(1:2));

    % Compute range residual
    res_RNG = Ym(k,3) - y_eci(3);

    % Collect residuals
    res = [res_AZEL, res_RNG];

    % Calculate the k-th residual (weighted)
    residuals(k,:) =  (W_m*(res)')';
end

end

function [kep, dkep] = car2kepEarth_v2(r,v)
%--------------------------------------------------------------------------
% This function allows to calculate:
% - keplerian elements (a,i) 
% - jacobian of the (a,i) with respect to cartesian state x = (r,v)
% Given the state (r,v) in cartesian inertial coordinates for the Earth
% -------------------------------------------------------------------------
% INPUT:
% - r            [3x1]: position in inertial coordinates
% - v            [3x1]: velocity in inertial coordinates
% OUTPUT:
% - kep          [2x1]: SMA (a) and inclination (i)
% - dkep         [2x6]: jacobian of a and i wrt x = (r,v)
% DEPENDECIES:
% - it requires to have the gm_de432.tpc kernel uploaded
%--------------------------------------------------------------------------

% Get the mu for the Earth
mu = cspice_bodvrd('EARTH', 'GM', 1);

% Allocate space for jacobian
dkep = zeros(2,6);

% Get norm of r and v
rnorm = norm(r);
vnorm = norm(v);

% calculate SMA, angular momentum and inclination
a = mu/(vnorm^2 - 2*mu/rnorm);
h = cross(r,v);
i = acos(h(3)/norm(h));

kep = [a,i];

% If requested calculate also the jacobian
if nargout > 1

    % Define components of position and velocity
    x = r(1); 
    y = r(2);
    z = r(3);
    vx = v(1);
    vy = v(2);
    vz = v(3);

    % First column of the jacobian (derivative of SMA wrt x)
    dkep(1,1:3) = (2*mu^2/(rnorm*vnorm^2 - 2*mu)^2 * r/rnorm)';
    dkep(1,4:6) = ((2*mu/((vnorm^2 - 2*mu/rnorm)^2))*v)';

    % Extract z component and norm from h
    hz = h(3);
    hnorm = norm(h);

    % Second row of the jacobian (derivative of i wrt x)
    dkep(2,:) = (1/sqrt(1 - (hz/hnorm)^2))/hnorm * (hz/hnorm^2 * [cross(v,h)'; -cross(r,h)'] - [vy; -vx; 0; -y; x; 0]);

% Symbolic expression for check of dkep(2,:) 
%     
%      grad = [ 
%                     -((vz*y - vy*z)*(z*vx^2 - vz*x*vx + z*vy^2 - vz*y*vy))/((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2);
%                      ((vz*x - vx*z)*(z*vx^2 - vz*x*vx + z*vy^2 - vz*y*vy))/((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2);
%                      ((2*vx*(vz*x - vx*z) + 2*vy*(vz*y - vy*z))*(vy*x - vx*y))/(2*((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2));
%                     -((vz*y - vy*z)*(vz*x^2 - vx*z*x + vz*y^2 - vy*z*y))/((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2);
%                      ((vz*x - vx*z)*(vz*x^2 - vx*z*x + vz*y^2 - vy*z*y))/((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2);
%                     -((2*x*(vz*x - vx*z) + 2*y*(vz*y - vy*z))*(vy*x - vx*y))/(2*((vy*x - vx*y)^2 + (vz*x - vx*z)^2 + (vz*y - vy*z)^2)^(3/2))
%             ];
% 
%     dkep(2,:) = -1/sqrt(1 - (h(3)/norm(h))^2) * grad';

end

end

function fontOptions(DefaultFontSize)
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaulttextinterpreter', 'latex');
    set(0, 'defaultAxesFontSize', DefaultFontSize)
end
