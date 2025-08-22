clc; clearvars; cspice_kclear; close all
%--------------------------------------------------------------------------
% Spacecraft Guidance & Navigation (2023/2024)
% Assignment #2, Exercise 1
% Author: Daniele Paternoster
%--------------------------------------------------------------------------
% Note: requires to have \kernels folder in the same position where it is 
% executed.
%--------------------------------------------------------------------------

% Upload kernels of interest for constants retrieval
cspice_furnsh('kernels\gm_de432.tpc'); %GM values
cspice_furnsh('kernels\pck00010.tpc'); % planetary constants (radii)

% Set default plot options
DEFAULT_FONT_SIZE = 12;
fontOptions(DEFAULT_FONT_SIZE);

% --- Data ---
% Reference initial state
r_i = [-0.011965533749906; -0.017025663128129];
v_i = [10.718855256727338;  0.116502348513671];

x_i = [r_i; v_i];

% Initial and final time
t_i = 1.282800225339865;
t_f = 9.595124551366348;

% Covariance of the initial state 
P_i = [1.041e-15 6.026e-17 5.647e-16 4.577e-15; 
       6.026e-17 4.287e-18 4.312e-17 1.855e-16;
       5.647e-16 4.312e-17 4.432e-16 1.455e-15;
       4.577e-15 1.855e-16 1.455e-15 2.822e-14];

% Extract Moon & Earth Kernel data
muE = cspice_bodvrd('399', 'GM', 1);
muM = cspice_bodvrd('301', 'GM', 1);

rE = cspice_bodvrd('399', 'RADII', 3);
rM = cspice_bodvrd('301', 'RADII', 3);

% compute mu parameter for the PBRFBP
mu = muM / (muE + muM);

% Problem Data
hi_d  = 167;
hi_a  = 100;
l_em  = 3.84405e5;
om_em = 2.66186135e-6;
om_s  = -9.25195985e-1;
rho   = 3.88811143e2;
m_s   = 3.28900541e5;

% Scaling Units
TU = 4.34256461 * 24 * 3600;
DU = l_em;
VU = 1.02454018;

% struct for PBRFBP
par.mu   = mu;
par.rho  = rho;
par.om_s = om_s;
par.m_s  = m_s;
par.rE   = rE(1);
par.rM   = rM(1);
par.DU   = DU;

%% 1): LinCov & UT propagation

% Elements for Propagation
n_el = 5;

% Define the time span of interest
tt = linspace(t_i, t_f, n_el);

% LinCov propagation
[xx, P] = LinCov_PBRFBP(x_i, P_i, tt, par);

% UT propagation
[Y_est, P_est] = UT_PBRFBP(x_i, P_i, tt, par);

%% Plot the covariance ellipses for the two methods (lincov and UT) at t_f

% Propagate reference trajectory (slightly more than t_f)
options = odeset('AbsTol',1e-12,'RelTol',1e-12);
[tt_plot, xplot] = ode78(@(t,x) PBRFBP_rhs(t,x,par), [t_i t_f + 0.0005], x_i, options);

figure

% LinCov ellipse
plot_ellipse(xx(1:2,end), P(1:2,1:2,end), 'r');
hold on

% UT ellipse                                           green color
plot_ellipse(Y_est(1:2,end), P_est(1:2,1:2,end), [0.4660 0.6740 0.1880]);

% Reference trajectory: blue color
plot(xplot(end-40:end,1), xplot(end-40:end,2),'LineStyle','--','LineWidth',1.4, 'Color', [0 0.4470 0.7410])

% Labels and legend
xlabel('$x$ [DU]'); ylabel('$y$ [DU]'); 
legend('Mean LinCov', 'Covariance LinCov', 'Mean UT', 'Covariance UT', 'Reference trajectory')

% Create a box near the mean points
axes('position',[.2 .2 .15 .15])
box on % put box around new pair of axes
plot(xplot(end-28:end-24,1), xplot(end-28:end-24,2), 'LineStyle','--','LineWidth',1.4, 'Color', [0 0.4470 0.7410])
hold on;
plot(xx(1,end), xx(2,end), '+', 'MarkerSize',7, 'LineWidth',2, 'Color','r');
plot(Y_est(1,end), Y_est(2,end), '+', 'MarkerSize',7, 'LineWidth',2, 'Color',[0.4660 0.6740 0.1880]);

grid minor
axis equal

%% 2): Monte carlo simulation ("truth") and comparison with previous methods

% Samples to simulate the MC
n_samples = 1000;

% Run the MonteCarlo
[mu_MC, P_MC, xx_samples] = monteCarlo(x_i, P_i, tt, par, n_samples);

%% Plot mean and ellipses at final time

% Create figure for the propagations covariance ellipses at t_f
figure

% LinCov ellipse
plot_ellipse(xx(1:2,end), P(1:2,1:2,end), 'r');
hold on

% UT ellipse
plot_ellipse(Y_est(1:2,end), P_est(1:2,1:2,end),  [0.4660 0.6740 0.1880]);

% MC ellipse 
plot_ellipse(mu_MC(1:2,end), P_MC(1:2,1:2,end), 'k', '--');

% MC samples
grayColor = [.7 .7 .7];
plot(squeeze(xx_samples(end, 1, :)), squeeze(xx_samples(end, 2, :)), ...
    'Marker','+', 'MarkerSize', 1.3, 'LineStyle', 'none', 'Color', grayColor);

xlabel('$x$ [DU]'); ylabel('$y$ [DU]'); 
legend('Mean LinCov', 'Covariance LinCov', 'Mean UT', 'Covariance UT', 'Mean MC', 'Covariance MC')

%% Max uncertainty on position and velocity
% The max uncertainty for position and velocity is calculated by taking the
% position 2x2 or velocity 2x2 submatrix, taking the square root of the 
% max eigenvalue and multiple it by try (it is an estimate of the 3sigma)

% Number of epochs considered
n_epochs = length(tt);

% Extract the 3sigma from each approach at each epoch

maxEigRV_LC = maxEig(cat(3, P_i, P), n_epochs);
maxEigRV_UT = maxEig(cat(3, P_i, P_est), n_epochs);
maxEigRV_MC = maxEig(cat(3, P_i, P_MC), n_epochs);

% Set marker size
mrkSize = 10;
lw = 1.5;

% Colors
blue = [0 0.4470 0.7410];
red  = [0.8500 0.3250 0.0980];
yellow = [0.9290 0.6940 0.1250];

% Plot the obtained results
figure; 

t1 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile
plot(maxEigRV_LC(1,:),"+", 'MarkerSize', mrkSize,'LineWidth',lw, 'LineStyle', '--','Color',blue)
hold on;
grid on;
grid minor;
plot(maxEigRV_UT(1,:),"x", 'MarkerSize', mrkSize,'LineWidth',lw,'LineStyle', '--','Color', red)
plot(maxEigRV_MC(1,:),"o", 'MarkerSize', mrkSize,'LineWidth',lw,'LineStyle', '--','Color', yellow)

xlabel('Epoch [-]');
ylabel('$3\sqrt{\max(\lambda_i(P_r))}$ [DU]')
xticks([1 2 3 4 5])
xticklabels({'$t_1$','$t_2$', '$t_3$', '$t_4$' , '$t_5$'})
legend("LC", "UT", "MC");

nexttile
plot(maxEigRV_LC(1,:)-maxEigRV_MC(1,:),"+", 'MarkerSize', mrkSize,'LineWidth',lw, 'LineStyle', '--','Color',blue)
hold on;
grid on;
plot(maxEigRV_UT(1,:)-maxEigRV_MC(1,:),"x", 'MarkerSize', mrkSize,'LineWidth',lw,'LineStyle', '--','Color', red)

xlabel('Epoch [-]');
ylabel('$3\sigma_r$ difference [DU]')
xticks([1 2 3 4 5])
xticklabels({'$t_1$','$t_2$', '$t_3$', '$t_4$' , '$t_5$'})
legend("LC-MC", "UT-MC");

%% Velocity uncertainty 

figure;

t2 = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
semilogy(tt,maxEigRV_LC(2,:),"+", 'MarkerSize', mrkSize, 'LineWidth',lw,'LineStyle', '--','Color',blue)
hold on
grid on;
grid minor;
semilogy(tt,maxEigRV_UT(2,:),"x", 'MarkerSize', mrkSize,'LineWidth',lw,'LineStyle', '--','Color',red)
semilogy(tt,maxEigRV_MC(2,:),"o", 'MarkerSize', mrkSize,'LineWidth',lw,'LineStyle', '--','Color',yellow)
% 
xticks([1 2 3 4 5])
xticklabels({'$t_1$','$t_2$', '$t_3$', '$t_4$' , '$t_5$'})

legend("LC", "UT", "MC");
xlabel("Epoch [-]")
ylabel('$3\sqrt{\max(\lambda_i(P_v))}$ [VU]')
 
nexttile
plot(maxEigRV_LC(2,:)-maxEigRV_MC(2,:),"+", 'MarkerSize', mrkSize,'LineWidth',lw, 'LineStyle', '--','Color',blue)
hold on;
grid on;
plot(maxEigRV_UT(2,:)-maxEigRV_MC(2,:),"x", 'MarkerSize', mrkSize,'LineWidth',lw,'LineStyle', '--','Color', red)

xlabel('Epoch [-]');
ylabel('$3\sigma_v$ difference [VU]')
xticks([1 2 3 4 5])
xticklabels({'$t_1$','$t_2$', '$t_3$', '$t_4$' , '$t_5$'})
legend("LC-MC", "UT-MC");


%% QQ plots - All components of the final samples with MC

xx_MCsamples_tf = squeeze(xx_samples(end,:,:));

DEFAULT_FONT_SIZE = 12;
fontOptions(DEFAULT_FONT_SIZE);
figure
subplot(2,2,1);
qqplot(xx_MCsamples_tf(1,:));
title('$r_x$')
grid on

subplot(2,2,2);
qqplot(xx_MCsamples_tf(2,:));
title('$r_y$')
grid on

subplot(2,2,3);
qqplot(xx_MCsamples_tf(3,:));
title('$v_x$')
grid on

subplot(2,2,4);
qqplot(xx_MCsamples_tf(4,:));
title('$v_y$')
grid on

%% Functions

function maxEigRV = maxEig(P, n_epochs) 
%------------------max eigenvalues of a series of covmatrices--------------
% This function calculates the value of 3 times the square root of the 
% maximum eigenvalue of the position associated covariance and velocity
% associated submatrices in P at the n epochs of interest
% INPUT:
% - P0          [4x4xn_epochs]: covariance on initial condition
% - tt          [Nx1]: temporal instants in which calculate quantities
% - par         [struct]: struct of parameters for the PBRFBP
% - n           [int]: number of samples to create
% OUTPUT:
% - mu_MC       [4xN]  : sample mean value from MC at each time instant
% - P_MC        [4x4xN]: sample covariance from MC at each time instant
% - xx_samples  [Nx4xn]: propagated samples at each time instant
% DEPENDENCIES:
% PBRFBP_rhs(t,x,par) for propagation of the dynamics
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Allocate memory for the [3sqrt(max(eigP_r)), 3sqrt(max(eigP_r))] at each
% epoch of interest (n_epoch),
maxEigRV = zeros(2,n_epochs);

% Loop over the epochs of interest
for k = 1:n_epochs
    % Extract submatrices at k-th poch
    P_r = P(1:2,1:2, k);
    P_v = P(3:4,3:4, k);
    
    % calculate eigenvalue of submatrices
    lambda_r = eig(P_r);
    lambda_v = eig(P_v);
    
    % Compute the estimator of the uncertainty
    maxEigRV(1,k) = 3*sqrt(max(lambda_r));
    maxEigRV(2,k) = 3*sqrt(max(lambda_v));
end

end 

function [mu_MC, P_MC, xx_samples] = monteCarlo(x0, P0, tt, par, n)
%------------------Monte Carlo simulation - algorithm----------------------
% This function performs a Monte Carlo simulation of the initial condition
% x0 with covariance P0 propagateed in the PBRFBP. The drawn samples are n.
% INPUT:
% - x0          [4x1]: initial condition
% - P0          [4x4]: covariance on initial condition
% - tt          [Nx1]: temporal instants in which calculate quantities
% - par         [struct]: struct of parameters for the PBRFBP
% - n           [int]: number of samples to create
% OUTPUT:
% - mu_MC       [4xN]  : sample mean value from MC at each time instant
% - P_MC        [4x4xN]: sample covariance from MC at each time instant
% - xx_samples  [Nx4xn]: propagated samples at each time instant
% DEPENDENCIES:
% PBRFBP_rhs(t,x,par) for propagation of the dynamics
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------
% Number of time instant of interest
n_times = length(tt)-1;

% Allocate space for variables of output
mu_MC   = zeros(4,n_times);
P_MC    = zeros(4,4,n_times);
xx_samples = zeros(n_times,4,n);

% Options for propagation
options = odeset('AbsTol',1e-12,'RelTol',1e-12); 

% Perform the propagation for each sample
for k = 1:n
    % Extract an initial condition from the distribution (Gaussian)
    x_i = mvnrnd(x0, P0, 1);
   
    % Propagate the initial condition
    [~, xx] = ode78(@(t,x) PBRFBP_rhs(t,x,par), tt, x_i, options);
   
    % Save the propagated state at each instant except t_i
    xx_samples(:,:,k) = xx(2:end,:);
end

%Compute sample means and sample covariances (at each epoch)
for k = 1:n_times
   mu_MC(:,k) = mean(xx_samples(k,:,:),3);
   P_MC(:,:,k) = cov(squeeze(xx_samples(k,:,:))');
end

end

function [xx, P] = LinCov_PBRFBP(x0, P0, tt, par)
%------------------LinCov propagation algorithm----------------------------
%
% This function propagates the initial state x0 and its initial covariance
% by means of Linearization, under the PBRFBP vector field assumption 
% for the dynamics adimensionalized as stated in "Topputo, 2013 - On optimal 
% two-impulse Earth–Moon transfers in a four-body model". 
% The function propagates the state along the time span tt: 
% - it can have only initial and final time --> return the values for tf, 
% - it can be a grid of n time epochs from ti to tf, in this case the values of state and uncertainty
% are returned for all epochs of interest except initial time ti
% 
%--------------------------------------------------------------------------
%
% INPUT: 
% x0         [4x1] initial state
% P0         [4x4] initial covariance of x0
% tt         [n] timespan of propagation (n>=2)
% par        [struct] struct of 7 fields for PBRFBP parameters
% 
% OUTPUT:
% xx         [4x(n-1)] propagated states along tt
% P          [4x4x(n-1)] propagated covariances along tt    
% 
% DEPENDENCIES:
% propagate_PBRFBP()
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Allocate space for the propagated covariances
P = NaN(4,4,length(tt)-1);

% Compute the STM at sample points (the output is N-1 values because t_i is
% not considered clearly)
[xx, STM] = propagate_PBRFBP(x0, tt, par);

% Loop over the sample points
for k = 1:length(tt)-1
    % Extract STM at time t_k 
    STM_k = STM(:,:,k);
    % Compute covariance at time t_k
    P(:,:,k) = STM_k * P0 * STM_k';
end


end

function [Y_est, P_est] = UT_PBRFBP(x0, P0, tt, par)

%----------------Uscented Transform propagation algorithm------------------
%
% This function propagates the initial state x0 and its initial covariance
% by means of the Unscented Transform, under the PBRFBP vector field assumption 
% for the dynamics adimensionalized as stated in "Topputo, 2013 - On optimal 
% two-impulse Earth–Moon transfers in a four-body model". 
% The function propagates the state along the time span tt: 
% - it can have only initial and final time --> return the values for tf, 
% - it can be a grid of  n time epochs from ti to tf, in this case the values 
% of state and uncertainty are returned for all epochs of interest.
% 
%--------------------------------------------------------------------------
%
% INPUT: 
% x0         [4x1] initial state
% P0         [4x4] initial covariance of x0
% tt         [1xn_t] timespan of propagation (n_t>=2), where  n_t: number 
%                    of time epochs.
% par        [struct] struct of 7 fields for PBRFBP parameters
% 
% OUTPUT:
% Y_est      [4x(n_t-1)]   sample mean along tt
% P_est      [4x4x(n_t-1)] sample covariances along tt    
% 
% DEPENDENCIES:
% propagate_PBRFBP()
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% ---- step 0 ---- Define UT main parameters

% UT parameters
alfa = 1;
beta = 2;
n = length(x0);
lambda = (alfa^2 - 1)*n;

% Allocation of variables (sigma points, evaluation of sigma points, weights)
sigma_pts = zeros(n,2*n+1);
Y_pts     = zeros(n,2*n+1,length(tt)-1); % propagated states to different final epochs
W_mean    = zeros(2*n+1,1);
W_cov     = W_mean;

% Allocation of output variables (sample mean and sample covariance
% at each time step, except ti)
Y_est = zeros(4,length(tt)-1);
P_est = zeros(4,4,length(tt)-1);

% Calculate square root of matrix
sqrtMAT   = sqrtm((n+lambda)*P0);

% ---- step 1 ---- Compute Sigma points and Weigths------------------------

% Definition of first sigma point - sigma0 (1)
sigma_pts(:,1) = x0;

% Definition of sigma points - sigma_i (2N) 
for k = 2:n+1
    sigma_pts(:,k)   = x0 + sqrtMAT(:,k-1);
    sigma_pts(:,n+k) = x0 - sqrtMAT(:,k-1);
end

W_mean(1) = lambda/(n+lambda);
W_cov(1)  = lambda/(n+lambda) + (1 - alfa^2 + beta);
W_mean(2:end) = 1/(2*(n+lambda));
W_cov(2:end) = W_mean(2:end);

% --- step 2 - propagate the sigma points through the dynamics-------------
 
options = odeset('AbsTol',1e-12,'RelTol',1e-12);
   
for k = 1:(2*n+1)
     [~, yy] = ode78(@(t,x) PBRFBP_rhs(t, x, par), tt, sigma_pts(:,k), options);
     % Extract all values except the one associated to t_i
     Y_pts(:,k,:) = yy(2:end,:)'; 
end

% --- step 3 - compute final weighted sample mean and covariance-----------

% For each epoch of interest
for k = 1:length(tt)-1
    % Weighted mean
    for l = 1:(2*n+1)
        Y_est(:,k) = Y_est(:,k) + W_mean(l) * Y_pts(:,l,k);
    end
    % Weighted covariance
    for l = 1:(2*n+1)
        P_est(:,:,k) = P_est(:,:,k) + W_cov(l) * (Y_pts(:,l,k) - Y_est(:,k))*(Y_pts(:,l,k) - Y_est(:,k))';
    end
end

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

[~, dOM4] = OM4(t,xx(1:2), par);
dotx = xx(3); doty = xx(4);

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
%
% OUTPUT:
% - rhs_var      [20x1]: evaluated RHS 
% 
% DEPENDENCIES:
% - PBRFBP_rhs(t,x(1:4),par)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Define position from augmented vector
x = yy(1); y = yy(2);

% Calculate the rhs for the PBRFBP dynamics
f_rhs =  PBRFBP_rhs(t, yy(1:4), par);

% Extract paramters from struct
om_s = par.om_s;  mu = par.mu; rho = par.rho; m_s = par.m_s;

% Hessian of the OM4 potential 
d2OM = [(mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - m_s/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) + (3*m_s*(2*x - 2*rho*cos(om_s*t))^2)/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*(2*mu + 2*x - 2)^2)/(4*((mu + x - 1)^2 + y^2)^(5/2)) - (3*(2*mu + 2*x)^2*(mu - 1))/(4*((mu + x)^2 + y^2)^(5/2)) + 1, (3*m_s*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)); 
       (3*m_s*(2*x - 2*rho*cos(om_s*t))*(2*y - 2*rho*sin(om_s*t)))/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) + (3*mu*y*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(5/2)) - (3*y*(2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(5/2)), (mu - 1)/((mu + x)^2 + y^2)^(3/2) - mu/((mu + x - 1)^2 + y^2)^(3/2) - m_s/((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2) + (3*m_s*(2*y - 2*rho*sin(om_s*t))^2)/(4*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(5/2)) - (3*y^2*(mu - 1))/((mu + x)^2 + y^2)^(5/2) + (3*mu*y^2)/((mu + x - 1)^2 + y^2)^(5/2) + 1
       ];

% Define A(t) matrix
AA = [zeros(2,2),     eye(2);
        d2OM   , [0,2;-2,0]];

% Calculate the rhs of the variational equation STM related
STM_dot = AA * reshape(yy(5:end), [4 4]);

% Assemble the full rhs and vectorize the STM_dot
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

% Extract parameters from the struct 
om_s = par.om_s;  mu = par.mu; rho = par.rho; m_s = par.m_s;

% Define position components from vector
x = xx(1); y = xx(2); 

% --- Define norm of distances
% S/C Earth distance
r1 = sqrt((x+mu)^2 + y^2);
% S/C Moon distance
r2 = sqrt((x+mu-1)^2 + y^2);
% S/C Sun distance
r3 = sqrt((x-rho*cos(om_s*t))^2 + (y-rho*sin(om_s*t))^2);

% Define the potential of the PBRFBP
OM4 = 0.5*(x^2 + y^2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu) + m_s/r3 - (m_s/rho^2) * (x*cos(om_s*t) + y*sin(om_s*t));

% Define the gradient of the PBRFBP
dOM4 = [ x - (mu*(2*mu + 2*x - 2))/(2*((mu + x - 1)^2 + y^2)^(3/2)) - ...
        (m_s*cos(om_s*t))/rho^2 + ((2*mu + 2*x)*(mu - 1))/(2*((mu + x)^2 + y^2)^(3/2)) - ...
        (m_s*(2*x - 2*rho*cos(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2));
         y - (m_s*sin(om_s*t))/rho^2 - (mu*y)/((mu + x - 1)^2 + y^2)^(3/2) - ...
        (m_s*(2*y - 2*rho*sin(om_s*t)))/(2*((x - rho*cos(om_s*t))^2 + (y - rho*sin(om_s*t))^2)^(3/2)) ...
        + (y*(mu - 1))/((mu + x)^2 + y^2)^(3/2)
        ];

end

function [xx, STM] = propagate_PBRFBP(x0, tt, par)

%------------------Propagate the PBRFBP dynamics---------------------------
%
% This function propagates the initial state x0 in PBRFBP vector field 
% assumption for the dynamics adimensionalized as stated in "Topputo, 2013 
% - On optimal two-impulse Earth–Moon transfers in a four-body model". The
% propagation outputs are the state and the STM at the epoch of interest
% The function propagates the state along the time span tt: 
% - it can have only initial and final time --> return the values for tf, 
% - it can be a grid of n time epochs from ti to tf, in this case the values 
% of state and uncertainty are returned for all epochs of interest except 
% initial time ti
% 
%--------------------------------------------------------------------------
%
% INPUT: 
% x0         [4x1] initial state
% tt         [n] timespan of propagation (n>=2)
% par        [struct] struct of 7 fields for PBRFBP parameters
% 
% OUTPUT:
% xx         [4x(n-1)]  orbital states along tt
% STM        [4x4x(n-1)] STMs along tt    
% 
% DEPENDENCIES:
% rhs_variational(t,y,par): propagation of the state/STM dynamics in PBRFBP
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------
    
% Output the STM and final state at all  times (expect first which is t_i)
STM = NaN(4,4,length(tt)-1);
xx  = NaN(4,length(tt)-1);

% Define the vectorized initial condition
y0 = [x0; reshape(eye(4), [16 1])];

% Set tolerances for propagation
options = odeset('AbsTol',1e-12,'RelTol',1e-12);

% Propagate the state/STM coupled dynamics
[~, yy] = ode78(@(t,y) rhs_variational(t,y,par), tt, y0, options);

% Extract the data of interest
for k = 1:length(tt)-1
    xx(:,k) = yy(k+1,1:4)';
    STM(:,:,k) = reshape(yy(k+1, 5:20), [4 4]);
end

end

function plot_ellipse(mu, P, color, ls)
%------------------Plot covariance ellipse routine-------------------------
%
% This function allows to plot the covariance ellipse of a bi-variate
% random variable with mean mu and covariance P assuming a gaussian
% bivariate distribution.
% 
%--------------------------------------------------------------------------
%
% INPUT: 
% mu         [2x1] mean value
% P          [2x2] covariance
% color      [char] expresses ellipse color
% ls         [char] expresses line style (OPTIONAL)
%
% OUTPUT:
%  plot of the ellipse 
%
% DEPENDENCIES:
% NONE
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Grid of angles for the ellipse
theta = linspace(0, 2*pi, 100);

% Eigen-decomposition of the covariance matrix
% V -> eigenvectors (ellipse direction of semiaxes in the x-y plane)
% D -> eigenvalues (size of the semimajor(minor axis)
[V, D] = eig(P);

% Define the length of the ellipse axes (3sigma)
axesLengths = 3*sqrt(diag(D));

% Define the ellipse points as product of eigenvectors (directions
% of axes) times the 3sigma lengths times the paramterization (cos(th)
% sin(th)
ellipsePoints = (V * diag(axesLengths)) * [cos(theta); sin(theta)];

% Shift ellipse to the mean
ellipsePoints(1, :) = ellipsePoints(1, :) + mu(1);
ellipsePoints(2, :) = ellipsePoints(2, :) + mu(2);

% If not in input, ls is set to normal solid line
if nargin < 4
    ls = '-';
end

% Plot the ellipse centered at its mean
plot(mu(1), mu(2), '+','MarkerSize',7, 'Color', color); % Plot mean
hold on
plot(ellipsePoints(1, :), ellipsePoints(2, :), 'Linestyle', ls, 'Color', color, 'LineWidth', 1);
axis equal;
grid on;

end

function fontOptions(DefaultFontSize)
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaulttextinterpreter', 'latex');
    set(0, 'defaultAxesFontSize', DefaultFontSize)
end