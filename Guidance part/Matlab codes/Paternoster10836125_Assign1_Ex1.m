clc; clearvars; close all
% -------------------------------------------------------------------------
% Spacecraft Guidance & Navigation (2023/2024)
% Assignment #1, Exercise 1 - Periodic Orbit
% Author: Daniele Paternoster
% -------------------------------------------------------------------------

% --- Data definition --- 

% Gravitational parameter for the Earth-Moon system
mu = 0.012150;

% Jacobi constant of the desired Halo
C = 3.09;

% Initial guess for the desired Halo 

xx0 = [1.068792441776
                    0 
       0.071093328515 
                    0
       0.319422926485
                    0];

% Plot labels font size
DEFAULT_FONT_SIZE = 13;

%% Point 1: determine coordinates of Li points in fixed rotational frame

% --- Collinear Libration points: L1, L2, L3 ---

% Graph the derivative wrt x coord. of the potential function U(x,0,0)
xx1 = linspace(-2,-mu-0.2,100);
xx2 = linspace(-mu+0.2,1+mu-0.05,100);
xx3 = linspace(1+mu+0.005,2,100);

% Define the dU(x,0,0)/dx as function handle
fun = @(x) dUdx00(x,mu);

% The graph of dU(x,0,0)/dx allows to understand where the zeros are located
figure
fontOptions(DEFAULT_FONT_SIZE);
plot(xx1, dUdx00(xx1, mu),"LineWidth",1.3,'Color',[0 0.4470 0.7410],'HandleVisibility','off');
hold on;
plot(xx2, dUdx00(xx2, mu),"LineWidth",1.3,'Color',[0 0.4470 0.7410],'HandleVisibility','off');
plot(xx3, dUdx00(xx3, mu),"LineWidth",1.3,'Color',[0 0.4470 0.7410],'HandleVisibility','off');
xline(1-mu,'LineStyle','--','HandleVisibility','off')
xline(-mu,'LineStyle','--','HandleVisibility','off')
xlabel('$x$ [-]')
ylabel('$\frac{dU}{dx}$','Rotation',0,'FontSize',20)
yline(0, 'Color', 'k','HandleVisibility','off')
ylim([-10 10])
plot(-mu,0,'o','MarkerSize',4,'LineWidth',4,'Color',[0.8500 0.3250 0.0980])
plot(1-mu,0,'o','MarkerSize',4,'LineWidth',4,'Color',[0.9290 0.6940 0.1250])
grid minor
legend('Earth', 'Moon')

% Set options for fzero, TolX is eps (default)
zero_options = optimset('Display', 'iter');

% Define zero-finding problem for zeros of dU(x,0,0)/dx, intervals defined 
% by graphical interpretation 
[xL1, fL1, flagL1] = fzero(fun, [0.6 1-mu-0.05],zero_options);
[xL2, fL2, flagL2] = fzero(fun, [1-mu+0.05 1.5],zero_options);
[xL3, fL3, flagL3] = fzero(fun, [-2 -mu-0.4],zero_options);

% --- Triangular points L4, L5 ---

% Use the analytical results of the triangular libration points to find
% their coordinates. The L4 and L5 points can be seen as two vertices of
% two equilateral triangles defined by the two attractors and the L4/L5 as
% vertices respectively.

% x-axis coordinate is the mean between the Earth-Moon location for both
% L4/L5
xL4 = ((1 - mu) - mu) / 2;
xL5 = xL4;

% y-axis coordinates are opposite, but same magnitude (equilateral
% triangles)
yL4 = sqrt(3)/2;
yL5 = -yL4;

% Define the L4 and L5 coordinates
L4 = [xL4, yL4];
L5 = [xL5, yL5];

% --- Post - processing ---

blue = [0 0.4470 0.7410];
red = [0.8500 0.3250 0.0980];
grey = [0.7 0.7 0.7];

figure
plot(-mu, 0, 'o', 'MarkerSize',4,'LineWidth',4,'Color',blue);
hold on
grid on
plot(1-mu, 0, 'o', 'MarkerSize',4,'LineWidth',4,'Color', grey);
plot(xL1, 0,'Marker', 'o', 'MarkerSize',4,'LineWidth',4,'Color',red);
plot(xL2, 0,'Marker', 'o', 'MarkerSize',4,'LineWidth',4,'Color',red);
plot(xL3, 0,'Marker', 'o', 'MarkerSize',4,'LineWidth',4,'Color',red);
plot(L4(1), L4(2),'Marker', 'o', 'MarkerSize',4,'LineWidth',4,'Color',red);
plot(L5(1), L5(2),'Marker', 'o', 'MarkerSize',4,'LineWidth',4,'Color',red);
fimplicit(@(x,y) (x+mu).^2 + y.^2 - 1, 'LineStyle', '-.', 'LineWidth', 0.3, 'Color','k');
fimplicit(@(x,y) (x-1+mu).^2 + y.^2 - 1, 'LineStyle','-.', 'LineWidth', 0.3, 'Color','k');
annotation("textarrow", [0.2 0.15], [0.4 0.51], 'String', 'L3')
annotation("textarrow", [0.5 0.6], [0.4 0.51], 'String', 'L1')
annotation("textarrow", [0.75 0.7], [0.4 0.51], 'String', 'L2')
annotation("textarrow", [0.3 0.5], [0.7 0.82], 'String', 'L4')
annotation("textarrow", [0.3 0.5], [0.36 0.24], 'String', 'L5')
axis equal
xlabel('x [-]');
ylabel('y [-]');
legend('Earth', 'Moon')

% Calculate Jacobi constant for L1 L2 L3 L4 L5

[~,J_L1] = dJ([xL1 0 0 0 0 0], mu);
[~,J_L2] = dJ([xL2 0 0 0 0 0], mu);
[~,J_L3] = dJ([xL3 0 0 0 0 0], mu);
[~,J_L4] = dJ([xL4 yL4 0 0 0 0], mu);
[~,J_L5] = dJ([xL5 yL5 0 0 0 0], mu);


%% Point 2: Find Halo orbit: given Jacobi constant C & initial guess

% Define acceptable Jacobi constant tolerance and x,z velocity components
% at (x,z) plane intersection
tol_r = 1e-10;
tol_v = 1e-10;

% Call the findHalo routine to find the Halo given: Jacobi constant C,
% initial guess for I.C., gravitational parameter mu, tolerances.
[x0_halo309, T_halo309, Jx0_halo309, xxt_f_309] = findHalo(xx0, C, mu, tol_v, tol_r, true);

% Propagate the found initial condition for nrev (C=3.09)
nrev = 1;
options  = odeset('reltol', 1e-12, 'abstol', 1e-12);
[~,xx309] = ode78(@(t,x) CRTBP_rhs(t,x,mu), [0 nrev*T_halo309], x0_halo309, options);

% Plot the propgated Halo (C=3.09) with I.C. point, L2 point and the Moon
figure;
plot3(xx309(:,1), xx309(:,2), xx309(:,3), 'LineWidth',2);
hold on;
plot3(xx309(1,1), xx309(1,2), xx309(1,3),'x', 'MarkerSize', 10, 'LineWidth',2);
plot3(xL2(1), 0 ,0, '+', 'MarkerSize', 10, 'LineWidth',2);
plot3(1-mu, 0 ,0, 'o', 'MarkerSize', 5,'Color',grey,'LineWidth',5);
grid on
axis square;
legend('Halo Orbit', 'Initial Condition', 'L2', 'Moon')
xlabel('$x \, [-]$'); ylabel('$y \, [-]$'); zlabel('$z \, [-]$');

% view(0,90)  % XY
% pause
% view(0,0)   % XZ
% pause
% view(90,0)  % YZ

%% Point 3: numerical continuation - family of Halo (C from 3.09 to 3.04)

% Physical data - 50 orbits with C 3.09 to 3.04 will be found
CC = linspace(3.09, 3.04, 10);

% Allocation of the interested variables
x0 = x0_halo309;
T  = T_halo309;
J = Jx0_halo309;
xxt_f_v = xxt_f_309;

% Define color map for graphical orbits representation
cmap = parula(length(CC));

figure;

plot3(xL2(1), 0 ,0,'o', 'MarkerSize', 5, 'LineWidth',5);
hold on;
plot3(1-mu, 0 ,0, 'o', 'MarkerSize', 5, 'LineWidth',5);
axis square;
grid on;

colorbar;
set(gca,'CLim',[3.040 3.09]);
xlabel('$x \, [-]$'); ylabel('$y \, [-]$'); zlabel('$z \, [-]$');
legend('L2', 'Moon','AutoUpdate','off');

plot3(xx309(:,1), xx309(:,2), xx309(:,3), 'LineWidth',2, 'Color', cmap(end,:));

for k = 2:length(CC)
    [x0_halo, T_halo, Jx0, xxt_f] = findHalo(x0(:,k-1), CC(k), mu, tol_v, tol_r, false);
    x0 = [x0, x0_halo];
    T = [T, T_halo];
    J = [J, Jx0];
    xxt_f_v = [xxt_f_v,xxt_f];
    [tt,xx] = ode78(@(t,x) CRTBP_rhs(t,x,mu), [0 T_halo], x0_halo, options);
    plot3(xx(:,1), xx(:,2), xx(:,3), 'LineWidth',2, 'Color', cmap(end-(k-1),:));
end

colorbar('Direction','reverse')

% view(0,90)  % XY
% pause
% view(0,0)   % XZ
% pause
% view(90,0)  % YZ


%% Functions 

function dU = dU(x, y, z, mu)
%------------------ Gradient of CRTBP potential ---------------------------
%
% This function calculates the gradient of the CRTBP potential expressed in
% the rotating frame.
% 
%--------------------------------------------------------------------------
%
% INPUT: 
% x          [1] x coordinate
% y          [1] y coordinate
% z          [1] z coordinate
% mu         [1] gravitational parameter
% 
% OUTPUT:
% dU         [3x1] gradient of U        
% 
% DEPENDENCIES:
% NONE
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------
    
% Allocate space
dU   = zeros(3,1);

% Calculate r1 and r2 distances (from E and M)
r1   = sqrt((x+mu)^2 + y^2 + z^2);
r2   = sqrt((x+mu-1)^2 + y^2 + z^2); 

% Calculate components 
dUdx = x - (mu+x)*(1-mu)/r1^3 + (1-mu-x)*mu/r2^3;
dUdy = y - y*(1-mu)/r1^3 - y*mu/r2^3;
dUdz = - z*(1-mu)/r1^3 - z*mu/r2^3;

% Build the gradient
dU = [dUdx; dUdy; dUdz];

end

function dUdx00 = dUdx00(x, mu)
%--------------x component of Gradient of CRTBP potential -----------------
%
% This function calculates the x component of the gradient on the x axis 
% of the CRTBP potential expressed inthe rotating frame.
% 
%--------------------------------------------------------------------------
%
% INPUT: 
% x          [n] x coordinate
% mu         [1] gravitational parameter
% 
% OUTPUT:
% dUdx00     [n] x component of gradient of U calculated along the x-axis      
% 
% DEPENDENCIES:
% NONE
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

%Partial derivative of potential wrt to x coordinate, evaluated at (x,0,0)
dUdx00 = x - (1-mu)*(x+mu)./abs(x+mu).^3 + mu*(1-x-mu)./abs(x+mu-1).^3;

end

function CRTBP_rhs = CRTBP_rhs(~,xx,mu)
% -------------------------------------------------------------------------
% This function allows to calculate the rhs of the CRTBP dynamics. All the 
% units are a-dimensional.
% INPUT:
% - xx           [6x1]: state vector
% - mu             [-]: adimensional mu parameter for the RTBP
%
% OUTPUT:
% - CRTBP_rhs    [6x1]: evaluated RHS of the state dynamics
%
% DEPENDENCIES:
% - dU(x,y,z,mu): gradient of the CRTBP potential
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Extract position and velocity components
x = xx(1); y = xx(2); z = xx(3); v_x = xx(4); v_y = xx(5); v_z = xx(6);

% Compute gradient of the CRTBP
du = dU(x,y,z,mu);
dUdx = du(1,1); dUdy = du(2,1); dUdz = du(3,1);

% Calculate rhs of the CRTBP
CRTBP_rhs = [v_x; v_y; v_z; 2*v_y + dUdx; -2*v_x + dUdy; dUdz];
    
end

function [CRTBP_STM, xx_end] = CRTBP_STM(tev,x0,mu)
% -------------------------------------------------------------------------
% This function allows to calculate the numerical STM of the CRTBP at time
% t_ev given a certain initial condition x0 and the state x propagated from 
% x0 to t_ev. The method to approximate the derivatives is the 
% forward-differences method.
% INPUT:
% - tev            [1]: time in which STM must be computed
% - x0           [6x1]: initial condition
% - mu             [-]: adimensional mu parameter for the CRTBP
%
% OUTPUT:
% - CRTBP_STM    [6x6]: approximated STM at t_ev with x0 as IC
% - xx_end       [6x1]: state 
% DEPENDENCIES:
% - CRTBP_rhs(t,x)
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Allocate space for the STM
CRTBP_STM = zeros(6,6);

% Propagate the state x0 from 0 to t_ev
options = odeset('AbsTol',1e-12, 'RelTol',1e-12);
[~,xx] = ode78(@(t,x) CRTBP_rhs(t,x,mu), [0 tev], x0, options);

% Extract solution at t_ev
xx_end = xx(end,:)';

% Loop over the six columns (each column is a perturbation along one
% component of the initial state)
for k = 1:6

    % Define the the perturbation vector
    eps_v = zeros(6,1);
    % Impose the perturbation to apply on the k-th component of x0
    eps_k = max(sqrt(eps), abs(xx_end(k))*sqrt(eps));
    % Define the k-th component of the pertubration vector
    eps_v(k) = eps_k; 

    % Propagate the initial condition perturbed along k-th component 
    [~,xxk] = ode78(@(t,x) CRTBP_rhs(t,x,mu), [0 tev], x0+eps_v, options);
    
    % calculate the k-th column of the STM with the forward scheme
    CRTBP_STM(:,k) = (xxk(end,:)' - xx_end) ./ eps_k;
end

end

function [value, isterminal, direction] = xz_plane_crossing(~, xx, isTerminal) 
% -------------------------------------------------------------------------
% Event function used in the differential correction scheme. This event
% function defines the if the solution is crossing of the xz plane.
% 
% INPUT:
% - xx           [6x1]: state to check
% - isTerminal  [bool]: flag which says if integration must stop
% 
% OUTPUT:
% - value          [1]: y component of the position from x
% - isterminal  [bool]: triggers when the xz crossing happens
% - direction   [bool]: direction to which is approached
% DEPENDENCIES:
% - NONE
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

value = xx(2);
isterminal = isTerminal; %1 --> stops, 0--> continue
direction = 0; % can be approached both ways

end

function [dJ,J] = dJ(xx, mu)
% -------------------------------------------------------------------------
% This function allows to calculate the analytical gradient of the Jacobi
% constant with respect to the orbital state. Optionally, it calculates
% also the Jacobi constant.
% INPUT:
% - xx           [6x1]: state
% - mu             [-]: adimensional mu parameter for the CRTBP
%
% OUTPUT:
% - dJ           [6x1]: gradient of the Jacobi constant
% - J              [1]: Jacobi constant (OPTIONAL)
% DEPENDENCIES:
% - NONE
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Extract the position and velocity components of the state
x = xx(1);
y = xx(2);
z = xx(3);
v_x = xx(4);
v_y = xx(5);
v_z = xx(6);

% Define distance of S/C from E and M
r1   = sqrt((x+mu)^2 + y^2 + z^2);
r2   = sqrt((x+mu-1)^2 + y^2 + z^2); 

% Velocity magnitude squared
v2   = v_x^2 + v_y^2 + v_z^2;

% CRTBP potential in rotating frame
OM = 0.5 * (x^2 + y^2) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu);

% Analytical gradient of the Jacobi constant
dJ = [2*x + ((2*mu + 2*x)*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (mu*(2*mu + 2*x - 2))/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
                          2*y - (2*mu*y)/((mu + x - 1)^2 + y^2 + z^2)^(3/2) + (2*y*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2);
                                (2*z*(mu - 1))/((mu + x)^2 + y^2 + z^2)^(3/2) - (2*mu*z)/((mu + x - 1)^2 + y^2 + z^2)^(3/2);
                                                                                                                     -2*v_x;
                                                                                                                     -2*v_y;
                                                                                                                     -2*v_z ];
% If requested, output also the Jacobi constant
if nargout > 1
    J = 2*OM - v2;
end

end

function [x0, T, Jx0, xxt_f] = findHalo(x0_guess, C, mu, tol_J, tol_v, dispIter)

% -------------------------------------------------------------------------
% This function allows to correct the initial guess on IC of an Halo orbit
% with presribed Jacobi constant. The initial guess must be not too far
% from the real condition and also it must be on the xz plane with
% perpendicular velocity. 
% INPUT:
% - x0_guess     [6x1]: initial condition first guess
% - C              [1]: Jacobi constant
% - mu             [1]: adimensional mu parameter for the CRTBP
% - tol_J          [1]: tolerance on the Jacobi constant
% - tol_v          [1]: tolerance on the x and z components of velocity
% - dispIter    [bool]: shows convergence results of the differential
%                       scheme
% OUTPUT:
% - x0           [6x1]: corrected initial condition of the Halo
% - T              [1]: period of the Halo
% - Jx0            [1]: Jacobi constant value at found solution
% - xxt_f        [6x1]: state at second crossing of xz plane
% DEPENDENCIES:
% - xz_plane_crossing(t,x,evtFlag): event function
% - CRTBP_rhs(t,x,mu): rhs of the CRTBP
% - CRTBP_STM(t,x,mu): numerical STM of the CRTBP at t
% - dJ(x,mu): jacobi constant gradient
%--------------------------------------------------------------------------
%
% AUTHOR: Daniele Paternoster
%
%--------------------------------------------------------------------------

% Function that finds the Halo Orbit associated to a certain mu, C jacobi
% constant and given intial guess on the xz plane.
% The method implemented is the STM differential correction 

% Event function flag set to 1 (stop integration when the xz event happens)
evtFlag = 1;

% Set solver options
options  = odeset('reltol', 1e-12, 'abstol', 1e-12,'Events',@(t,x) xz_plane_crossing(t,x,evtFlag));

% Propagate the initial guess until xz plane crossing
[tt,xx] = ode78(@(t,x) CRTBP_rhs(t,x,mu), [0 4], x0_guess, options);

% Define a first guess for the half period of the orbit.
t_fg = tt(end); %first guess

% Define state at final 
xxt_f = xx(end,:);

% Define variables that will be updated in the loop of the differential
% scheme
x0 = x0_guess;
t_f = t_fg;
k = 0;
[~,Jx0] = dJ(x0, mu);

% Set max iterations
Nmax = 30;

% Loop until all tolerances are below a certain threshold and the
% iterations are < Nmax
while (abs(xxt_f(2)) > 1e-10 || abs(xxt_f(4))>tol_v || abs(xxt_f(6))>tol_v || abs(C - Jx0)>tol_J) && k < Nmax
  
    % Propagate the new guess until new propagation time to find STM and
    % final state
    [STM, xxt_f]  = CRTBP_STM(t_f, x0, mu);
    xxt_f = xxt_f';
    % rhs of the CRTBP at final state
    f  = CRTBP_rhs(0,xxt_f, mu);

    % Gradient of Jacobi constant at updated initial condition
    [dJdx0, Jx0]  = dJ(x0,mu);

    % Define the augmented matrix of the system of eqs to solve
    augM(1:3,1:3) = STM([2 4 6], [1 3 5]);
    augM(4,4)     = 0;
    augM(4,1:3)   = dJdx0([1 3 5])';
    augM(1:3,4)   = f([2 4 6]);

    % define the RHS of the system of eqs to solve
    dx = [-xxt_f(2); -xxt_f(4); -xxt_f(6)];
    dX = [dx; C - Jx0];

    % Solve the augmented system to find new correction
    dX0 = augM \ dX;
    
    % Extract new corrections on initial guess (x,z,vy) components and
    % correction on final time for the propagation
    dx0 = dX0(1:3);
    dt = dX0(end);

    % Update the initial guess
    x0([1 3 5]) = x0([1 3 5]) + dx0;
    
    % Update the propagation time
    t_f         = t_f + dt;
   
    % Update the index of the iterations
    k = k+1;

end

% The total period is two times the propagation time (which was the one
% relative to the reachment of the xz plane).
T = 2*t_f;

% Display number of iterations 
if dispIter == true
    fprintf('The Halo orbit with Jacobi constant C = %f is found in %o iterations\n', C, k);
end

end 

function fontOptions(DefaultFontSize)
    set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
    set(groot, 'defaultLegendInterpreter', 'latex');
    set(groot, 'defaulttextinterpreter', 'latex');
    set(0, 'defaultAxesFontSize', DefaultFontSize)
end

