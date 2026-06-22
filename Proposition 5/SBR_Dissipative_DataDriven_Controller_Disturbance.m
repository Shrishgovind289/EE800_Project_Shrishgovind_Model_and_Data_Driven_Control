clc;
clear;
close all;

%% =========================================
% SELF-BALANCING ROBOT
% DATA-DRIVEN DISSIPATIVITY CONTROLLER
% USING INPUT-OUTPUT DATA MATRICES
% =========================================
%
% CATEGORY:
%   Data-Driven Dissipativity Controller
%
% Meaning:
%   1) Load data matrices: X0, X1, U0, Y0, Y1
%   2) Do NOT recover A_est, B_est
%   3) Do NOT use Ad, Bd for controller design
%   4) Design controller directly from data matrices
%
% Convention:
%   u(k) = K_data*x(k)
%
% Data-driven closed-loop map:
%   x(k+1) = Acl_dd*x(k)
%
% where:
%   Acl_dd = X1c*Q*(X0c*Q)^(-1)
%   K_data = U0c*Q*(X0c*Q)^(-1)

%% =========================================
% LOAD INPUT-OUTPUT DATA MATRICES
% =========================================

load('C:\Users\shris\Desktop\CPE 800 Project\Proposition 3\SBR_Dissipative_Data_Matrices.mat');

disp('================ LOADED DISSIPATIVE DATA MATRICES ================');

fprintf('size(X0) = [%d %d]\n', size(X0,1), size(X0,2));
fprintf('size(X1) = [%d %d]\n', size(X1,1), size(X1,2));
fprintf('size(U0) = [%d %d]\n', size(U0,1), size(U0,2));

if exist('Y0','var')
    fprintf('size(Y0) = [%d %d]\n', size(Y0,1), size(Y0,2));
end

if exist('Y1','var')
    fprintf('size(Y1) = [%d %d]\n', size(Y1,1), size(Y1,2));
end

%% =========================================
% BASIC DIMENSIONS
% =========================================

n = size(X0,1);      % number of states
m = size(U0,1);      % number of inputs
Tdata = size(X0,2);  % number of samples

% Sampling time fallback
if ~exist('Ts','var')
    Ts = 0.01;
end

%% =========================================
% DATA RANK CHECK
% =========================================

rank_data = rank([U0; X0]);
required_rank = n + m;

disp('================ DATA RANK CHECK =================');
fprintf('rank([U0; X0]) = %d\n', rank_data);
fprintf('required rank  = %d\n', required_rank);

if rank_data < required_rank
    error('Data is NOT rich enough for data-driven dissipativity control.');
else
    disp('Data is rich enough for data-driven dissipativity control.');
end

%% =========================================
% COMPRESS DATA TO WELL-CONDITIONED COLUMNS
% =========================================
%
% IMPORTANT:
%   Z has size (n+m) x Tdata = 5 x 1000.
%
%   We want to select useful COLUMNS from Z.
%   Therefore use:
%
%       qr(Z,'vector')
%
%   NOT:
%
%       qr(Z','vector')
%
%   Using qr(Z','vector') only gives 5 pivot columns and forces Tc = 5.
%   Correct version should give Tc = 20.

Z = [U0; X0];

[~,~,pivot_cols] = qr(Z,'vector');   % Correct: no transpose

Tc_desired = 20;

if length(pivot_cols) < Tc_desired
    Tc_desired = length(pivot_cols);
end

idx = pivot_cols(1:Tc_desired);

X0c = X0(:,idx);
X1c = X1(:,idx);
U0c = U0(:,idx);

Tc = size(X0c,2);

disp('================ COMPRESSED DATA =================');
fprintf('Original Tdata = %d\n', Tdata);
fprintf('Compressed Tc  = %d\n', Tc);
fprintf('rank([U0c; X0c]) = %d\n', rank([U0c; X0c]));

if rank([U0c; X0c]) < n+m
    error('Compressed data lost rank. Increase Tc_desired.');
end

%% =========================================
% DATA-DRIVEN DISSIPATIVITY CONTROLLER DESIGN
% =========================================

disp('====== DATA-DRIVEN DISSIPATIVITY CONTROLLER ======');

if ~exist('sdpvar','file')
    error('YALMIP not found. Install YALMIP and SDPT3/SeDuMi/MOSEK.');
end

% Decision variables
Q = sdpvar(Tc,n,'full');
P = sdpvar(n,n,'symmetric');
gamma = sdpvar(1,1);

% Data-driven closed-loop product
% This represents Acl*P directly from measured data
AclP = X1c * Q;

% Controller product
UP = U0c * Q;

% Disturbance input matrix
Ew = 0.02 * eye(n);
nw = size(Ew,2);

% Performance output weighting
% Higher weight on body angle phi
Cz = diag([0.1 0.02 1 0.1]);
nz = size(Cz,1);

% Data-driven dissipativity / bounded-real LMI
BRL = [ P,             AclP,          Ew,              zeros(n,nz);
        AclP',         P,             zeros(nw,n),     P*Cz';
        Ew',           zeros(nw,n),   gamma*eye(nw),   zeros(nw,nz);
        zeros(nz,n),   Cz*P,          zeros(nz,nw),    gamma*eye(nz) ];

constraints = [];

% Link Q to Lyapunov matrix
constraints = [constraints, X0c*Q == P];

% Positive definite P
constraints = [constraints, P >= 1e-5*eye(n)];

% Normalize P to avoid SDPT3 scaling issues
constraints = [constraints, trace(P) == 1];

% Gamma lower bound
constraints = [constraints, gamma >= 1e-4];

% Main dissipativity LMI
constraints = [constraints, BRL >= 1e-6*eye(size(BRL,1))];

% Keep controller product bounded
K_bound = 500;
constraints = [constraints, -K_bound <= UP <= K_bound];

% Objective
objective = gamma;

%% =========================================
% SOLVE DATA-DRIVEN DISSIPATIVITY LMI
% ACCEPT NEAR-FEASIBLE / NUMERICAL SDPT3 RESULT
% =========================================

ops = sdpsettings('solver','sdpt3', ...
                  'verbose',1, ...
                  'sdpt3.maxit',200, ...
                  'sdpt3.rmdepconstr',1, ...
                  'sdpt3.gaptol',1e-6, ...
                  'sdpt3.inftol',1e-4);

sol = optimize(constraints, objective, ops);

disp('================ SOLVER STATUS =================');
fprintf('YALMIP problem code = %d\n', sol.problem);
disp(sol.info);

solver_info_lower = lower(sol.info);

solver_ok = false;

% 0 = clean solve
if sol.problem == 0
    solver_ok = true;
end

% 4 = numerical problems, but often usable with SDPT3
if sol.problem == 4
    solver_ok = true;
end

% Accept common SDPT3 near-feasible exits
if contains(solver_info_lower, 'numerical problems')
    solver_ok = true;
end

if contains(solver_info_lower, 'lack of progress')
    solver_ok = true;
end

if contains(solver_info_lower, 'successfully solved')
    solver_ok = true;
end

if ~solver_ok
    error('Data-driven dissipativity LMI failed: %s', sol.info);
else
    disp('Solver result accepted. Extracting controller and checking stability...');
end

Q_val = value(Q);
P_val = value(P);
gamma_val = value(gamma);

% Safety check
if any(isnan(Q_val(:))) || any(isnan(P_val(:))) || isnan(gamma_val)
    error('Solver returned NaN values. Cannot continue.');
end

% Symmetrize P to remove tiny numerical asymmetry
P_val = 0.5*(P_val + P_val');

disp('Eigenvalues of P:');
disp(eig(P_val));

if min(eig(P_val)) <= 0
    warning('P is not strictly positive definite numerically. Continuing to eigenvalue check.');
end

% Data-driven controller gain
% u = Kx
K_data = U0c * Q_val / P_val;

% Data-driven closed-loop matrix
Acl_dd = X1c * Q_val / P_val;

disp('Data-driven dissipativity gain K_data:');
disp(K_data);

fprintf('Gamma value = %.6f\n', gamma_val);

%% =========================================
% CLOSED-LOOP STABILITY CHECK
% =========================================

disp('============ DATA-DRIVEN CLOSED-LOOP CHECK =========');

cl_eigs = eig(Acl_dd);

disp('Closed-loop eigenvalues from data:');
disp(cl_eigs);

max_eig_mag = max(abs(cl_eigs));
fprintf('Max eigenvalue magnitude = %.6f\n', max_eig_mag);

if max_eig_mag < 1
    disp('Data-driven dissipativity closed-loop system is stable.');
else
    error('Data-driven dissipativity controller is NOT stable. Do not use this result.');
end

%% =========================================
% DISTURBANCE RESPONSE TEST
% =========================================
%
% No A_est, B_est, Ad, or Bd are used here.
% Simulation uses the data-driven closed-loop map:
%
%   x(k+1) = Acl_dd*x(k)

Ntest = 2000;              % 20 seconds
u_sat_limit = 20;

disturbance_time = 10;     % seconds
dist_step = round(disturbance_time/Ts) + 1;

dist_angle_deg = -20;      % angle disturbance
dist_angvel_deg = -80;     % angular velocity disturbance

angle_threshold_deg = 1;   % recovery threshold

x_test = zeros(n, Ntest+1);
u_test = zeros(m, Ntest);
u_cmd_test = zeros(m, Ntest);

% Initial condition
x_test(:,1) = [0;
               0;
               deg2rad(45);
               0];

for k = 1:Ntest

    % u = Kx
    u_cmd = K_data * x_test(:,k);

    % Saturation for displayed control signal
    u_sat = max(min(u_cmd, u_sat_limit), -u_sat_limit);

    % Store
    u_cmd_test(:,k) = u_cmd;
    u_test(:,k) = u_sat;

    % Pure data-driven closed-loop update
    x_test(:,k+1) = Acl_dd * x_test(:,k);

    % Inject disturbance at 10 seconds
    if k == dist_step
        x_test(3,k+1) = x_test(3,k+1) + deg2rad(dist_angle_deg);
        x_test(4,k+1) = x_test(4,k+1) + deg2rad(dist_angvel_deg);

        fprintf('\n*** Disturbance injected at t = %.2f s ***\n', disturbance_time);
        fprintf('Added angle disturbance = %.2f deg\n', dist_angle_deg);
        fprintf('Added angular velocity disturbance = %.2f deg/s\n', dist_angvel_deg);
    end
end

%% =========================================
% RECOVERY TIME CALCULATION
% =========================================

angle_deg = rad2deg(x_test(3,:));

post_dist_idx = dist_step:length(angle_deg);
recovery_idx_rel = find(abs(angle_deg(post_dist_idx)) <= angle_threshold_deg, 1, 'first');

if isempty(recovery_idx_rel)
    recovery_time = NaN;
    disp('Recovery threshold was not reached within simulation window.');
else
    recovery_idx = post_dist_idx(recovery_idx_rel);
    recovery_time = (recovery_idx-1)*Ts - disturbance_time;
    fprintf('Recovery time after disturbance = %.4f s\n', recovery_time);
end

disp('================ DISTURBANCE TEST RESPONSE ===================');
fprintf('Initial angle (deg)                 = %.4f\n', rad2deg(x_test(3,1)));
fprintf('Disturbance angle added (deg)        = %.4f\n', dist_angle_deg);
fprintf('Disturbance angular velocity (deg/s) = %.4f\n', dist_angvel_deg);
fprintf('Final angle (deg)                   = %.6f\n', rad2deg(x_test(3,end)));
fprintf('Final state norm                    = %.6e\n', norm(x_test(:,end)));
fprintf('Max |u_cmd|                         = %.6f\n', max(abs(u_cmd_test(:))));
fprintf('Max |u_sat|                         = %.6f\n', max(abs(u_test(:))));

%% =========================================
% PLOTS
% =========================================

t = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

figure('Name','Data-Driven Dissipative Cart Position');
plot(t, x_test(1,:), 'LineWidth', 2);
xline(disturbance_time, '--r', 'Disturbance');
grid on;
xlabel('Time (s)');
ylabel('x (m)');
title('Data-Driven Dissipativity Controller: Cart Position');

figure('Name','Data-Driven Dissipative Linear Velocity');
plot(t, x_test(2,:), 'LineWidth', 2);
xline(disturbance_time, '--r', 'Disturbance');
grid on;
xlabel('Time (s)');
ylabel('xdot (m/s)');
title('Data-Driven Dissipativity Controller: Linear Velocity');

figure('Name','Data-Driven Dissipative Body Angle');
plot(t, rad2deg(x_test(3,:)), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--r', 'Disturbance');
yline(angle_threshold_deg, '--k', '+1 deg');
yline(-angle_threshold_deg, '--k', '-1 deg');
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
title('Data-Driven Dissipativity Controller: Body Angle with Disturbance');

figure('Name','Data-Driven Dissipative Angular Velocity');
plot(t, rad2deg(x_test(4,:)), 'LineWidth', 2);
xline(disturbance_time, '--r', 'Disturbance');
grid on;
xlabel('Time (s)');
ylabel('\phidot (deg/s)');
title('Data-Driven Dissipativity Controller: Angular Velocity');

figure('Name','Data-Driven Dissipative Control Input');
plot(t_input, u_test, 'LineWidth', 2);
xline(disturbance_time, '--r', 'Disturbance');
grid on;
xlabel('Time (s)');
ylabel('Control input u');
title('Data-Driven Dissipativity Controller: Control Input');

%% =========================================
% EXPORT ANIMATION VIDEO ONLY
% =========================================

disp('========== EXPORTING DATA-DRIVEN DISSIPATIVITY VIDEO ==========');

video_name = 'data_driven_dissipativity_robot_animation.mp4';

v = VideoWriter(video_name, 'MPEG-4');
v.FrameRate = 20;
open(v);

fig_anim = figure('Name','Data-Driven Dissipativity Robot Animation Video', ...
                  'NumberTitle','off', ...
                  'WindowStyle','normal');

set(fig_anim, 'Position', [100 100 960 540]);

wheel_radius = 0.05;
body_length  = 0.5;

for k = 1:10:length(x_test)

    clf(fig_anim);
    hold on;
    grid on;
    axis equal;

    % Fixed view
    xlim([-1.5 1.5]);
    ylim([-0.1 1.0]);

    % States
    x_pos = x_test(1,k);
    phi   = x_test(3,k);

    % Geometry
    wheel_x = x_pos;
    wheel_y = wheel_radius;

    body_x = wheel_x + body_length*sin(phi);
    body_y = wheel_y + body_length*cos(phi);

    % Ground
    plot([-5 5], [0 0], 'k', 'LineWidth', 2);

    % Wheel
    theta = linspace(0, 2*pi, 80);
    plot(wheel_x + wheel_radius*cos(theta), ...
         wheel_y + wheel_radius*sin(theta), ...
         'b', 'LineWidth', 2);

    % Body
    plot([wheel_x body_x], [wheel_y body_y], ...
         'r', 'LineWidth', 4);

    % Center of mass
    plot(body_x, body_y, 'ko', 'MarkerFaceColor', 'k');

    % Disturbance marker
    if (k-1)*Ts >= disturbance_time
        text(-1.35, 0.9, 'Disturbance applied', ...
             'Color', 'r', 'FontWeight', 'bold');
    end

    title(sprintf('Data-Driven Dissipativity Animation | Time = %.2f s | Angle = %.2f deg', ...
          (k-1)*Ts, rad2deg(phi)));

    xlabel('Position x (m)');
    ylabel('Height (m)');

    drawnow;

    frame = getframe(fig_anim);
    writeVideo(v, frame);
end

close(v);

disp('Video saved successfully:');
disp(fullfile(pwd, video_name));

%% =========================================
% PLAY SAVED VIDEO
% =========================================

disp('Playing saved video...');

video_path = fullfile(pwd, video_name);

if exist(video_path, 'file')
    try
        implay(video_path);
    catch
        open(video_path);
    end
else
    warning('Video file was not found.');
end

%% =========================================
% SAVE RESULTS
% =========================================

save('SBR_DataDriven_Dissipativity_Controller_Result.mat', ...
    'X0','X1','U0','Y0','Y1', ...
    'X0c','X1c','U0c', ...
    'Q_val','P_val', ...
    'K_data', ...
    'Acl_dd','cl_eigs', ...
    'x_test','u_test','u_cmd_test', ...
    'disturbance_time','dist_step', ...
    'dist_angle_deg','dist_angvel_deg', ...
    'angle_threshold_deg','recovery_time', ...
    'Ts','Ew','Cz','gamma_val','u_sat_limit', ...
    'rank_data','required_rank','video_name');

disp('================ FINISHED ========================');
disp('Saved results to SBR_DataDriven_Dissipativity_Controller_Result.mat');
disp('Data-driven dissipativity controller completed.');