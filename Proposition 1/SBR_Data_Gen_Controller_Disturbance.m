clc;
clear;
close all;

%% =========================================
% SELF-BALANCING ROBOT
% PURE MODEL-BASED LQR CONTROLLER
% + DATA GENERATION FOR DATA-DRIVEN CONTROLLER
% + DISTURBANCE TEST
% =========================================
%
% IMPORTANT CLASSIFICATION:
%   This script is PURE MODEL-BASED for controller design.
%
%   Controller design uses ONLY:
%       Ad, Bd
%
%   Controller design does NOT use:
%       X0, X1, U0
%       A_est, B_est
%       pinv([U0; X0])
%
%   Data generation is included ONLY to create:
%       SBR_Data_Matrices.mat
%   for the separate data-driven controller script.
%
% State vector:
%   x = [ cart position;
%         cart velocity;
%         body angle phi;
%         body angular velocity phidot ]
%
% Control convention:
%   u(k) = -K_lqr*x(k)
%   x(k+1) = Ad*x(k) + Bd*u(k)
%          = (Ad - Bd*K_lqr)*x(k)

%% =========================================
% PHYSICAL PARAMETERS
% =========================================
M = 0.5;         % chassis mass (kg)
m_body = 0.2;    % body mass (kg)
l = 0.3;         % COM distance from wheel axle (m)
I = 0.006;       % body inertia (kg*m^2)
b = 0.1;         % viscous friction coefficient
g = 9.81;        % gravity (m/s^2)

%% =========================================
% MODEL CONSTANTS
% =========================================
a = M + m_body;
c = m_body * l;
d = I + m_body * l^2;
Delta = a*d + c^2;

%% =========================================
% CONTINUOUS-TIME LINEARIZED MODEL
% x = [x; xdot; phi; phidot]
% =========================================
A = [0, 1, 0, 0;
     0, -(d*b)/Delta, -(c*m_body*g*l)/Delta, 0;
     0, 0, 0, 1;
     0, -(c*b)/Delta, (a*m_body*g*l)/Delta, 0];

B = [0;
     d/Delta;
     0;
     c/Delta];

C = eye(4);
D = zeros(4,1);

fprintf('\n================ PURE MODEL-BASED LQR + DATA GENERATION ================\n');
fprintf('Controller is designed directly from Ad and Bd.\n');
fprintf('Data matrices are generated only for the separate data-driven controller.\n');

fprintf('\n================ CONTINUOUS MODEL ================\n');
disp('A matrix:');
disp(A);
disp('B matrix:');
disp(B);

%% =========================================
% DISCRETIZATION
% =========================================
Ts = 0.01;       % sampling time (s)

n = size(A,1);
m = size(B,2);
p = size(C,1);

% Exact zero-order-hold discretization using matrix exponential
Maug = [A B;
        zeros(m, n+m)];

Md = expm(Maug * Ts);

Ad = Md(1:n,1:n);
Bd = Md(1:n,n+1:n+m);
Cd = C;
Dd = D;

fprintf('\n================ DISCRETE MODEL ==================\n');
disp('Ad matrix:');
disp(Ad);
disp('Bd matrix:');
disp(Bd);

%% =========================================
% PURE MODEL-BASED DISCRETE LQR DESIGN
% =========================================
% Larger Q weight on phi makes the robot prioritize upright balance.
Q_lqr = diag([10 1 300 1]);
R_lqr = 1;

fprintf('\n================ PURE MODEL-BASED LQR DESIGN ===============\n');

if exist('dlqr','file') == 2
    disp('Using dlqr from Control System Toolbox...');
    K_lqr = dlqr(Ad, Bd, Q_lqr, R_lqr);
else
    disp('dlqr not found. Using iterative discrete Riccati solution...');

    P = Q_lqr;
    max_iter = 5000;
    tol = 1e-10;

    for iter = 1:max_iter
        P_next = Ad' * P * Ad ...
            - Ad' * P * Bd * ((R_lqr + Bd' * P * Bd) \ (Bd' * P * Ad)) ...
            + Q_lqr;

        if norm(P_next - P, 'fro') < tol
            P = P_next;
            fprintf('Riccati iteration converged in %d iterations.\n', iter);
            break;
        end

        P = P_next;
    end

    if iter == max_iter
        warning('Riccati iteration did not fully converge.');
    end

    K_lqr = (R_lqr + Bd' * P * Bd) \ (Bd' * P * Ad);
end

fprintf('\nPure model-based LQR gain K_lqr:\n');
disp(K_lqr);

%% =========================================
% CLOSED-LOOP STABILITY CHECK
% =========================================
Acl = Ad - Bd*K_lqr;
cl_eigs = eig(Acl);

fprintf('\n============ CLOSED-LOOP STABILITY CHECK =========\n');
disp('Acl = Ad - Bd*K_lqr:');
disp(Acl);

disp('Closed-loop eigenvalues:');
disp(cl_eigs);

if all(abs(cl_eigs) < 1)
    disp('Closed-loop discrete-time system is stable.');
else
    warning('Closed-loop discrete-time system is NOT stable.');
end

%% =========================================
% DATA GENERATION FOR DATA-DRIVEN CONTROLLER
% =========================================
% This section creates the data matrices used by:
%   SBR_DataDriven_Controller_Disturbance.m
%
% It does NOT affect K_lqr.
% It can be deleted and K_lqr will remain unchanged.

fprintf('\n================ DATA GENERATION FOR DATA-DRIVEN CONTROL ================\n');

Tdata = 1000;   % number of data samples

x0_data = [0;
           0;
           deg2rad(45);
           0];

% Temporary stabilizing feedback used only to collect safe, bounded data.
% This is NOT the final controller.
Ktemp = [-1.0  -2.0  18.0  3.5];
Acl_temp = Ad - Bd*Ktemp;

fprintf('Temporary data-collection closed-loop eigenvalues:\n');
disp(eig(Acl_temp));

% Persistently exciting input signal
u_amp = 0.15;
hold_steps = 8;

num_blocks = ceil(Tdata/hold_steps);
u_blocks = u_amp * sign(randn(1,num_blocks));
u_exc = repelem(u_blocks, hold_steps);
u_exc = u_exc(1:Tdata);

% Storage
x_data = zeros(n, Tdata+1);
y_data = zeros(p, Tdata);
u_data = zeros(m, Tdata);

x_data(:,1) = x0_data;

% Simulate data collection using the true model Ad, Bd
for k = 1:Tdata
    u_data(:,k) = -Ktemp*x_data(:,k) + u_exc(k);
    y_data(:,k) = Cd*x_data(:,k) + Dd*u_data(:,k);
    x_data(:,k+1) = Ad*x_data(:,k) + Bd*u_data(:,k);
end

% Build data matrices
X0 = x_data(:,1:Tdata);
X1 = x_data(:,2:Tdata+1);
U0 = u_data(:,1:Tdata);
Y0 = y_data(:,1:Tdata);

Y1 = zeros(p,Tdata);
Y1(:,1:Tdata-1) = y_data(:,2:Tdata);
Y1(:,Tdata) = Cd*x_data(:,Tdata+1) + Dd*u_data(:,Tdata);

fprintf('\n================ DATA MATRIX SIZES ===============\n');
fprintf('size(X0) = [%d %d]\n', size(X0,1), size(X0,2));
fprintf('size(X1) = [%d %d]\n', size(X1,1), size(X1,2));
fprintf('size(U0) = [%d %d]\n', size(U0,1), size(U0,2));
fprintf('size(Y0) = [%d %d]\n', size(Y0,1), size(Y0,2));
fprintf('size(Y1) = [%d %d]\n', size(Y1,1), size(Y1,2));

% Rank check for data-driven control
rank_data = rank([U0; X0]);
required_rank = n + m;

fprintf('\n================ DATA RANK CHECK =================\n');
fprintf('rank([U0; X0]) = %d\n', rank_data);
fprintf('required rank  = %d\n', required_rank);

if rank_data < required_rank
    warning('Data is NOT rich enough for pure data-driven control. Regenerate or increase excitation.');
else
    disp('Data is rich enough for pure data-driven control.');
end

% Save data matrices needed by the data-driven controller.
% Ad/Bd are saved only for optional comparison/debugging.
% The pure data-driven controller should NOT use Ad/Bd for design.
save('SBR_Data_Matrices.mat', ...
    'X0','X1','U0','Y0','Y1', ...
    'Ts','Tdata','n','m','p', ...
    'Ad','Bd','Cd','Dd', ...
    'A','B','C','D', ...
    'Ktemp','u_amp','hold_steps', ...
    'M','m_body','l','I','b','g','a','c','d','Delta');

fprintf('\nSaved data matrices to SBR_Data_Matrices.mat\n');

%% =========================================
% OPTIONAL MODEL RECOVERY CHECK ONLY
% =========================================
% This section is only a diagnostic to verify the generated data.
% It is NOT used for K_lqr design.

AB_check = X1 * pinv([U0; X0]);
B_check = AB_check(:,1:m);
A_check = AB_check(:,m+1:m+n);

fprintf('\n================ OPTIONAL DATA QUALITY CHECK ================\n');
fprintf('||A_check - Ad|| = %.6e\n', norm(A_check - Ad));
fprintf('||B_check - Bd|| = %.6e\n', norm(B_check - Bd));
fprintf('NOTE: A_check/B_check are NOT used for controller design.\n');

%% =========================================
% ROBUSTNESS TEST WITH PARAMETER VARIATION
% =========================================
% This does not affect the controller design.
% It only checks how the fixed model-based K_lqr performs
% if the real plant has different friction/inertia.

fprintf('\n================ ROBUSTNESS TEST =================\n');

b_values = [0.9*b, 1.1*b];
I_values = [0.9*I, 1.1*I];

for i = 1:length(b_values)
    for j = 1:length(I_values)

        b_test = b_values(i);
        I_test = I_values(j);

        a_test = M + m_body;
        c_test = m_body * l;
        d_test = I_test + m_body * l^2;
        Delta_test = a_test*d_test + c_test^2;

        A_test = [0, 1, 0, 0;
                  0, -(d_test*b_test)/Delta_test, -(c_test*m_body*g*l)/Delta_test, 0;
                  0, 0, 0, 1;
                  0, -(c_test*b_test)/Delta_test, (a_test*m_body*g*l)/Delta_test, 0];

        B_test = [0;
                  d_test/Delta_test;
                  0;
                  c_test/Delta_test];

        Maug_test = [A_test B_test;
                     zeros(m, n+m)];

        Md_test = expm(Maug_test * Ts);

        Ad_test = Md_test(1:n,1:n);
        Bd_test = Md_test(1:n,n+1:n+m);

        Acl_test = Ad_test - Bd_test*K_lqr;
        eig_test = eig(Acl_test);

        fprintf('\n--- Test Case ---\n');
        fprintf('b = %.4f, I = %.6f\n', b_test, I_test);
        fprintf('Stable? %d\n', all(abs(eig_test) < 1));
        fprintf('Max eigenvalue magnitude = %.4f\n', max(abs(eig_test)));

        % Quick linear simulation on varied model
        x_sim = zeros(n, 301);
        x_sim(:,1) = [0; 0; deg2rad(45); 0];

        for k = 1:300
            u_sim = -K_lqr*x_sim(:,k);
            x_sim(:,k+1) = Ad_test*x_sim(:,k) + Bd_test*u_sim;
        end

        fprintf('Final angle after 3 s (deg) = %.4f\n', rad2deg(x_sim(3,end)));
    end
end

%% =========================================
% CLOSED-LOOP TEST SIMULATION WITH DISTURBANCE
% =========================================
Ntest = 2000;          % 20 seconds total simulation
u_sat_limit = 20;      % actuator saturation limit

% Disturbance settings
% At 10 seconds, apply:
%   phi    <- phi - 20 deg
%   phidot <- phidot - 80 deg/s
disturbance_time = 10;
dist_step = round(disturbance_time/Ts);
dist_angle_deg = -20;
dist_angvel_deg = -80;

dist_vec = zeros(n,1);
dist_vec(3) = deg2rad(dist_angle_deg);
dist_vec(4) = deg2rad(dist_angvel_deg);
disturbance_applied = false;

x_test = zeros(n, Ntest+1);
u_test = zeros(m, Ntest);
u_cmd_test = zeros(m, Ntest);

% Initial condition: 45 deg tilt
x_test(:,1) = [0;
               0;
               deg2rad(45);
               0];

for k = 1:Ntest

    % Model-based LQR law
    u_cmd = -K_lqr*x_test(:,k);

    % Saturate the actuator command for realistic testing
    u_sat = max(min(u_cmd, u_sat_limit), -u_sat_limit);

    % Store input history
    u_cmd_test(:,k) = u_cmd;
    u_test(:,k) = u_sat;

    % Pure model-based plant update using Ad and Bd
    x_test(:,k+1) = Ad*x_test(:,k) + Bd*u_sat;

    % Apply impulse-like state disturbance at 10 seconds
    if k == dist_step
        x_test(:,k+1) = x_test(:,k+1) + dist_vec;
        disturbance_applied = true;
        fprintf('\n*** Disturbance applied at t = %.2f s: phi %+g deg, phidot %+g deg/s ***\n', ...
            k*Ts, dist_angle_deg, dist_angvel_deg);
    end
end

fprintf('\n================ TEST RESPONSE ===================\n');
fprintf('Initial angle (deg) = %.4f\n', rad2deg(x_test(3,1)));
fprintf('Final angle (deg)   = %.6f\n', rad2deg(x_test(3,end)));
fprintf('Final state norm    = %.6e\n', norm(x_test(:,end)));
fprintf('Max |u_cmd|         = %.6f\n', max(abs(u_cmd_test(:))));
fprintf('Max |u_sat|         = %.6f\n', max(abs(u_test(:))));
fprintf('Disturbance applied = %d\n', disturbance_applied);

% Recovery time after disturbance: first time |phi| <= 1 deg after disturbance
angle_threshold_deg = 1;
idx_after_dist = (dist_step+1):length(x_test);
recover_idx_local = find(abs(rad2deg(x_test(3,idx_after_dist))) <= angle_threshold_deg, 1, 'first');

if ~isempty(recover_idx_local)
    recover_idx = idx_after_dist(recover_idx_local);
    recovery_time = (recover_idx-1)*Ts - disturbance_time;
    fprintf('Recovery time to |phi| <= %.1f deg = %.4f s\n', angle_threshold_deg, recovery_time);
else
    recovery_time = NaN;
    fprintf('Recovery time to |phi| <= %.1f deg = NOT reached\n', angle_threshold_deg);
end

%% =========================================
% SAVE MODEL-BASED RESULTS
% =========================================
save('SBR_ModelBased_LQR_With_DataGeneration_Disturbance_Result.mat', ...
    'A','B','C','D', ...
    'Ad','Bd','Cd','Dd', ...
    'K_lqr','Q_lqr','R_lqr', ...
    'Acl','cl_eigs', ...
    'x_test','u_test','u_cmd_test', ...
    'disturbance_time','dist_step','dist_angle_deg','dist_angvel_deg', ...
    'dist_vec','recovery_time', ...
    'Ts','M','m_body','l','I','b','g','a','c','d','Delta');

fprintf('\n================ RESULTS SAVED ========================\n');
disp('Saved results to SBR_ModelBased_LQR_With_DataGeneration_Disturbance_Result.mat');

%% =========================================
% PLOTS: ALL STATES + CONTROL INPUT
% =========================================
t_state = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

% 1) Robot / Cart Position
figure;
plot(t_state, x_test(1,:), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('x (m)');
title('Pure Model-Based LQR: Robot / Cart Position Response');

% 2) Linear Velocity
figure;
plot(t_state, x_test(2,:), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('xdot (m/s)');
title('Pure Model-Based LQR: Linear Velocity Response');

% 3) Body Angle
figure;
plot(t_state, rad2deg(x_test(3,:)), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
yline(angle_threshold_deg, ':', '+1 deg');
yline(-angle_threshold_deg, ':', '-1 deg');
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
title('Pure Model-Based LQR: Body Angle Response');

% 4) Angular Velocity
figure;
plot(t_state, rad2deg(x_test(4,:)), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('\phidot (deg/s)');
title('Pure Model-Based LQR: Angular Velocity Response');

% 5) Control Input
figure;
plot(t_input, u_cmd_test, '--', 'LineWidth', 1.2);
hold on;
plot(t_input, u_test, 'LineWidth', 2);
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Control Input');
legend('u command','u saturated','Location','best');
title('Pure Model-Based LQR: Control Input Response');

%% =========================================
% SAFE VIDEO GENERATION
% =========================================
opengl software;

video_name = 'pure_model_based_lqr_with_data_generation_disturbance.mp4';

v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = 10;
open(v);

fig = figure('Visible','on', 'Name','Pure Model-Based LQR Animation');

wheel_radius = 0.05;
body_length = 0.5;
frame_count = 0;

for k = 1:5:length(x_test)

    clf(fig);
    hold on;
    grid on;
    axis equal;
    xlim([-1 1]);
    ylim([-0.1 1]);

    x_pos = x_test(1,k);
    phi = x_test(3,k);

    wheel_x = x_pos;
    wheel_y = wheel_radius;

    body_x = wheel_x + body_length*sin(phi);
    body_y = wheel_y + body_length*cos(phi);

    % Ground
    plot([-5 5],[0 0],'k','LineWidth',2);

    % Wheel
    theta = linspace(0,2*pi,50);
    plot(wheel_x + wheel_radius*cos(theta), ...
         wheel_y + wheel_radius*sin(theta), ...
         'b','LineWidth',2);

    % Body
    plot([wheel_x body_x], [wheel_y body_y], ...
         'r','LineWidth',4);

    % COM
    plot(body_x, body_y, 'ko','MarkerFaceColor','k');

    title(sprintf('Pure Model-Based LQR + Disturbance | t = %.2f s | Angle = %.2f deg', ...
          (k-1)*Ts, rad2deg(phi)));

    xlabel('Position x (m)');
    ylabel('Height (m)');

    drawnow;

    frame = getframe(fig);

    if ~isempty(frame.cdata)
        writeVideo(v, frame);
        frame_count = frame_count + 1;
    end
end

close(v);

fprintf('Frames written: %d\n', frame_count);

if frame_count > 0
    disp('Video successfully saved:');
    disp(fullfile(pwd, video_name));
else
    warning('No frames captured — video not created.');
end

fprintf('\n================ SCRIPT COMPLETE ========================\n');