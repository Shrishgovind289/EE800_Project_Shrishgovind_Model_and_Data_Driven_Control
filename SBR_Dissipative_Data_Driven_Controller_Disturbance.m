clc;
clear;
close all;

%% =========================================
% PURE DATA-DRIVEN CONTROLLER USING SDPT3
%
% Control convention:
%   u(k) = K_data*x(k)
%
% No A_est, B_est, Ad, Bd, A, or B are used.
%
% Includes:
%   1) Pure data-driven LMI controller
%   2) State disturbance at 10 seconds
%   3) Recovery time calculation
%   4) Plots
%   5) Animation video
% =========================================

%% Load data
load('SBR_Dissipative_Data_Matrices.mat');

disp('================ LOADED DATA =====================');
fprintf('size(X0) = [%d %d]\n', size(X0,1), size(X0,2));
fprintf('size(X1) = [%d %d]\n', size(X1,1), size(X1,2));
fprintf('size(U0) = [%d %d]\n', size(U0,1), size(U0,2));

%% Dimensions
n = size(X0,1);
T = size(X0,2);
m = size(U0,1);

%% Rank check
rank_data = rank([U0; X0]);
required_rank = n + m;

disp('================ RANK CHECK ======================');
fprintf('rank([U0; X0]) = %d\n', rank_data);
fprintf('required rank  = %d\n', required_rank);

if rank_data < required_rank
    error('Data is not rich enough for data-driven control.');
else
    disp('Data is rich enough for pure data-driven control.');
end

%% Check YALMIP
if ~exist('sdpvar','file')
    error('YALMIP not found.');
end

%% =========================================
% COMPRESS DATA TO MINIMUM INDEPENDENT COLUMNS
% =========================================
Z = [U0; X0];

[~,~,pivot_cols] = qr(Z','vector');
idx = pivot_cols(1:n+m);

X0c = X0(:,idx);
X1c = X1(:,idx);
U0c = U0(:,idx);

Tc = size(X0c,2);

disp('================ COMPRESSED DATA =================');
fprintf('Original T = %d\n', T);
fprintf('Compressed T = %d\n', Tc);
fprintf('rank([U0c; X0c]) = %d\n', rank([U0c; X0c]));

%% =========================================
% DATA-DRIVEN STABILIZATION LMI
% convention: u = Kx
% =========================================
Q = sdpvar(Tc,n,'full');

eps_lmi = 1e-9;

% Pole-speed constraint
% 0.99 is safer. 0.98 may fail depending on the data.
alpha = 0.99;

P_data = X0c * Q;
X1Q    = X1c * Q;

P_sym = 0.5*(P_data + P_data');

constraints = [];
constraints = [constraints, P_data - P_data' == 0];
constraints = [constraints, P_sym >= eps_lmi*eye(n)];

constraints = [constraints, ...
    [alpha*P_sym, X1Q;
     X1Q',        alpha*P_sym] >= eps_lmi*eye(2*n)];

objective = [];

ops = sdpsettings('solver','sdpt3', ...
                  'verbose',1, ...
                  'sdpt3.maxit',300);

sol = optimize(constraints, objective, ops);

if sol.problem ~= 0
    error('Data-driven LMI failed using SDPT3: %s', sol.info);
end

%% Extract controller
Q_val = value(Q);
P_val = value(P_sym);

K_data = U0c * Q_val / P_val;
Acl_dd = X1c * Q_val / P_val;

disp('================ DATA-DRIVEN GAIN ================');
disp('K_data for u = K_data*x:');
disp(K_data);

%% Closed-loop check
cl_eigs = eig(Acl_dd);

disp('============ DATA-DRIVEN CLOSED-LOOP CHECK =======');
disp('Closed-loop eigenvalues from data:');
disp(cl_eigs);

fprintf('Max eigenvalue magnitude = %.6f\n', max(abs(cl_eigs)));

if all(abs(cl_eigs) < 1)
    disp('Data-driven closed-loop system is stable.');
else
    disp('Data-driven closed-loop system is NOT stable.');
end

%% =========================================
% CLOSED-LOOP SIMULATION WITH DISTURBANCE
% x(k+1) = Acl_dd*x(k)
% =========================================
Ntest = 2000;     % 20 seconds if Ts = 0.01

x_test = zeros(n,Ntest+1);
u_test = zeros(m,Ntest);

x_test(:,1) = [0;
               0;
               deg2rad(45);
               0];

%% Disturbance settings
disturbance_time = 10;                 % seconds
disturbance_step = round(disturbance_time/Ts);

dist_angle_deg = -20;       % sudden body angle disturbance
dist_angvel_deg = -80;      % sudden angular velocity disturbance deg/s

angle_disturbance = deg2rad(dist_angle_deg);       % sudden 10 degree tilt
velocity_disturbance = deg2rad(dist_angvel_deg);    % angular velocity kick

for k = 1:Ntest

    % Apply disturbance at 10 seconds
    if k == disturbance_step
        x_test(3,k) = x_test(3,k) + angle_disturbance;
        x_test(4,k) = x_test(4,k) + velocity_disturbance;

        disp('================ DISTURBANCE APPLIED =============');
        fprintf('Disturbance time = %.2f s\n', disturbance_time);
        fprintf('Angle disturbance = %.2f deg\n', rad2deg(angle_disturbance));
        fprintf('Angular velocity disturbance = %.2f deg/s\n', rad2deg(velocity_disturbance));
    end

    u_test(:,k) = K_data * x_test(:,k);
    x_test(:,k+1) = Acl_dd * x_test(:,k);
end

t = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

%% =========================================
% RECOVERY TIME CALCULATION
% =========================================
angle_threshold_deg = 1;
angle_threshold = deg2rad(angle_threshold_deg);

post_disturbance_angle = abs(x_test(3,disturbance_step:end));

recovery_index = find(post_disturbance_angle <= angle_threshold, 1, 'first');

if isempty(recovery_index)
    recovery_time = NaN;
    disp('Robot did not recover within simulation time.');
else
    recovery_time = (recovery_index - 1)*Ts;
    fprintf('Recovery time after disturbance = %.4f seconds\n', recovery_time);
end

%% Test response output
disp('================ TEST RESPONSE ===================');
fprintf('Initial angle (deg) = %.4f\n', rad2deg(x_test(3,1)));
fprintf('Final angle (deg)   = %.6f\n', rad2deg(x_test(3,end)));
fprintf('Final state norm    = %.6e\n', norm(x_test(:,end)));
fprintf('Max |u|             = %.6f\n', max(abs(u_test(:))));

%% =========================================
% PLOTS
% =========================================
figure('Name','Data-Driven Cart Position');
plot(t, x_test(1,:), 'LineWidth', 2);
xline(disturbance_time, '--r', 'Disturbance');
grid on;
xlabel('Time (s)');
ylabel('x (m)');
title('Data-Driven Controller: Cart Position');

figure('Name','Data-Driven Linear Velocity');
plot(t, x_test(2,:), 'LineWidth', 2);
xline(disturbance_time, '--r', 'Disturbance');
grid on;
xlabel('Time (s)');
ylabel('xdot (m/s)');
title('Data-Driven Controller: Linear Velocity');

figure('Name','Data-Driven Body Angle');
plot(t, rad2deg(x_test(3,:)), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--r', 'Disturbance');
yline(angle_threshold_deg, '--k', '+1 deg');
yline(-angle_threshold_deg, '--k', '-1 deg');
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
title('Data-Driven Controller: Body Angle with Disturbance');

figure('Name','Data-Driven Angular Velocity');
plot(t, rad2deg(x_test(4,:)), 'LineWidth', 2);
xline(disturbance_time, '--r', 'Disturbance');
grid on;
xlabel('Time (s)');
ylabel('\phidot (deg/s)');
title('Data-Driven Controller: Angular Velocity');

figure('Name','Data-Driven Control Input');
plot(t_input, u_test, 'LineWidth', 2);
xline(disturbance_time, '--r', 'Disturbance');
grid on;
xlabel('Time (s)');
ylabel('Control input u');
title('Data-Driven Controller: Control Input');

%% ========================================
% Save results
% =========================================
save('SBR_Dissipative_DataDriven_Controller_Result.mat', ...
    'K_data', ...
    'Q_val', 'P_val', ...
    'Acl_dd', 'cl_eigs', ...
    'x_test', 'u_test', ...
    'Ts', 'alpha', ...
    'disturbance_time', ...
    'angle_disturbance', ...
    'velocity_disturbance', ...
    'recovery_time');

disp('================ RESULTS SAVED ===================');
disp('Saved result as SBR_Dissipative_DataDriven_Controller_Result.mat');

%% =========================================
% SHOW + SAVE DATA-DRIVEN ROBOT ANIMATION
% =========================================
disp('========== GENERATING AND SAVING ANIMATION ==========');

opengl software

video_name = 'dissipative_data_driven_robot_disturbance.mp4';

v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = 10;
open(v);

fig = figure('Name','Data-Driven Robot Animation with Disturbance', ...
             'Visible','on');

wheel_radius = 0.05;
body_length = 0.5;

frame_count = 0;

for k = 1:10:length(x_test)

    clf(fig);
    hold on;
    grid on;
    axis equal;

    xlim([-1 1]);
    ylim([-0.1 1]);

    x_pos = x_test(1,k);
    phi   = x_test(3,k);

    wheel_x = x_pos;
    wheel_y = wheel_radius;

    body_x = wheel_x + body_length*sin(phi);
    body_y = wheel_y + body_length*cos(phi);

    % Ground
    plot([-5 5],[0 0],'k','LineWidth',2);

    % Wheel
    theta = linspace(0,2*pi,60);
    plot(wheel_x + wheel_radius*cos(theta), ...
         wheel_y + wheel_radius*sin(theta), ...
         'b','LineWidth',2);

    % Body
    plot([wheel_x body_x], [wheel_y body_y], ...
         'r','LineWidth',4);

    % COM
    plot(body_x, body_y, 'ko','MarkerFaceColor','k');

    current_time = (k-1)*Ts;

    if abs(current_time - disturbance_time) < 0.05
        title_text = sprintf('DISTURBANCE APPLIED | Time = %.2f s | Angle = %.2f deg', ...
                             current_time, rad2deg(phi));
    else
        title_text = sprintf('Time = %.2f s | Angle = %.2f deg', ...
                             current_time, rad2deg(phi));
    end

    title(title_text);

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
    disp('Animation saved successfully:');
    disp(fullfile(pwd, video_name));
else
    warning('No frames captured.');
end