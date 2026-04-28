clc;
clear;
close all;

%% =========================================
% PURE DATA-DRIVEN CONTROLLER DESIGN
% Uses only SBR_Data_Matrices.mat
% Convention: u(k) = K_data*x(k)
% =========================================

%% Load data
load('SBR_Data_Matrices.mat');

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
    error('Data is not rich enough. Need persistently exciting data.');
else
    disp('Data is rich enough for pure data-driven control.');
end

%% Check YALMIP
if ~exist('sdpvar','file')
    error('YALMIP not found. Install YALMIP first.');
end

%% =========================================
% COMPRESS DATA TO INDEPENDENT COLUMNS
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
%
% LMI:
% [X0Q   X1Q
%  Q'X1' X0Q] > 0
%
% Controller:
% K_data = U0Q(X0Q)^(-1)
% =========================================

Q = sdpvar(Tc,n,'full');

eps_lmi = 1e-8;

P_data = X0c * Q;
X1Q = X1c * Q;

P_sym = 0.5*(P_data + P_data');

constraints = [];
constraints = [constraints, P_data - P_data' == 0];
constraints = [constraints, P_sym >= eps_lmi*eye(n)];

constraints = [constraints, ...
    [P_sym, X1Q;
     X1Q',  P_sym] >= eps_lmi*eye(2*n)];

objective = [];

ops = sdpsettings('solver','sdpt3','verbose',1,'sdpt3.maxit',200);
sol = optimize(constraints, objective, ops);

if sol.problem ~= 0
    error('Data-driven LMI failed: %s', sol.info);
end

%% Extract controller
Q_val = value(Q);
P_val = value(P_sym);

K_data = U0c * Q_val / P_val;

Acl_dd = X1c * Q_val / P_val;

disp('================ DATA-DRIVEN GAIN ================');
disp('K_data for convention u = K_data*x:');
disp(K_data);

disp('============ DATA-DRIVEN CLOSED-LOOP CHECK =======');
cl_eigs = eig(Acl_dd);
disp('Closed-loop eigenvalues from data:');
disp(cl_eigs);

if all(abs(cl_eigs) < 1)
    disp('Data-driven closed-loop system is stable.');
else
    disp('Data-driven closed-loop system is NOT stable.');
end

%% =========================================
% CLOSED-LOOP RESPONSE WITH DISTURBANCE AT 10s
% =========================================

Ntest = 2000;          % 20 seconds for Ts = 0.01
u_sat_limit = 20;      % actuator saturation

dist_time = 10;        % seconds
dist_step = round(dist_time/Ts);

dist_angle_deg = -20;   % disturbance added to body angle
dist_vec = zeros(n,1);
dist_vec(3) = deg2rad(dist_angle_deg);

x_test = zeros(n,Ntest+1);
u_test = zeros(m,Ntest);
u_cmd_test = zeros(m,Ntest);

x_test(:,1) = [0;
               0;
               deg2rad(45);
               0];

for k = 1:Ntest

    % Control law: u = Kx
    u_cmd = K_data * x_test(:,k);

    % Saturation
    u_sat = max(min(u_cmd, u_sat_limit), -u_sat_limit);

    u_cmd_test(:,k) = u_cmd;
    u_test(:,k) = u_sat;

    % Closed-loop update using data-driven map
    x_test(:,k+1) = Acl_dd * x_test(:,k);

    % Apply impulse angle disturbance at 10s
    if k == dist_step
        x_test(:,k+1) = x_test(:,k+1) + dist_vec;
    end
end

t = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

%% =========================================
% RECOVERY TIME CALCULATION
% =========================================

angle_deg = rad2deg(x_test(3,:));

recovery_threshold_deg = 1;
hold_time = 0.5;
hold_steps = round(hold_time/Ts);

recovery_time = NaN;

for k = dist_step:length(angle_deg)-hold_steps
    if all(abs(angle_deg(k:k+hold_steps)) < recovery_threshold_deg)
        recovery_time = (k - dist_step)*Ts;
        break;
    end
end

disp('================ DISTURBANCE TEST =================');
fprintf('Disturbance applied at %.2f s\n', dist_time);
fprintf('Disturbance magnitude = %.2f deg added to phi\n', dist_angle_deg);

if isnan(recovery_time)
    disp('Recovery time: not recovered within simulation time.');
else
    fprintf('Recovery time = %.4f seconds\n', recovery_time);
end

disp('================ TEST RESPONSE ===================');
fprintf('Initial angle (deg) = %.4f\n', rad2deg(x_test(3,1)));
fprintf('Final angle (deg)   = %.6f\n', rad2deg(x_test(3,end)));
fprintf('Final state norm    = %.6e\n', norm(x_test(:,end)));
fprintf('Max |u_cmd|         = %.6f\n', max(abs(u_cmd_test(:))));
fprintf('Max |u_sat|         = %.6f\n', max(abs(u_test(:))));

%% =========================================
% PLOTS
% =========================================

figure('Name','Data-Driven Cart Position');
plot(t, x_test(1,:), 'LineWidth', 2);
xline(dist_time,'--','Disturbance','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('x (m)');
title('Data-Driven Cart Position');

figure('Name','Data-Driven Linear Velocity');
plot(t, x_test(2,:), 'LineWidth', 2);
xline(dist_time,'--','Disturbance','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('xdot (m/s)');
title('Data-Driven Linear Velocity');

figure('Name','Data-Driven Body Angle With Disturbance');
plot(t, rad2deg(x_test(3,:)), 'LineWidth', 2);
hold on;
xline(dist_time,'--','Disturbance','LineWidth',1.5);
yline(recovery_threshold_deg,':','+1 deg');
yline(-recovery_threshold_deg,':','-1 deg');
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
title('Data-Driven Body Angle Response With Disturbance');

figure('Name','Data-Driven Angular Velocity');
plot(t, rad2deg(x_test(4,:)), 'LineWidth', 2);
xline(dist_time,'--','Disturbance','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('\phidot (deg/s)');
title('Data-Driven Angular Velocity');

figure('Name','Data-Driven Control Input');
plot(t_input, u_cmd_test, '--', 'LineWidth', 1.2);
hold on;
plot(t_input, u_test, 'LineWidth', 2);
xline(dist_time,'--','Disturbance','LineWidth',1.5);
grid on;
xlabel('Time (s)');
ylabel('Control input u');
legend('u command','u saturated','Disturbance','Location','best');
title('Data-Driven Control Input');

%% =========================================
% SAVE RESULTS
% =========================================

save('SBR_DataDriven_Controller_Result.mat', ...
    'K_data', ...
    'Q_val', 'P_val', ...
    'Acl_dd', 'cl_eigs', ...
    'x_test', 'u_test', 'u_cmd_test', ...
    'Ts', ...
    'dist_time', 'dist_angle_deg', ...
    'recovery_time');

disp('================ RESULTS SAVED ===================');
disp('Saved result as SBR_DataDriven_Controller_Result.mat');

%% =========================================
% SHOW + SAVE DATA-DRIVEN ROBOT ANIMATION
% FIXED AXIS
% =========================================

disp('========== GENERATING AND SAVING ANIMATION ==========');

opengl software

video_name = 'data_driven_robot_disturbance_animation.mp4';

v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = 10;
open(v);

fig = figure('Name','Data-Driven Self-Balancing Robot Animation', ...
             'Visible','on');

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
    phi   = x_test(3,k);

    wheel_x = x_pos;
    wheel_y = wheel_radius;

    body_x = wheel_x + body_length*sin(phi);
    body_y = wheel_y + body_length*cos(phi);

    plot([-5 5],[0 0],'k','LineWidth',2);

    theta = linspace(0,2*pi,60);
    plot(wheel_x + wheel_radius*cos(theta), ...
         wheel_y + wheel_radius*sin(theta), ...
         'b','LineWidth',2);

    plot([wheel_x body_x], [wheel_y body_y], ...
         'r','LineWidth',4);

    plot(body_x, body_y, 'ko','MarkerFaceColor','k');

    if abs((k-1)*Ts - dist_time) < 0.05
        text(-0.85,0.9,'DISTURBANCE!', ...
            'FontSize',12,'FontWeight','bold','Color','r');
    end

    title(sprintf('Data-Driven Robot | Time = %.2f s | Angle = %.2f deg', ...
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
    disp('Animation saved successfully:');
    disp(fullfile(pwd, video_name));
else
    warning('No frames captured. Video was not saved correctly.');
end