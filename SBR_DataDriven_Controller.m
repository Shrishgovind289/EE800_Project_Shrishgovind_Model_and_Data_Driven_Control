clc;
clear;
close all;

%% =========================================
% PURE DATA-DRIVEN CONTROLLER DESIGN
% Uses only SBR_Data_Matrices.mat
% No A_est, B_est, Ad, Bd, A, or B are used
% =========================================

%% Load data
load('SBR_Data_Matrices.mat');

disp('================ LOADED DATA =====================');
fprintf('size(X0) = [%d %d]\n', size(X0,1), size(X0,2));
fprintf('size(X1) = [%d %d]\n', size(X1,1), size(X1,2));
fprintf('size(U0) = [%d %d]\n', size(U0,1), size(U0,2));

%% Dimensions
n = size(X0,1);      % number of states
T = size(X0,2);      % number of data samples
m = size(U0,1);      % number of inputs

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
% convention: u = Kx
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

ops = sdpsettings('solver','sdpt3','verbose',1,'sdpt3.maxit',200); %sdpsettings('solver','sdpt3','verbose',1);
sol = optimize(constraints, objective, ops);

if sol.problem ~= 0
    error('Data-driven LMI failed: %s', sol.info);
end

%% Extract solution
Q_val = value(Q);
P_val = value(P_sym);

K_data = U0c * Q_val / P_val;

%% =========================================
% DATA-DRIVEN CLOSED-LOOP MATRIX
%
% This is not A_est or B_est.
% It is the closed-loop map directly from data.
% =========================================

Acl_dd = X1c * Q_val / P_val;

disp('============ DATA-DRIVEN CLOSED-LOOP CHECK =======');
disp('Closed-loop eigenvalues from data:');
cl_eigs = eig(Acl_dd);
disp(cl_eigs);

if all(abs(cl_eigs) < 1)
    disp('Data-driven closed-loop system is stable.');
else
    disp('Data-driven closed-loop system is NOT stable.');
end

%% =========================================
% CLOSED-LOOP RESPONSE USING DATA-DRIVEN MAP
% x(k+1) = Acl_dd*x(k)
% =========================================

Ntest = 1000;

x_test = zeros(n,Ntest+1);
u_test = zeros(m,Ntest);

x_test(:,1) = [0;
               0;
               deg2rad(45);
               0];

for k = 1:Ntest
    u_test(:,k) = K_data * x_test(:,k);
    x_test(:,k+1) = Acl_dd * x_test(:,k);
end

t = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

disp('================ TEST RESPONSE ===================');
fprintf('Initial angle (deg) = %.4f\n', rad2deg(x_test(3,1)));
fprintf('Final angle (deg)   = %.6f\n', rad2deg(x_test(3,end)));
fprintf('Final state norm    = %.6e\n', norm(x_test(:,end)));
fprintf('Max |u|             = %.6f\n', max(abs(u_test(:))));

%% =========================================
% PLOTS
% =========================================

figure('Name','Data-Driven State Response');
plot(t, x_test(1,:), 'LineWidth', 2);
grid on;
ylabel('x (m)');
title('Data-Driven Cart Position');

figure('Name','Data-Driven State Response');
plot(t, x_test(2,:), 'LineWidth', 2);
grid on;
ylabel('xdot (m/s)');
title('Data-Driven Linear Velocity');

figure('Name','Data-Driven State Response');
plot(t, rad2deg(x_test(3,:)), 'LineWidth', 2);
grid on;
ylabel('\phi (deg)');
title('Data-Driven Body Angle');

figure('Name','Data-Driven State Response');
plot(t, rad2deg(x_test(4,:)), 'LineWidth', 2);
grid on;
ylabel('\phidot (deg/s)');
xlabel('Time (s)');
title('Data-Driven Angular Velocity');

figure('Name','Data-Driven State Response');
figure('Name','Data-Driven Control Input');
plot(t_input, u_test, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Control input u');
title('Data-Driven Control Input');

%% =========================================
% SAVE RESULTS
% =========================================

save('SBR_DataDriven_Controller_Result.mat', ...
    'K_data', ...
    'Q_val', 'P_val', ...
    'Acl_dd', 'cl_eigs', ...
    'x_test', 'u_test', ...
    'Ts');

disp('================ RESULTS SAVED ===================');
disp('Saved result as SBR_DataDriven_Controller_Result.mat');

%% =========================================
% SHOW + SAVE DATA-DRIVEN ROBOT ANIMATION (FIXED AXIS)
% =========================================

disp('========== GENERATING AND SAVING ANIMATION ==========');

opengl software   % safer rendering

video_name = 'data_driven_robot_animation.mp4';

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

    % FIXED VIEW
    xlim([-1 1]);
    ylim([-0.1 1]);

    % States
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