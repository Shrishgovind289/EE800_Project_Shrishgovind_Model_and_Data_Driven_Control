clc;
clear;
close all;

%% =========================================
% SELF-BALANCING ROBOT
% MODEL-BASED DISSIPATIVITY CONTROLLER
% + DISTURBANCE VALIDATION (-20 deg, -80 deg/s at 10 s)
% =========================================
%
% CATEGORY:
%   Model-Based Controller with Dissipativity
%
% IMPORTANT:
%   This script does NOT use data matrices.
%   This script does NOT recover A_est or B_est from data.
%   The dissipativity controller is designed directly from Ad and Bd.
%
% REQUIREMENTS:
%   - YALMIP installed
%   - SDP solver installed (SDPT3 / SeDuMi / MOSEK)
%
% Closed-loop convention:
%   u(k) = -K_diss * x(k)
%   x(k+1) = Ad*x(k) + Bd*u(k)
%          = (Ad - Bd*K_diss)*x(k)

%% =========================================
% PHYSICAL PARAMETERS
% =========================================
M = 0.5;         % chassis mass
m_body = 0.2;    % body mass
l = 0.3;         % COM distance
I = 0.006;       % body inertia
b = 0.1;         % viscous friction
g = 9.81;        % gravity

%% =========================================
% MODEL CONSTANTS
% =========================================
a = M + m_body;
c = m_body * l;
d = I + m_body * l^2;
Delta = a*d + c^2;

%% =========================================
% CONTINUOUS-TIME MODEL
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

disp('================ CONTINUOUS MODEL ================');
disp('A matrix:');
disp(A);
disp('B matrix:');
disp(B);

%% =========================================
% DISCRETIZATION
% =========================================
Ts = 0.01;

n = size(A,1);
m = size(B,2);
p = size(C,1);

Maug = [A B;
        zeros(m, n+m)];

Md = expm(Maug * Ts);

Ad = Md(1:n,1:n);
Bd = Md(1:n,n+1:n+m);

Cd = C;
Dd = D;

disp('================ DISCRETE MODEL ==================');
disp('Ad matrix:');
disp(Ad);
disp('Bd matrix:');
disp(Bd);

%% =========================================
% MODEL-BASED DISSIPATIVITY CONTROLLER DESIGN
% =========================================
disp('====== MODEL-BASED DISSIPATIVITY CONTROLLER ======');

if ~exist('sdpvar','file')
    error('YALMIP not found. Install YALMIP and SDPT3/SeDuMi/MOSEK.');
end

% Disturbance matrix
Ew = 0.05 * eye(n);
nw = size(Ew,2);

% Performance output: focus mainly on body angle
Cz = diag([0.2 0.05 1 0.2]);
nz = size(Cz,1);

% YALMIP variables
X = sdpvar(n,n,'symmetric');
L = sdpvar(m,n,'full');
gamma = sdpvar(1,1);

% Closed-loop term for u = -Kx
% MODEL-BASED: uses Ad and Bd directly
AclX = Ad*X - Bd*L;

% Discrete-time bounded-real / dissipativity-style LMI
BRL = [ X,            AclX,         Ew,               zeros(n,nz);
        AclX',        X,            zeros(nw,n),      X*Cz';
        Ew',          zeros(nw,n),  gamma*eye(nw),    zeros(nw,nz);
        zeros(nz,n),  Cz*X,         zeros(nz,nw),     gamma*eye(nz) ];

constraints = [];
constraints = [constraints, X >= 1e-5*eye(n)];
constraints = [constraints, gamma >= 1e-4];
constraints = [constraints, BRL >= 1e-5*eye(size(BRL,1))];

% Simple linear bounds on L instead of norm(L,'fro')
Lmax = 5;
constraints = [constraints, -Lmax <= L <= Lmax];

objective = gamma;

ops = sdpsettings('solver','sdpt3','verbose',1);
sol = optimize(constraints, objective, ops);

if sol.problem ~= 0
    error('SDPT3 failed: %s', sol.info);
end

X_val = value(X);
L_val = value(L);
gamma_val = value(gamma);

K_diss = L_val / X_val;

disp('Model-based dissipativity gain K_diss:');
disp(K_diss);
fprintf('Optimal gamma = %.6f\n', gamma_val);

%% =========================================
% CLOSED-LOOP STABILITY CHECK
% =========================================
Acl = Ad - Bd*K_diss;


disp('============ CLOSED-LOOP STABILITY CHECK =========');
disp('Acl = Ad - Bd*K_diss:');
disp(Acl);

cl_eigs = eig(Acl);

disp('Closed-loop eigenvalues:');
disp(cl_eigs);

if all(abs(cl_eigs) < 1)
    disp('Closed-loop discrete-time model-based system is stable.');
else
    disp('Closed-loop discrete-time model-based system is NOT stable.');
end

%% =========================================
% ROBUSTNESS TEST (±10% variation in b and I)
% =========================================
disp('================ ROBUSTNESS TEST =================');

b_values = [0.9*b, 1.1*b];
I_values = [0.9*I, 1.1*I];

u_sat_limit = 20;   % actuator saturation for realistic testing

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

        Acl_test = Ad_test - Bd_test*K_diss;
        eig_test = eig(Acl_test);

        fprintf('\n--- Test Case ---\n');
        fprintf('b = %.4f, I = %.6f\n', b_test, I_test);
        fprintf('Linear closed-loop stable? %d\n', all(abs(eig_test) < 1));
        fprintf('Max eigenvalue magnitude = %.4f\n', max(abs(eig_test)));

        x_sim = zeros(n, 1001);
        x_sim(:,1) = [0;0;deg2rad(45);0];

        u_sim_hist = zeros(m,1000);

        for k = 1:300
            u_cmd = -K_diss * x_sim(:,k);
            u_sat = max(min(u_cmd, u_sat_limit), -u_sat_limit);
            u_sim_hist(:,k) = u_sat;
            x_sim(:,k+1) = Ad_test * x_sim(:,k) + Bd_test * u_sat;
        end

        fprintf('Final angle (deg) = %.4f\n', rad2deg(x_sim(3,end)));
        fprintf('Max |u| with saturation = %.4f\n', max(abs(u_sim_hist(:))));
    end
end

%% =========================================
% DISTURBANCE RESPONSE CHECK
% =========================================
% Disturbance is only injected during validation.

Ntest = 2000;              % 20 seconds, so recovery after 10 s is visible
u_sat_limit = 20;

disturbance_time = 10;     % seconds
dist_step = round(disturbance_time/Ts) + 1;
dist_angle_deg = -20;     % body angle disturbance
dist_angvel_deg = -80;    % angular velocity disturbance
angle_threshold_deg = 1;   % recovery threshold for body angle

x_test = zeros(n, Ntest+1);
u_test = zeros(m, Ntest);
u_cmd_test = zeros(m, Ntest);

x_test(:,1) = [0;
               0;
               deg2rad(45);
               0];

for k = 1:Ntest
    u_cmd = -K_diss * x_test(:,k);
    u_sat = max(min(u_cmd, u_sat_limit), -u_sat_limit);

    u_cmd_test(:,k) = u_cmd;
    u_test(:,k) = u_sat;

    % MODEL-BASED simulation uses Ad and Bd directly
    x_test(:,k+1) = Ad * x_test(:,k) + Bd * u_sat;

    % Apply impulse-like state disturbance at 10 seconds
    if k == dist_step
        x_test(3,k+1) = x_test(3,k+1) + deg2rad(dist_angle_deg);
        x_test(4,k+1) = x_test(4,k+1) + deg2rad(dist_angvel_deg);

        fprintf('\n*** Disturbance injected at t = %.2f s ***\n', disturbance_time);
        fprintf('Added angle disturbance = %.2f deg\n', dist_angle_deg);
        fprintf('Added angular velocity disturbance = %.2f deg/s\n', dist_angvel_deg);
    end
end

% Recovery time calculation: first time after disturbance when |phi| <= threshold
angle_deg = rad2deg(x_test(3,:));
post_dist_idx = dist_step:length(angle_deg);
recovery_idx_rel = find(abs(angle_deg(post_dist_idx)) <= angle_threshold_deg, 1, 'first');

if isempty(recovery_idx_rel)
    recovery_time = NaN;
    disp('Recovery threshold was not reached within the simulation window.');
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
% PLOTS: ALL STATES + CONTROL INPUT
% =========================================

t_state = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

plot_folder = 'Model_Based_Dissipative_Controller_Disturbance_Plots';
if ~exist(plot_folder,'dir')
    mkdir(plot_folder);
end

% 1) Robot / Cart Position
fig1 = figure('Name','Model-Based Dissipative Robot Position');
plot(t_state, x_test(1,:), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('x (m)');
title('Model-Based Dissipativity Controller: Robot / Cart Position');
exportgraphics(fig1, fullfile(plot_folder,'01_robot_position.png'), 'Resolution', 300);

% 2) Linear Velocity
fig2 = figure('Name','Model-Based Dissipative Linear Velocity');
plot(t_state, x_test(2,:), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('xdot (m/s)');
title('Model-Based Dissipativity Controller: Linear Velocity');
exportgraphics(fig2, fullfile(plot_folder,'02_linear_velocity.png'), 'Resolution', 300);

% 3) Body Angle
fig3 = figure('Name','Model-Based Dissipative Body Angle With Disturbance');
plot(t_state, rad2deg(x_test(3,:)), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
yline(angle_threshold_deg, ':', '+1 deg');
yline(-angle_threshold_deg, ':', '-1 deg');
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
title('Model-Based Dissipativity Controller: Body Angle Response With Disturbance');
exportgraphics(fig3, fullfile(plot_folder,'03_body_angle.png'), 'Resolution', 300);

% 4) Angular Velocity
fig4 = figure('Name','Model-Based Dissipative Angular Velocity');
plot(t_state, rad2deg(x_test(4,:)), 'LineWidth', 2);
hold on;
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('\phidot (deg/s)');
title('Model-Based Dissipativity Controller: Angular Velocity');
exportgraphics(fig4, fullfile(plot_folder,'04_angular_velocity.png'), 'Resolution', 300);

% 5) Control Input
fig5 = figure('Name','Model-Based Dissipative Control Input');
plot(t_input, u_cmd_test, '--', 'LineWidth', 1.2);
hold on;
plot(t_input, u_test, 'LineWidth', 2);
xline(disturbance_time, '--', 'Disturbance', 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Control Input');
legend('u command','u saturated','Location','best');
title('Model-Based Dissipativity Controller: Control Input Response');
exportgraphics(fig5, fullfile(plot_folder,'05_control_input.png'), 'Resolution', 300);

fprintf('Saved disturbance plots to folder: %s\n', plot_folder);

%% =========================================
% SAVE EVERYTHING
% =========================================
save('SBR_ModelBased_Dissipativity_Controller_Disturbance_Result.mat', ...
    'A','B','C','D', ...
    'Ad','Bd','Cd','Dd', ...
    'K_diss', ...
    'Acl','cl_eigs', ...
    'x_test','u_test','u_cmd_test', ...
    'disturbance_time','dist_step','dist_angle_deg','dist_angvel_deg', ...
    'angle_threshold_deg','recovery_time','plot_folder', ...
    'Ts','M','m_body','l','I','b','g','a','c','d','Delta', ...
    'Ew','Cz','gamma_val','u_sat_limit');

disp('================ FINISHED ========================');
disp('Saved results to SBR_ModelBased_Dissipativity_Controller_Disturbance_Result.mat');

%% =========================================
% SMOOTHER VIDEO GENERATION
% =========================================
opengl software

video_name = 'self_balancing_robot_model_based_dissipativity_disturbance.mp4';
v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = 100;
open(v);

fig = figure('Visible','on', 'Position', [100 100 900 600]);

wheel_radius = 0.05;
body_length = 0.5;

frame_count = 0;

for k = 1:length(x_test)

    clf(fig); hold on; grid on;
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
    th = linspace(0,2*pi,60);
    plot(wheel_x + wheel_radius*cos(th), ...
         wheel_y + wheel_radius*sin(th), 'b','LineWidth',2);

    % Body
    plot([wheel_x body_x], [wheel_y body_y], 'r','LineWidth',4);

    % COM
    plot(body_x, body_y, 'ko','MarkerFaceColor','k');

    sim_time = (k-1)*Ts;
    if abs(sim_time - disturbance_time) < Ts/2
        plot(wheel_x, wheel_y + 0.18, 'rv', 'MarkerFaceColor','r', 'MarkerSize',8);
    end

    title(sprintf('Model-Based Dissipativity Controller | t = %.2f s | Angle = %.2f deg', sim_time, rad2deg(phi)));
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
fprintf('Video frame rate: %d fps\n', 100);

if frame_count > 0
    disp('Video successfully saved!');
    disp(fullfile(pwd, video_name));
else
    warning('No frames captured — video not created.');
end
