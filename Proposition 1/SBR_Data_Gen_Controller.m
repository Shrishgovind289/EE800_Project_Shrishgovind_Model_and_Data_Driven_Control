clc;
clear;
close all;

%% =========================================
% SELF-BALANCING ROBOT
% DATA GENERATION + DATA-BASED MODEL RECOVERY
% + CONTROLLER DESIGN IN ONE SCRIPT
% =========================================

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
% DATA COLLECTION SETTINGS
% =========================================
T = 1000;   % number of samples

x0 = [0;
      0;
      deg2rad(45);
      0];

%% =========================================
% TEMPORARY STABILIZING FEEDBACK FOR SAFE DATA COLLECTION
% u = Ktemp*x + excitation
% =========================================
Ktemp = [1.0  2.0  -18.0  -3.5];

Acl_temp = Ad + Bd*Ktemp;

disp('========= DATA COLLECTION CLOSED-LOOP EIGS ========');
disp(eig(Acl_temp));

%% =========================================
% PERSISTENTLY EXCITING INPUT
% =========================================
u_amp = 0.15;
hold_steps = 8;

num_blocks = ceil(T/hold_steps);
u_blocks = u_amp * sign(randn(1,num_blocks));
u_exc = repelem(u_blocks, hold_steps);
u_exc = u_exc(1:T);

%% =========================================
% SIMULATION STORAGE
% =========================================
x = zeros(n, T+1);
y = zeros(p, T);
u = zeros(m, T);

x(:,1) = x0;

%% =========================================
% SIMULATE DATA COLLECTION
% =========================================
for k = 1:T
    u(:,k) = Ktemp * x(:,k) + u_exc(k);
    y(:,k) = Cd * x(:,k) + Dd * u(:,k);
    x(:,k+1) = Ad * x(:,k) + Bd * u(:,k);
end

%% =========================================
% BUILD DATA MATRICES
% =========================================
X0 = x(:,1:T);
X1 = x(:,2:T+1);
U0 = u(:,1:T);
Y0 = y(:,1:T);

Y1 = zeros(p,T);
Y1(:,1:T-1) = y(:,2:T);
Y1(:,T) = Cd*x(:,T+1) + Dd*u(:,T);

disp('================ DATA MATRIX SIZES ===============');
fprintf('size(X0) = [%d %d]\n', size(X0,1), size(X0,2));
fprintf('size(X1) = [%d %d]\n', size(X1,1), size(X1,2));
fprintf('size(U0) = [%d %d]\n', size(U0,1), size(U0,2));
fprintf('size(Y0) = [%d %d]\n', size(Y0,1), size(Y0,2));

%% =========================================
% RANK CHECK
% =========================================
rank_data = rank([U0; X0]);
required_rank = n + m;

disp('================ RANK CHECK ======================');
fprintf('rank([U0; X0]) = %d\n', rank_data);
fprintf('required rank  = %d\n', required_rank);

if rank_data < required_rank
    error('Data is NOT rich enough for data-driven control.');
else
    disp('Data matrix has full rank. Good for data-driven control.');
end

%% =========================================
% SAVE DATA MATRICES FOR PURE DATA-DRIVEN CONTROLLER
% =========================================
% These matrices will be used in the pure data-driven controller script.
% The next controller file should load this .mat file and use X0, X1, U0
% directly, without using A_est, B_est, Ad, or Bd.

data_file_name = 'SBR_Data_Matrices.mat';

save(data_file_name, ...
    'X0', 'X1', 'U0', ...
    'Y0', 'Y1', ...
    'Ts', ...
    'n', 'm', 'p', ...
    'x0', ...
    'Ktemp', ...
    'u_exc');

disp('================ DATA SAVED ======================');
fprintf('Saved data matrices to: %s\n', data_file_name);
disp('This file can now be used for pure data-driven controller design.');
%% =========================================
% RECOVER MODEL FROM DATA
% X1 = [B A] * [U0; X0]
% =========================================
AB_est = X1 * pinv([U0; X0]);

B_est = AB_est(:,1:m);
A_est = AB_est(:,m+1:m+n);

disp('================ RECOVERED MODEL =================');
disp('Recovered A_est from data:');
disp(A_est);

disp('Recovered B_est from data:');
disp(B_est);

fprintf('||A_est - Ad|| = %g\n', norm(A_est - Ad));
fprintf('||B_est - Bd|| = %g\n', norm(B_est - Bd));

%% =========================================
% DISCRETE LQR DESIGN
% If dlqr exists, use it. Otherwise use iterative Riccati.
% =========================================
Q_lqr = diag([10 1 300 1]);
R_lqr = 1;

disp('================ CONTROLLER DESIGN ===============');

if exist('dlqr','file') == 2
    disp('Using dlqr from Control System Toolbox...');
    K_lqr_neg = dlqr(A_est, B_est, Q_lqr, R_lqr);
    K_lqr = -K_lqr_neg;   % convert to u = Kx
else
    disp('dlqr not found. Using iterative Riccati solution...');

    P = Q_lqr;
    max_iter = 1000;
    tol = 1e-9;

    for k = 1:max_iter
        P_next = A_est' * P * A_est ...
            - A_est' * P * B_est * inv(R_lqr + B_est' * P * B_est) * B_est' * P * A_est ...
            + Q_lqr;

        if norm(P_next - P, 'fro') < tol
            P = P_next;
            fprintf('Riccati converged in %d iterations.\n', k);
            break;
        end

        P = P_next;
    end

    if k == max_iter
        warning('Riccati iteration did not fully converge.');
    end

    K_lqr = inv(R_lqr + B_est' * P * B_est) * (B_est' * P * A_est);
end

disp('LQR gain K_lqr:');
disp(K_lqr);

disp('LQR gain K_lqr:');
disp(K_lqr);

%% =========================================
% ROBUSTNESS TEST (±5% variation in b and I)
% =========================================
disp('================ ROBUSTNESS TEST =================');

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

        Acl_test = Ad_test + Bd_test*K_lqr;
        eig_test = eig(Acl_test);

        fprintf('\n--- Test Case ---\n');
        fprintf('b = %.4f, I = %.6f\n', b_test, I_test);
        fprintf('Stable? %d\n', all(abs(eig_test) < 1));
        fprintf('Max eigenvalue magnitude = %.4f\n', max(abs(eig_test)));

        % quick simulation
        x_sim = zeros(n, 1001);
        x_sim(:,1) = [0;0;deg2rad(45);0];

        for k = 1:300
            u_sim = K_lqr * x_sim(:,k);
            x_sim(:,k+1) = Ad_test * x_sim(:,k) + Bd_test * u_sim;
        end

        fprintf('Final angle (deg) = %.4f\n', rad2deg(x_sim(3,end)));
    end
end

%% =========================================
% CLOSED-LOOP STABILITY CHECK
% Standard convention: u = K*x
% =========================================
Acl = A_est + B_est*K_lqr;

disp('============ CLOSED-LOOP STABILITY CHECK =========');
disp('Acl = A_est + B_est*K_lqr:');
disp(Acl);

cl_eigs = eig(Acl);

disp('Closed-loop eigenvalues:');
disp(cl_eigs);

if all(abs(cl_eigs) < 1)
    disp('Closed-loop discrete-time system is stable.');
else
    disp('Closed-loop discrete-time system is NOT stable.');
end

%% =========================================
% OPTIONAL COMPARISON WITH TRUE DISCRETE MODEL
% =========================================
disp('======== EIGENVALUES USING TRUE DISCRETE MODEL ========');
disp('eig(Ad + Bd*K_lqr):');
disp(eig(Ad + Bd*K_lqr));

%% =========================================
% CLOSED-LOOP TEST SIMULATION
% =========================================
Ntest = 1000;
u_sat_limit = 20;

x_test = zeros(n, Ntest+1);
u_test = zeros(m, Ntest);
u_cmd_test = zeros(m, Ntest);

% Initial condition (45 deg tilt)
x_test(:,1) = [0;
               0;
               deg2rad(45);
               0];

for k = 1:Ntest
    % Raw control
    u_cmd = K_lqr * x_test(:,k);

    % Saturation (real actuator)
    u_sat = max(min(u_cmd, u_sat_limit), -u_sat_limit);

    % Store
    u_cmd_test(:,k) = u_cmd;
    u_test(:,k) = u_sat;

    % System update
    x_test(:,k+1) = A_est * x_test(:,k) + B_est * u_sat;
end

disp('================ TEST RESPONSE ===================');
fprintf('Initial angle (deg) = %.4f\n', rad2deg(x_test(3,1)));
fprintf('Final angle (deg)   = %.6f\n', rad2deg(x_test(3,end)));
fprintf('Final state norm    = %.6e\n', norm(x_test(:,end)));
fprintf('Max |u_cmd|         = %.6f\n', max(abs(u_cmd_test(:))));
fprintf('Max |u_sat|         = %.6f\n', max(abs(u_test(:))));
%% =========================================
% SAVE EVERYTHING
% =========================================
save('SBR_DataDriven_AllInOne_Result.mat', ...
    'A','B','C','D', ...
    'Ad','Bd','Cd','Dd', ...
    'X0','X1','U0','Y0','Y1', ...
    'A_est','B_est', ...
    'Ktemp','K_lqr', ...
    'Acl','cl_eigs', ...
    'x','y','u','u_exc','x_test','u_test', ...
    'Ts','M','m_body','l','I','b','g','a','c','d','Delta');

disp('================ FINISHED ========================');
disp('Saved results to SBR_DataDriven_AllInOne_Result.mat');

%% =========================================
% PLOTS: BODY ANGLE + CONTROL INPUT
% =========================================

t_state = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

figure;
plot(t_state, rad2deg(x_test(3,:)), 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
title('Body Angle Response');

figure;
plot(t_input, u_cmd_test, '--', 'LineWidth', 1.2); hold on;
plot(t_input, u_test, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Control Input');
legend('u command','u saturated','Location','best');
title('Control Input Response');

%% =========================================
% SAFE VIDEO GENERATION (WILL WORK)
% =========================================

opengl software   % force rendering

v = VideoWriter('self_balancing_robot.mp4','MPEG-4');
v.FrameRate = 10;
open(v);

fig = figure('Visible','on');

wheel_radius = 0.05;
body_length = 0.5;

frame_count = 0;

for k = 1:5:length(x_test)

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
    theta = linspace(0,2*pi,50);
    plot(wheel_x + wheel_radius*cos(theta), ...
         wheel_y + wheel_radius*sin(theta), 'b','LineWidth',2);

    % Body
    plot([wheel_x body_x], [wheel_y body_y], 'r','LineWidth',4);

    % COM
    plot(body_x, body_y, 'ko','MarkerFaceColor','k');

    title(sprintf('Time step = %d | Angle = %.2f deg', ...
          k, rad2deg(phi)));

    drawnow;   % REQUIRED

    frame = getframe(fig);

    % Safety check
    if ~isempty(frame.cdata)
        writeVideo(v, frame);
        frame_count = frame_count + 1;
    end
end

close(v);

fprintf('Frames written: %d\n', frame_count);

if frame_count > 0
    disp('Video successfully saved!');
    disp(fullfile(pwd, 'self_balancing_robot.mp4'));
else
    warning('No frames captured — video not created.');
end
