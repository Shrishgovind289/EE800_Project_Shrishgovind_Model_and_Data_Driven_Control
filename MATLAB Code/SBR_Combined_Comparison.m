
clc;
clear;
close all;

%% ============================================================
% COMBINED COMPARISON SCRIPT
% Compares:
%   1) Data-driven LQR controller
%   2) Dissipativity-based LMI controller
%
% Based on your uploaded scripts:
%   - SBR_Data_Gen_Controller.m
%   - SBR_Dissipative_Controller.m
%
% Main idea:
%   - Use ONE shared plant model
%   - Use ONE shared data-collection step
%   - Recover A_est, B_est once
%   - Design BOTH controllers from same recovered model
%   - Simulate BOTH on the same initial condition
%   - Compare angle, position, control effort, and robustness
%
% Requirements for dissipativity part:
%   - YALMIP
%   - SDPT3 (or another SDP solver)
%
% If YALMIP is not available, the script still runs the LQR part.
%% ============================================================

%% =========================
% PHYSICAL PARAMETERS
% =========================
M = 0.5;          % chassis mass
m_body = 0.2;     % body mass
l = 0.3;          % COM distance
I = 0.006;        % body inertia
b = 0.1;          % viscous friction
g = 9.81;         % gravity

%% =========================
% MODEL CONSTANTS
% =========================
a = M + m_body;
c = m_body * l;
d = I + m_body * l^2;
Delta = a*d + c^2;

%% =========================
% CONTINUOUS MODEL
% x = [x; xdot; phi; phidot]
% =========================
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
disp('A matrix:'); disp(A);
disp('B matrix:'); disp(B);

%% =========================
% DISCRETIZATION
% =========================
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
disp('Ad matrix:'); disp(Ad);
disp('Bd matrix:'); disp(Bd);

%% =========================
% DATA COLLECTION SETTINGS
% =========================
T = 1000;
x0_data = [0; 0; deg2rad(45); 0];

%% =========================
% SAFE DATA COLLECTION FEEDBACK
% u = -Ktemp*x + excitation
% =========================
Ktemp = [-1.0  -2.0  18.0  3.5];
Acl_temp = Ad - Bd*Ktemp;

disp('========= DATA COLLECTION CLOSED-LOOP EIGS ========');
disp(eig(Acl_temp));

%% =========================
% PERSISTENTLY EXCITING INPUT
% =========================
rng(1);  % reproducible
u_amp = 0.15;
hold_steps = 8;

num_blocks = ceil(T/hold_steps);
u_blocks = u_amp * sign(randn(1,num_blocks));
u_exc = repelem(u_blocks, hold_steps);
u_exc = u_exc(1:T);

%% =========================
% SIMULATE DATA COLLECTION
% =========================
x = zeros(n, T+1);
y = zeros(p, T);
u = zeros(m, T);
x(:,1) = x0_data;

for k = 1:T
    u(:,k) = -Ktemp * x(:,k) + u_exc(k);
    y(:,k) = Cd * x(:,k) + Dd * u(:,k);
    x(:,k+1) = Ad * x(:,k) + Bd * u(:,k);
end

%% =========================
% BUILD DATA MATRICES
% =========================
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

%% =========================
% RANK CHECK
% =========================
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

%% =========================
% RECOVER MODEL FROM DATA
% X1 = [B A] * [U0; X0]
% =========================
AB_est = X1 * pinv([U0; X0]);
B_est = AB_est(:,1:m);
A_est = AB_est(:,m+1:m+n);

disp('================ RECOVERED MODEL =================');
disp('Recovered A_est from data:'); disp(A_est);
disp('Recovered B_est from data:'); disp(B_est);
fprintf('||A_est - Ad|| = %g\n', norm(A_est - Ad));
fprintf('||B_est - Bd|| = %g\n', norm(B_est - Bd));

%% ============================================================
% CONTROLLER 1: DATA-DRIVEN LQR
%% ============================================================
Q_lqr = diag([10 1 300 1]);
R_lqr = 1;

disp('================ LQR DESIGN ======================');

if exist('dlqr','file') == 2
    disp('Using dlqr from Control System Toolbox...');
    K_lqr = dlqr(A_est, B_est, Q_lqr, R_lqr);
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

Acl_lqr = A_est - B_est*K_lqr;
eig_lqr = eig(Acl_lqr);

disp('K_lqr ='); disp(K_lqr);
disp('eig(A_est - B_est*K_lqr) ='); disp(eig_lqr);

%% ============================================================
% CONTROLLER 2: DISSIPATIVITY-BASED CONTROLLER
%% ============================================================
disp('============= DISSIPATIVITY DESIGN ===============');

dissipativity_available = exist('sdpvar','file') == 2;
K_diss = [];
gamma_val = NaN;
eig_diss = NaN(n,1);

if ~dissipativity_available
    warning(['YALMIP not found. Dissipativity controller will be skipped. ', ...
             'Install YALMIP + SDPT3 to run this part.']);
else
    Ew = 0.05 * eye(n);            % smaller disturbance scaling
    nw = size(Ew,2);

    % same as uploaded code: emphasize angle regulation
    Cz = diag([0.2 0.05 1.0 0.2]);
    nz = size(Cz,1);

    X = sdpvar(n,n,'symmetric');
    L = sdpvar(m,n,'full');
    gamma = sdpvar(1,1);

    AclX = A_est*X - B_est*L;

    BRL = [ X,            AclX,         Ew,               zeros(n,nz);
            AclX',        X,            zeros(nw,n),      X*Cz';
            Ew',          zeros(nw,n),  gamma*eye(nw),    zeros(nw,nz);
            zeros(nz,n),  Cz*X,         zeros(nz,nw),     gamma*eye(nz) ];

    constraints = [];
    constraints = [constraints, X >= 1e-5*eye(n)];
    constraints = [constraints, gamma >= 1e-4];
    constraints = [constraints, BRL >= 1e-5*eye(size(BRL,1))];

    Lmax = 5;
    constraints = [constraints, -Lmax <= L <= Lmax];

    objective = gamma;
    ops = sdpsettings('solver','sdpt3','verbose',1);
    sol = optimize(constraints, objective, ops);

    if sol.problem ~= 0
        warning('Dissipativity controller failed: %s', sol.info);
    else
        X_val = value(X);
        L_val = value(L);
        gamma_val = value(gamma);

        K_diss = L_val / X_val;
        Acl_diss = A_est - B_est*K_diss;
        eig_diss = eig(Acl_diss);

        disp('K_diss ='); disp(K_diss);
        fprintf('Optimal gamma = %.6f\n', gamma_val);
        disp('eig(A_est - B_est*K_diss) ='); disp(eig_diss);
    end
end

%% ============================================================
% SAME TEST CONDITIONS FOR BOTH CONTROLLERS
%% ============================================================
Ntest = 1000;
u_sat_limit = 200;
x0_test = [0; 0; deg2rad(45); 0];

results = struct();

% ---------- LQR SIMULATION ----------
x_lqr = zeros(n, Ntest+1);
u_cmd_lqr = zeros(m, Ntest);
u_sat_lqr = zeros(m, Ntest);
x_lqr(:,1) = x0_test;

for k = 1:Ntest
    u_cmd = -K_lqr * x_lqr(:,k);
    u_sat = max(min(u_cmd, u_sat_limit), -u_sat_limit);

    u_cmd_lqr(:,k) = u_cmd;
    u_sat_lqr(:,k) = u_sat;

    x_lqr(:,k+1) = A_est*x_lqr(:,k) + B_est*u_sat;
end

results.lqr.final_angle_deg = rad2deg(x_lqr(3,end));
results.lqr.final_state_norm = norm(x_lqr(:,end));
results.lqr.max_u_cmd = max(abs(u_cmd_lqr(:)));
results.lqr.max_u_sat = max(abs(u_sat_lqr(:)));

% settling time estimate on angle
phi_lqr_deg = rad2deg(x_lqr(3,:));
results.lqr.settle_time = estimate_settling_time(phi_lqr_deg, Ts, 0.5);

% ---------- DISSIPATIVITY SIMULATION ----------
if ~isempty(K_diss)
    x_diss = zeros(n, Ntest+1);
    u_cmd_diss = zeros(m, Ntest);
    u_sat_diss = zeros(m, Ntest);
    x_diss(:,1) = x0_test;

    for k = 1:Ntest
        u_cmd = -K_diss * x_diss(:,k);
        u_sat = max(min(u_cmd, u_sat_limit), -u_sat_limit);

        u_cmd_diss(:,k) = u_cmd;
        u_sat_diss(:,k) = u_sat;

        x_diss(:,k+1) = A_est*x_diss(:,k) + B_est*u_sat;
    end

    results.diss.final_angle_deg = rad2deg(x_diss(3,end));
    results.diss.final_state_norm = norm(x_diss(:,end));
    results.diss.max_u_cmd = max(abs(u_cmd_diss(:)));
    results.diss.max_u_sat = max(abs(u_sat_diss(:)));
    phi_diss_deg = rad2deg(x_diss(3,:));
    results.diss.settle_time = estimate_settling_time(phi_diss_deg, Ts, 0.5);
else
    x_diss = [];
    u_cmd_diss = [];
    u_sat_diss = [];
end

%% ============================================================
% PRINT COMPARISON
%% ============================================================
disp('================ FINAL COMPARISON =================');
fprintf('\nLQR controller:\n');
fprintf('  Final angle (deg) = %.6f\n', results.lqr.final_angle_deg);
fprintf('  Final state norm  = %.6e\n', results.lqr.final_state_norm);
fprintf('  Max |u_cmd|       = %.6f\n', results.lqr.max_u_cmd);
fprintf('  Max |u_sat|       = %.6f\n', results.lqr.max_u_sat);
fprintf('  Settling time (s) = %.4f\n', results.lqr.settle_time);

if ~isempty(K_diss)
    fprintf('\nDissipativity controller:\n');
    fprintf('  Final angle (deg) = %.6f\n', results.diss.final_angle_deg);
    fprintf('  Final state norm  = %.6e\n', results.diss.final_state_norm);
    fprintf('  Max |u_cmd|       = %.6f\n', results.diss.max_u_cmd);
    fprintf('  Max |u_sat|       = %.6f\n', results.diss.max_u_sat);
    fprintf('  Settling time (s) = %.4f\n', results.diss.settle_time);
    fprintf('  Gamma             = %.6f\n', gamma_val);
end

%% ============================================================
% ROBUSTNESS TEST FOR BOTH
%% ============================================================
disp('================ ROBUSTNESS TEST ==================');

b_values = [0.9*b, 1.1*b];
I_values = [0.9*I, 1.1*I];

robust_table = [];

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

        % LQR
        Acl_test_lqr = Ad_test - Bd_test*K_lqr;
        eig_lqr_test = eig(Acl_test_lqr);
        stable_lqr = all(abs(eig_lqr_test) < 1);

        xsim = zeros(n,301);
        xsim(:,1) = x0_test;
        for k = 1:300
            ucmd = -K_lqr*xsim(:,k);
            usat = max(min(ucmd, u_sat_limit), -u_sat_limit);
            xsim(:,k+1) = Ad_test*xsim(:,k) + Bd_test*usat;
        end
        final_angle_lqr = rad2deg(xsim(3,end));

        if ~isempty(K_diss)
            Acl_test_diss = Ad_test - Bd_test*K_diss;
            eig_diss_test = eig(Acl_test_diss);
            stable_diss = all(abs(eig_diss_test) < 1);

            xsim2 = zeros(n,301);
            xsim2(:,1) = x0_test;
            for k = 1:300
                ucmd = -K_diss*xsim2(:,k);
                usat = max(min(ucmd, u_sat_limit), -u_sat_limit);
                xsim2(:,k+1) = Ad_test*xsim2(:,k) + Bd_test*usat;
            end
            final_angle_diss = rad2deg(xsim2(3,end));
            maxeig_diss = max(abs(eig_diss_test));
        else
            stable_diss = NaN;
            final_angle_diss = NaN;
            maxeig_diss = NaN;
        end

        maxeig_lqr = max(abs(eig_lqr_test));

        robust_table = [robust_table;
            b_test, I_test, stable_lqr, maxeig_lqr, final_angle_lqr, stable_diss, maxeig_diss, final_angle_diss]; %#ok<AGROW>
    end
end

disp('Columns = [b, I, stable_lqr, maxeig_lqr, final_angle_lqr, stable_diss, maxeig_diss, final_angle_diss]');
disp(robust_table);

%% ============================================================
% PLOTS
%% ============================================================
t_state = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

figure('Name','Angle Comparison');
plot(t_state, rad2deg(x_lqr(3,:)), 'LineWidth', 2); hold on;
if ~isempty(K_diss)
    plot(t_state, rad2deg(x_diss(3,:)), 'LineWidth', 2);
    legend('LQR','Dissipativity','Location','best');
else
    legend('LQR','Location','best');
end
grid on;
xlabel('Time (s)');
ylabel('\phi (deg)');
title('Body Angle Comparison');

figure('Name','Position Comparison');
plot(t_state, x_lqr(1,:), 'LineWidth', 2); hold on;
if ~isempty(K_diss)
    plot(t_state, x_diss(1,:), 'LineWidth', 2);
    legend('LQR','Dissipativity','Location','best');
else
    legend('LQR','Location','best');
end
grid on;
xlabel('Time (s)');
ylabel('Cart / body position x (m)');
title('Position Comparison');

figure('Name','Control Comparison');

plot(t_input, u_sat_lqr, 'LineWidth', 2); hold on;

if ~isempty(K_diss)
    plot(t_input, u_sat_diss, 'LineWidth', 2);
    legend('LQR','Dissipativity','Location','best');
    u_all = [u_sat_lqr, u_sat_diss];
else
    legend('LQR','Location','best');
    u_all = u_sat_lqr;
end

grid on;
xlabel('Time (s)');
ylabel('Saturated control input');
title('Control Input Comparison');

% Auto-scale with margin
ylim([min(u_all(:)) max(u_all(:))] * 1.1);

figure('Name','Angular Velocity Comparison');
plot(t_state, rad2deg(x_lqr(4,:)), 'LineWidth', 2); hold on;
if ~isempty(K_diss)
    plot(t_state, rad2deg(x_diss(4,:)), 'LineWidth', 2);
    legend('LQR','Dissipativity','Location','best');
else
    legend('LQR','Location','best');
end
grid on;
xlabel('Time (s)');
ylabel('\dot{\phi} (deg/s)');
title('Angular Velocity Comparison');

%% ============================================================
% SAVE
%% ============================================================
save('SBR_Combined_Comparison.mat', ...
    'A','B','C','D', ...
    'Ad','Bd','Cd','Dd', ...
    'A_est','B_est', ...
    'X0','X1','U0','Y0','Y1', ...
    'Ktemp','K_lqr','K_diss', ...
    'eig_lqr','eig_diss','gamma_val', ...
    'x_lqr','u_cmd_lqr','u_sat_lqr', ...
    'x_diss','u_cmd_diss','u_sat_diss', ...
    'results','robust_table', ...
    'Ts','M','m_body','l','I','b','g');

disp('================ DONE ============================');
disp('Saved: SBR_Combined_Comparison.mat');

%% ============================================================
% LOCAL FUNCTION
%% ============================================================
function t_settle = estimate_settling_time(signal, Ts, tol)
    % signal: vector
    % tol: band around final value
    final_val = signal(end);
    idx = length(signal);

    for k = 1:length(signal)
        tail = signal(k:end);
        if all(abs(tail - final_val) <= tol)
            idx = k;
            break;
        end
    end

    t_settle = (idx-1) * Ts;
end
