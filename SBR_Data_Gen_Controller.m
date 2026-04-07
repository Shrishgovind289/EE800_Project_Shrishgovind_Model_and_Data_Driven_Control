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
T = 800;   % number of samples

x0 = [0;
      0;
      deg2rad(5);
      0];

%% =========================================
% TEMPORARY STABILIZING FEEDBACK FOR SAFE DATA COLLECTION
% u = -Ktemp*x + excitation
% =========================================
Ktemp = [-1.0  -2.0  18.0  3.5];

Acl_temp = Ad - Bd*Ktemp;

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
    u(:,k) = -Ktemp * x(:,k) + u_exc(k);
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
Q_lqr = diag([10 1 100 1]);
R_lqr = 1;

disp('================ CONTROLLER DESIGN ===============');

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

disp('LQR gain K_lqr:');
disp(K_lqr);

%% =========================================
% CLOSED-LOOP STABILITY CHECK
% Standard convention: u = -K*x
% =========================================
Acl = A_est - B_est*K_lqr;

disp('============ CLOSED-LOOP STABILITY CHECK =========');
disp('Acl = A_est - B_est*K_lqr:');
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
disp('eig(Ad - Bd*K_lqr):');
disp(eig(Ad - Bd*K_lqr));

%% =========================================
% SIMPLE NUMERICAL RESPONSE CHECK
% NO PLOTS, ONLY FINAL VALUES
% =========================================
x_test = zeros(n, 201);
u_test = zeros(m, 200);

x_test(:,1) = [0;
               0;
               deg2rad(10);
               0];

for k = 1:200
    u_test(:,k) = -K_lqr * x_test(:,k);
    x_test(:,k+1) = A_est * x_test(:,k) + B_est * u_test(:,k);
end

disp('================ TEST RESPONSE ===================');
fprintf('Initial angle (deg) = %.4f\n', rad2deg(x_test(3,1)));
fprintf('Final angle (deg)   = %.6f\n', rad2deg(x_test(3,end)));
fprintf('Final state norm    = %.6e\n', norm(x_test(:,end)));
fprintf('Max |u| during test = %.6f\n', max(abs(u_test(:))));

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