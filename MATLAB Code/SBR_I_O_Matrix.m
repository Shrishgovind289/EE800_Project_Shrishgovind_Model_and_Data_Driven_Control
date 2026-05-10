clc;
clear;
close all;

%% =========================================
% PHYSICAL PARAMETERS
% =========================================
M = 0.5;         % chassis mass
m_body = 0.2;    % body mass
l = 0.3;         % center of mass distance
I = 0.006;       % body inertia
b = 0.1;         % viscous friction
g = 9.81;        % gravity

%% =========================================
% MODEL CONSTANTS
% =========================================
a = M + m_body;
c = m_body*l;
d = I + m_body*l^2;
Delta = a*d + c^2;

%% =========================================
% CONTINUOUS-TIME STATE-SPACE MODEL
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

C = eye(4);              % full-state output
D = zeros(4,1);

disp('A matrix:');
disp(A);

disp('B matrix:');
disp(B);

%% =========================================
% DISCRETIZATION
% =========================================
Ts = 0.01;               % sampling time

n = size(A,1);
m = size(B,2);
p = size(C,1);

Maug = [A B;
        zeros(m,n+m)];

Md = expm(Maug*Ts);

Ad = Md(1:n,1:n);
Bd = Md(1:n,n+1:n+m);

Cd = C;
Dd = D;

disp('Discrete-time Ad:');
disp(Ad);

disp('Discrete-time Bd:');
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
% TEMPORARY STABILIZING FEEDBACK
% Used only to safely collect rich data
% =========================================
Ktemp = [-1.0  -2.0  18.0  3.5];

Acl_temp = Ad - Bd*Ktemp;

disp('Closed-loop eigenvalues during data collection:');
disp(eig(Acl_temp));

%% =========================================
% PERSISTENTLY EXCITING INPUT
% Piecewise-random signal
% =========================================
u_amp = 0.15;
hold_steps = 8;

num_blocks = ceil(T/hold_steps);
u_blocks = u_amp * sign(randn(1,num_blocks));
u_exc = repelem(u_blocks, hold_steps);
u_exc = u_exc(1:T);

%% =========================================
% PREALLOCATE
% =========================================
x = zeros(n,T+1);
y = zeros(p,T);
u = zeros(m,T);

x(:,1) = x0;

%% =========================================
% SIMULATION
% u(k) = -Ktemp*x(k) + excitation
% =========================================
for k = 1:T
    u(:,k) = -Ktemp*x(:,k) + u_exc(k);
    y(:,k) = Cd*x(:,k) + Dd*u(:,k);
    x(:,k+1) = Ad*x(:,k) + Bd*u(:,k);
end

%% =========================================
% BUILD DATA MATRICES
% =========================================
X0 = x(:,1:T);        % x(0) ... x(T-1)
X1 = x(:,2:T+1);      % x(1) ... x(T)
U0 = u(:,1:T);        % u(0) ... u(T-1)
Y0 = y(:,1:T);        % y(0) ... y(T-1)

Y1 = zeros(p,T);
Y1(:,1:T-1) = y(:,2:T);
Y1(:,T) = Cd*x(:,T+1) + Dd*u(:,T);

%% =========================================
% CHECK DATA RANK
% =========================================
rank_data = rank([U0; X0]);
required_rank = n + m;

fprintf('rank([U0; X0]) = %d\n', rank_data);
fprintf('required rank  = %d\n', required_rank);

if rank_data < required_rank
    warning('Data is NOT rich enough.');
else
    disp('Data matrix has full rank. Good for data-driven control.');
end

%% =========================================
% SAVE EVERYTHING
% =========================================
save('ddc_data_self_balancing.mat', ...
    'A','B','C','D', ...
    'Ad','Bd','Cd','Dd', ...
    'X0','X1','U0','Y0','Y1', ...
    'x','y','u','Ts', ...
    'Ktemp','u_exc', ...
    'M','m_body','l','I','b','g','a','c','d','Delta');

disp('Saved data to ddc_data_self_balancing.mat');