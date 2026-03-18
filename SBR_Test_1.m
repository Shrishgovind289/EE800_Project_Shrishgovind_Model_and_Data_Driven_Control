clc;              
clear;            
close all;        

%% ============================================================
% LMI-BASED STATE FEEDBACK CONTROL FOR A SELF-BALANCING ROBOT
%
% State vector:
%   x1 = translational position of the robot base
%   x2 = translational velocity of the robot base
%   x3 = body tilt angle (phi)
%   x4 = body angular velocity
%
% Control law used in this formulation:
%   u = Kx
%
% Closed-loop system:
%   x_dot = (A + B*K)x
%
% Objective:
%   1. Build the linearized state-space model
%   2. Design a stabilizing feedback gain K using an LMI
%   3. Verify closed-loop stability
%   4. Simulate the closed-loop response using ode45
%   5. Plot the state trajectories and control input
%% ============================================================

%% 1) Physical parameters
M = 0.5;      
% Mass of the base / wheel assembly (kg)

m = 0.2;      
% Mass of the robot body (kg)

b = 0.1;      
% Viscous friction coefficient

I = 0.006;    
% Moment of inertia of the body about its center of mass (kg.m^2)

g = 9.81;     
% Gravitational acceleration (m/s^2)

l = 0.3;      
% Distance from wheel axle to center of mass (m)

%% 2) Common denominator
Delta = (M + m)*(I + m*l^2) - (m*l)^2;
% This term appears when solving the coupled translational and rotational
% equations for x_ddot and phi_ddot.

%% 3) Linearized state-space model
% State vector:
%   x = [position; velocity; angle; angular_velocity]

A = [ 0, 1, 0, 0;
      0, -(I + m*l^2)*b/Delta,  (m^2*g*l^2)/Delta, 0;
      0, 0, 0, 1;
      0, -(m*l*b)/Delta,        (m*g*l*(M + m))/Delta, 0 ];
% A is the system matrix

B = [ 0;
      (I + m*l^2)/Delta;
      0;
      (m*l)/Delta ];
% B is the input matrix

%% 4) Open-loop poles
disp('Open-loop poles:');
disp(eig(A));
% These eigenvalues describe the stability of the uncontrolled system

%% 5) Dimensions of the system
n = size(A,1);         
% Number of states

m_in = size(B,2);      
% Number of control inputs

%% 6) YALMIP decision variables
P = sdpvar(n,n,'symmetric');   
% Lyapunov matrix P > 0

Y = sdpvar(m_in,n,'full');     
% Auxiliary variable Y = K*P

%% 7) Small positive scalar for strict inequalities
eps_val = 1e-6;

%% 8) LMI for the control law u = Kx
% Closed-loop system:
%   x_dot = (A + B*K)x
%
% Lyapunov condition:
%   (A + B*K)'P + P(A + B*K) < 0
%
% Using Y = K*P, this becomes:
%   A*P + P*A' + B*Y + Y'*B' < 0

LMI = A*P + P*A' + B*Y + Y'*B';

%% 9) Constraints
Constraints = [];

Constraints = [Constraints, P >= eps_val*eye(n)];
% Enforces P > 0

Constraints = [Constraints, LMI <= -eps_val*eye(n)];
% Enforces the Lyapunov inequality

%% 10) Solver settings
ops = sdpsettings('solver','sdpt3','verbose',1);

%% 11) Solve the LMI feasibility problem
sol = optimize(Constraints, [], ops);

%% 12) Check solver status
if sol.problem ~= 0
    error('LMI problem not solved successfully.');
end

%% 13) Recover numerical values of P and Y
P_val = value(P);
Y_val = value(Y);

%% 14) Recover gain matrix K
% Since Y = K*P, then K = Y*inv(P)
K = Y_val / P_val;

disp('LMI-based feedback gain K = ');
disp(K);

%% 15) Closed-loop system matrix for u = Kx
Acl = A + B*K;

disp('Closed-loop poles:');
disp(eig(Acl));
% If all real parts are negative, the closed-loop system is stable

%% 16) Initial condition
x0 = [0;
      0;
      deg2rad(5);
      0];
% Initial state:
%   position = 0
%   velocity = 0
%   angle = 5 degrees
%   angular velocity = 0

%% 17) Simulation interval
tspan = [0 10];
% Simulate from 0 to 10 seconds

%% 18) Closed-loop dynamics
odefun = @(t, x) Acl*x;
% Since the closed-loop linear system is:
%   x_dot = Acl*x

%% 19) Simulate with ode45
[t, x] = ode45(odefun, tspan, x0);

%% 20) Compute control input over time
u = zeros(length(t),1);

for k = 1:length(t)
    u(k) = K * x(k,:)';
end
% Here the control law is u = Kx, so no negative sign is used

%% 21) Plot all states and control input
figure;

subplot(5,1,1);
plot(t, x(:,1), 'LineWidth', 1.5);
grid on;
ylabel('x (m)');
title('Closed-Loop Response of Self-Balancing Robot');

subplot(5,1,2);
plot(t, x(:,2), 'LineWidth', 1.5);
grid on;
ylabel('x dot (m/s)');

subplot(5,1,3);
plot(t, rad2deg(x(:,3)), 'LineWidth', 1.5);
grid on;
ylabel('\phi (deg)');

subplot(5,1,4);
plot(t, rad2deg(x(:,4)), 'LineWidth', 1.5);
grid on;
ylabel('\phi dot (deg/s)');

subplot(5,1,5);
plot(t, u, 'LineWidth', 1.5);
grid on;
ylabel('u');
xlabel('Time (s)');

%% 22) Separate tilt-angle plot
figure;
plot(t, rad2deg(x(:,3)), 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Tilt Angle \phi (deg)');
title('Tilt Angle Response');

%% 23) Separate position plot
figure;
plot(t, x(:,1), 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Position x (m)');
title('Position Response');

%% 24) Separate control input plot
figure;
plot(t, u, 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Control Input u');
title('Control Effort');