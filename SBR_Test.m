%% ============================================================
% LMI-BASED STATE FEEDBACK CONTROL FOR A SELF-BALANCING ROBOT
%
% State vector:
%   x1 = translational position of the robot base
%   x2 = translational velocity of the robot base
%   x3 = body tilt angle (phi)
%   x4 = body angular velocity
%
% Control law:
%   u = -Kx
%
% Objective:
%   Compute a stabilizing state-feedback gain K using an LMI
%   derived from Lyapunov stability theory.
% ============================================================

%% 1) Physical parameters of the robot
M = 0.5;      % Mass of the base / wheel assembly (kg)
m = 0.2;      % Mass of the robot body (kg)
b = 0.1;      % Viscous friction coefficient (N.s/m)
I = 0.006;    % Moment of inertia of the body about its center of mass (kg.m^2)
g = 9.81;     % Gravitational acceleration (m/s^2)
l = 0.3;      % Distance from wheel axle to center of mass of the body (m)

%% 2) Compute the common denominator used in the linearized equations
% This term appears when solving the coupled translational and rotational
% equations of motion for x_ddot and phi_ddot.
Delta = (M + m)*(I + m*l^2) - (m*l)^2;

%% 3) Construct the state-space matrices A and B
% The linearized system is written as:
%
%   x_dot = A*x + B*u
%
% where the state vector is:
%
%   x = [position; velocity; angle; angular_velocity]
%
% Matrix A describes the internal system dynamics.
A = [ 0, 1, 0, 0;
      0, -(I + m*l^2)*b/Delta,  (m^2*g*l^2)/Delta, 0;
      0, 0, 0, 1;
      0, -(m*l*b)/Delta,        (m*g*l*(M + m))/Delta, 0 ];

% Matrix B describes how the control input u affects the system states.
B = [ 0;
      (I + m*l^2)/Delta;
      0;
      (m*l)/Delta ];

%% 4) Check the open-loop poles of the system
% The eigenvalues of A indicate whether the uncontrolled system is stable.
% For a self-balancing robot linearized around the upright position,
% the open-loop system is typically unstable.
disp('Open-loop poles:');
disp(eig(A));

%% 5) Define dimensions of the state-space model
% n = number of states
% m_in = number of control inputs
n = size(A,1);
m_in = size(B,2);

%% 6) Define LMI decision variables using YALMIP
% P is the Lyapunov matrix and must be symmetric positive definite.
P = sdpvar(n,n,'symmetric');

% Y is an auxiliary matrix introduced using the substitution:
%   Y = K*P
% This converts the nonlinear Lyapunov inequality into a linear matrix inequality.
Y = sdpvar(m_in,n,'full');

%% 7) Define a small positive scalar for numerical strictness
% Instead of enforcing P > 0 and LMI < 0 exactly,
% we use a small epsilon to impose:
%
%   P >= epsilon*I
%   LMI <= -epsilon*I
%
% This improves numerical robustness in the solver.
eps_val = 1e-6;

%% 8) Construct the LMI expression
% Starting from the closed-loop system:
%
%   x_dot = (A - B*K)x
%
% and choosing a Lyapunov function:
%
%   V(x) = x^T P x
%
% the stability condition becomes:
%
%   (A - B*K)'*P + P*(A - B*K) < 0
%
% After substituting Y = K*P, this becomes:
%
%   A*P + P*A' - B*Y - Y'*B' < 0
%
% which is linear in P and Y.
LMI = A*P + P*A' - B*Y - Y'*B';

%% 9) Define the optimization constraints
Constraints = [];

% Enforce P to be positive definite
Constraints = [Constraints, P >= eps_val*eye(n)];

% Enforce the Lyapunov inequality to be negative definite
Constraints = [Constraints, LMI <= -eps_val*eye(n)];

%% 10) Define solver settings
% SDPT3 is selected as the semidefinite programming solver.
% 'verbose',1 tells MATLAB to display solver progress in the command window.
ops = sdpsettings('solver','sdpt3','verbose',1);

%% 11) Solve the LMI feasibility problem
% Since the goal is only to find a stabilizing controller,
% no objective function is required.
% Therefore, the second argument is left as [].
sol = optimize(Constraints, [], ops);

%% 12) Check whether the solver succeeded
% If sol.problem is not zero, the LMI was not solved successfully.
if sol.problem ~= 0
    error('LMI problem not solved successfully. Check YALMIP/solver installation or model formulation.');
end

%% 13) Extract the numerical values of the decision variables
% value(P) and value(Y) convert YALMIP symbolic variables into numeric matrices.
P_val = value(P);
Y_val = value(Y);

%% 14) Recover the state-feedback gain matrix K
% Since Y = K*P, the controller gain is obtained as:
%
%   K = Y*inv(P)
%
% In MATLAB, right division by P is equivalent to multiplying by inv(P).
K = Y_val / P_val;

%% 15) Display the computed feedback gain
disp('LMI-based feedback gain K = ');
disp(K);

%% 16) Form the closed-loop system matrix
% With the state-feedback law u = -Kx, the closed-loop dynamics become:
%
%   x_dot = (A - B*K)x
Acl = A - B*K;

%% 17) Display the closed-loop poles
% If all eigenvalues of Acl have negative real parts,
% the closed-loop system is asymptotically stable.
disp('Closed-loop poles:');
disp(eig(Acl));