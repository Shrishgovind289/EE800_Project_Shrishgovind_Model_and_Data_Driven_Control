clc;
clear;
close all;

%% =========================================
% LOAD SAVED DATA
% =========================================
load('ddc_data_self_balancing.mat');

disp('Loaded variables:');
whos

%% =========================================
% DIMENSIONS
% =========================================
n = size(X0,1);   % number of states
m = size(U0,1);   % number of inputs

fprintf('n = %d\n', n);
fprintf('m = %d\n', m);

%% =========================================
% RECOVER MODEL FROM DATA
% X1 = [B A] * [U0; X0]
% =========================================
AB_est = X1 * pinv([U0; X0]);

B_est = AB_est(:,1:m);
A_est = AB_est(:,m+1:m+n);

disp('Recovered A_est from data:');
disp(A_est);

disp('Recovered B_est from data:');
disp(B_est);

%% =========================================
% COMPARE WITH TRUE DISCRETE MODEL
% =========================================
if exist('Ad','var') && exist('Bd','var')
    fprintf('||A_est - Ad|| = %g\n', norm(A_est - Ad));
    fprintf('||B_est - Bd|| = %g\n', norm(B_est - Bd));
end

%% =========================================
% DESIGN CONTROLLER USING DISCRETE LQR
% =========================================
Q_lqr = diag([10 1 100 1]);   % state weights
R_lqr = 1;                    % input weight

K_lqr = dlqr(A_est, B_est, Q_lqr, R_lqr);

disp('LQR gain from data-recovered model:');
disp(K_lqr);

%% =========================================
% CHECK CLOSED-LOOP STABILITY
% Standard convention: u = -Kx
% =========================================
Acl = A_est - B_est*K_lqr;

disp('Closed-loop matrix Acl = A_est - B_est*K_lqr:');
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
% OPTIONAL CHECK AGAINST TRUE MODEL
% =========================================
if exist('Ad','var') && exist('Bd','var')
    disp('Eigenvalues of Ad - Bd*K_lqr:');
    disp(eig(Ad - Bd*K_lqr));
end

%% =========================================
% SAVE RESULT
% =========================================
save('controller_from_data_result.mat', 'A_est', 'B_est', 'K_lqr', 'Acl', 'cl_eigs');

disp('Saved controller to controller_from_data_result.mat');