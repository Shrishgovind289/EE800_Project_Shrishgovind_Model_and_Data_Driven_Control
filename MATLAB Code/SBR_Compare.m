clc;
clear;
close all;

%% =========================================
% LOAD CONTROLLER RESULTS
% =========================================

load('SBR_DataDriven_Controller_Disturbance_Result.mat');      
% Expected variables:
% K_data, Acl_dd, Ts

load('SBR_Dissipative_DataDriven_Controller_Disturbance_Result.mat');     
% Expected variables:
% K_diss, Acl_diss or Acl_diss_dd

%% =========================================
% BASIC SETTINGS
% =========================================

Ntest = 1500;              % simulation samples
t = 0:Ts:Ntest*Ts;
t_input = 0:Ts:(Ntest-1)*Ts;

n = 4;
m = 1;

dist_time = 10;            % disturbance at 10 seconds
dist_idx = find(t_input >= dist_time, 1, 'first');

%% =========================================
% INITIAL CONDITIONS
% =========================================

x0 = [0;
      0;
      deg2rad(45);
      0];

x_data = zeros(n, Ntest+1);
x_diss = zeros(n, Ntest+1);

u_data = zeros(m, Ntest);
u_diss = zeros(m, Ntest);

x_data(:,1) = x0;
x_diss(:,1) = x0;

%% =========================================
% DISTURBANCE SETTINGS
% =========================================
% Disturbance added to angle and angular velocity

angle_disturbance = deg2rad(-20);       % angle kick
angvel_disturbance = deg2rad(-80);      % angular velocity kick

disturbance = zeros(n, Ntest);

disturbance(3, dist_idx) = angle_disturbance;
disturbance(4, dist_idx) = angvel_disturbance;

%% =========================================
% SELECT CLOSED-LOOP MATRICES
% =========================================
% Pure data-driven closed-loop matrix

Acl_data = Acl_dd;

% Dissipative closed-loop matrix
% Use whichever variable exists in your saved dissipative result file

if exist('Acl_diss_dd','var')
    Acl_diss = Acl_diss_dd;
elseif exist('Acl_diss','var')
    Acl_diss = Acl_diss;
elseif exist('Acl','var')
    Acl_diss = Acl;
else
    error('No dissipative closed-loop matrix found. Save Acl_diss or Acl_diss_dd in your dissipative controller file.');
end

%% =========================================
% SIMULATE BOTH CONTROLLERS
% =========================================
% Important:
% If your convention is u = Kx, keep as shown.
% If your convention is u = -Kx, change both lines to -K*x.

for k = 1:Ntest

    % -------------------------
    % Pure data-driven controller
    % -------------------------
    u_data(:,k) = K_data * x_data(:,k);
    x_data(:,k+1) = Acl_data * x_data(:,k) + disturbance(:,k);

    % -------------------------
    % Dissipative data-driven controller
    % -------------------------
    u_diss(:,k) = K_diss * x_diss(:,k);
    x_diss(:,k+1) = Acl_diss * x_diss(:,k) + disturbance(:,k);

end

%% =========================================
% CONVERT STATES FOR PLOTTING
% =========================================

phi_data_deg = rad2deg(x_data(3,:));
phi_diss_deg = rad2deg(x_diss(3,:));

phidot_data_deg = rad2deg(x_data(4,:));
phidot_diss_deg = rad2deg(x_diss(4,:));

%% =========================================
% CHECK IF DISSIPATIVE PART IS ACTUALLY MOVING
% =========================================

fprintf('\n=========== SIMULATION CHECK ===========\n');
fprintf('Max |x_data| = %.6f\n', max(abs(x_data(:))));
fprintf('Max |x_diss| = %.6f\n', max(abs(x_diss(:))));
fprintf('Max |u_data| = %.6f\n', max(abs(u_data(:))));
fprintf('Max |u_diss| = %.6f\n', max(abs(u_diss(:))));
fprintf('Norm difference between u_data and u_diss = %.6f\n', norm(u_data - u_diss));

if norm(x_diss(:)) < 1e-8
    warning('x_diss is almost zero. Dissipative simulation is still not updating correctly.');
end

if norm(u_data - u_diss) < 1e-8
    warning('u_data and u_diss are almost identical. Check if the same gain/state is being used accidentally.');
end

%% =========================================
% DISSIPATIVITY VISUAL METRIC
% z(k) = Cz*x(k)
% =========================================

Cz = diag([0.2 0.05 1 0.2]);

z_data = Cz * x_data;
z_diss = Cz * x_diss;

z_data_norm = vecnorm(z_data, 2, 1);
z_diss_norm = vecnorm(z_diss, 2, 1);

E_data = sum(z_data_norm(dist_idx:end).^2);
E_diss = sum(z_diss_norm(dist_idx:end).^2);

fprintf('\n=========== DISSIPATIVITY METRIC ==========\n');
fprintf('Post-disturbance energy Data-Driven       = %.6f\n', E_data);
fprintf('Post-disturbance energy Dissipative       = %.6f\n', E_diss);

%% =========================================
% PLOT SETTINGS
% =========================================

tick_interval = 0.5;
line_width = 2;
dist_line_width = 1.5;

%% =========================================
% 1. BODY ANGLE COMPARISON
% =========================================

figure('Name','Body Angle Comparison','Color','w');
plot(t, phi_data_deg, 'LineWidth', line_width); hold on;
plot(t, phi_diss_deg, 'LineWidth', line_width);
xline(dist_time, '--k', 'Disturbance', ...
    'LineWidth', dist_line_width, ...
    'LabelOrientation','aligned');

grid on;
xlabel('Time (s)');
ylabel('\phi (deg)', 'Interpreter','tex');
title('Body Angle Comparison Under Disturbance');
legend('Data-Driven','Dissipative','Location','best');
xlim([0 max(t)]);
xticks(0:tick_interval:max(t));

%% =========================================
% 2. ANGULAR VELOCITY COMPARISON
% =========================================

figure('Name','Angular Velocity Comparison','Color','w');
plot(t, phidot_data_deg, 'LineWidth', line_width); hold on;
plot(t, phidot_diss_deg, 'LineWidth', line_width);
xline(dist_time, '--k', 'Disturbance', ...
    'LineWidth', dist_line_width, ...
    'LabelOrientation','aligned');

grid on;
xlabel('Time (s)');
ylabel('$\dot{\phi}$ (deg/s)', 'Interpreter','latex');
title('Angular Velocity Comparison Under Disturbance');
legend('Data-Driven','Dissipative','Location','best');
xlim([0 max(t)]);
xticks(0:tick_interval:max(t));

%% =========================================
% 3. CART POSITION COMPARISON
% =========================================

figure('Name','Robot Position Comparison','Color','w');
plot(t, x_data(1,:), 'LineWidth', line_width); hold on;
plot(t, x_diss(1,:), 'LineWidth', line_width);
xline(dist_time, '--k', 'Disturbance', ...
    'LineWidth', dist_line_width, ...
    'LabelOrientation','aligned');

grid on;
xlabel('Time (s)');
ylabel('x (m)');
title('Robot Position Comparison Under Disturbance');
legend('Data-Driven','Dissipative','Location','best');
xlim([0 max(t)]);
xticks(0:tick_interval:max(t));

%% =========================================
% 4. CART VELOCITY COMPARISON
% =========================================

figure('Name','Robot Velocity Comparison','Color','w');
plot(t, x_data(2,:), 'LineWidth', line_width); hold on;
plot(t, x_diss(2,:), 'LineWidth', line_width);
xline(dist_time, '--k', 'Disturbance', ...
    'LineWidth', dist_line_width, ...
    'LabelOrientation','aligned');

grid on;
xlabel('Time (s)');
ylabel('$\dot{x}$ (m/s)', 'Interpreter','latex');
title('Robot Velocity Comparison Under Disturbance');
legend('Data-Driven','Dissipative','Location','best');
xlim([0 max(t)]);
xticks(0:tick_interval:max(t));

%% =========================================
% 5. CONTROL INPUT COMPARISON
% =========================================

figure('Name','Control Input Comparison','Color','w');
plot(t_input, u_data, 'LineWidth', line_width); hold on;
plot(t_input, u_diss, 'LineWidth', line_width);
xline(dist_time, '--k', 'Disturbance', ...
    'LineWidth', dist_line_width, ...
    'LabelOrientation','aligned');

grid on;
xlabel('Time (s)');
ylabel('Control input u');
title('Control Input Comparison Under Disturbance');
legend('Data-Driven','Dissipative','Location','best');
xlim([0 max(t_input)]);
xticks(0:tick_interval:max(t_input));

%% =========================================
% 6. DISSIPATIVITY PERFORMANCE OUTPUT NORM
% =========================================

figure('Name','Dissipativity Performance Output Norm','Color','w');
plot(t, z_data_norm, 'LineWidth', line_width); hold on;
plot(t, z_diss_norm, 'LineWidth', line_width);
xline(dist_time, '--k', 'Disturbance', ...
    'LineWidth', dist_line_width, ...
    'LabelOrientation','aligned');

grid on;
xlabel('Time (s)');
ylabel('||z(k)||');
title('Dissipativity Performance Output Norm');
legend('Data-Driven','Dissipative','Location','best');
xlim([0 max(t)]);
xticks(0:tick_interval:max(t));

%% =========================================
% 7. POST-DISTURBANCE CUMULATIVE ENERGY
% =========================================

t_energy = t(dist_idx:end);

energy_data = cumsum(z_data_norm(dist_idx:end).^2);
energy_diss = cumsum(z_diss_norm(dist_idx:end).^2);

figure('Name','Post-Disturbance Dissipativity Energy','Color','w');
plot(t_energy, energy_data, 'LineWidth', line_width); hold on;
plot(t_energy, energy_diss, 'LineWidth', line_width);

grid on;
xlabel('Time (s)');
ylabel('Cumulative \Sigma ||z(k)||^2');
title('Post-Disturbance Dissipativity Energy Comparison');
legend('Data-Driven','Dissipative','Location','best');
xlim([dist_time max(t)]);
xticks(dist_time:tick_interval:max(t));

%% =========================================
% 8. BODY ANGLE + DISSIPATIVITY METRIC TOGETHER
% =========================================

figure('Name','Body Angle with Dissipativity Metric','Color','w');

yyaxis left
h1 = plot(t, phi_data_deg, 'LineWidth', line_width); hold on;
h2 = plot(t, phi_diss_deg, 'LineWidth', line_width);
ylabel('\phi (deg)', 'Interpreter','tex');

yyaxis right
h3 = plot(t, z_data_norm, ':', 'LineWidth', line_width);
h4 = plot(t, z_diss_norm, '--', 'LineWidth', line_width);
ylabel('Dissipativity metric ||z(k)||');

xline(dist_time, '--k', 'Disturbance', ...
    'LineWidth', dist_line_width, ...
    'LabelOrientation','aligned');

grid on;
xlabel('Time (s)');
title('Body Angle Response with Dissipativity Metric');

legend([h1 h2 h3 h4], ...
    '\phi Data-Driven', ...
    '\phi Dissipative', ...
    '||z|| Data-Driven', ...
    '||z|| Dissipative', ...
    'Location','best');

xlim([0 max(t)]);
xticks(0:tick_interval:max(t));

txt = sprintf(['Post-disturbance energy:\n', ...
               'Data-Driven = %.4f\n', ...
               'Dissipative = %.4f'], E_data, E_diss);

annotation('textbox',[0.58 0.18 0.28 0.16], ...
    'String',txt, ...
    'FitBoxToText','on', ...
    'BackgroundColor','white');

%% =========================================
% SAVE RESULT
% =========================================

save('SBR_Controller_Comparison_Result_FIXED.mat', ...
    'x_data','x_diss', ...
    'u_data','u_diss', ...
    'z_data','z_diss', ...
    'z_data_norm','z_diss_norm', ...
    'energy_data','energy_diss', ...
    'E_data','E_diss', ...
    't','t_input','Ts','dist_time');

disp('Fixed comparison result saved as SBR_Controller_Comparison_Result_FIXED.mat');