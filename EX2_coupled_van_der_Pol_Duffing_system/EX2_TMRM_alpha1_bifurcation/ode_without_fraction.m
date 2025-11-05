clear; clc; close all;  % Clear workspace, command window, and close all figures
tic  % Start timing
global mass_matrix damping_matrix stiffness_matrix eta beta omega_cos N_2 omega_R1
% 定义参数
A = 12;
xi = 0.003;
lambda = 0.31831;
St = 0.2;
epsilon_omega = 0.3;
omega_R1 = 0;
omega_R2 = 1;
alpha = 0.5;
eta = 0.8;
beta = 1;
omega_cos = 0.2;
delta = 0.01194;
Omegas = 0.6;

m_21 = -A;
c_11 = 2 * xi + lambda * Omegas / (pi * St);
c_22 = -epsilon_omega * Omegas;
d_11 = omega_R1^2;
k_11_1 = omega_R2^2;
k_12 = -delta * Omegas^2 / (pi^2 * St^2);
k_22 = Omegas^2;
N_2 = 0.75 * epsilon_omega * Omegas;

mass_matrix = [1, 0; m_21, 1];
damping_matrix = [c_11, 0; 0, c_22];
stiffness_matrix = [k_11_1, k_12; 0, k_22];

% 时间参数
time_step = 0.01;
time_data = 0:time_step:2000;
y0 = [1; 2; 0; 0]; % 初始条件 [x1, x2, v1, v2]

% 使用 ode45 求解微分方程
options=odeset('RelTol',1e-8,'AbsTol',1e-8);
% [t,num]=ode45('odehomo',tt,odex,options);

[t, y] = ode45(@(t, y) odefun(t, y), time_data, y0 ,options);
displacement = y(:, 1:2);
velocity = y(:, 3:4);

% Plot displacement response over time
figure;
plot(time_data(180000:end), displacement(180000:end, 1), 'r-', 'LineWidth', 1);
hold on;
plot(time_data(180000:end), displacement(180000:end, 2), 'b-', 'LineWidth', 1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
legend('$$x_1$$', '$$x_2$$', 'Interpreter', 'latex', 'FontSize', 15);

% Phase plot for both degrees of freedom
figure;
plot(displacement(180000:end, 1), velocity(180000:end, 1), 'r-', 'LineWidth', 1);
hold on;
plot(displacement(180000:end, 2), velocity(180000:end, 2), 'b-', 'LineWidth', 1);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
legend('$$x_1$$', '$$x_2$$', 'Interpreter', 'latex', 'FontSize', 15);

function dydt = odefun(t, y)
    % 状态变量
    global mass_matrix damping_matrix stiffness_matrix eta beta omega_cos N_2 omega_R1
    x1 = y(1);
    x2 = y(2);
    v1 = y(3);
    v2 = y(4);

    % 时变刚度项
    k_11 = stiffness_matrix(1, 1) + eta * beta * cos(omega_cos * t)+omega_R1^2;
    stiffness_matrix1=stiffness_matrix;
    % 定义刚度矩阵
    stiffness_matrix1(1, 1)= k_11;

    % 状态空间表示
    dydt = zeros(4, 1);
    dydt(1) = v1;
    dydt(2) = v2;

    % 计算加速度
    A = mass_matrix \ (-stiffness_matrix1 * [x1; x2] - damping_matrix * [v1; v2] - [0; N_2 * x2^2 * v2]);
    dydt(3) = A(1);
    dydt(4) = A(2);
end
