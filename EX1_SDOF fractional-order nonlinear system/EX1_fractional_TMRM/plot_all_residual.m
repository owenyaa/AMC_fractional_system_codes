clear; clc; close all;
load N_5_x_1_nonliear_stiffness_matrix_0_epsilon_fu0.8_alpha_0.5_omega_2_1_parameter_constant_1.1.mat;
time_step = 0.01;
time_data = 0:time_step:10;
% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals over time
subplot(3,1,1);
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

clear; clc; 
load N_10_x_1_nonliear_stiffness_matrix_0_epsilon_fu0.8_alpha_0.5_omega_2_1_parameter_constant_1.1.mat;
time_step = 0.01;
time_data = 0:time_step:10;
% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals over time
subplot(3,1,2);
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

clear; clc; 
load N_15_x_1_nonliear_stiffness_matrix_0_epsilon_fu0.8_alpha_0.5_omega_2_1_parameter_constant_1.1.mat;
time_step = 0.01;
time_data = 0:time_step:10;
% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals over time
subplot(3,1,3);
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);





clear; clc; %close all;
load N_5_x_1_nonliear_stiffness_matrix_1_epsilon_fu0.8_alpha_0.5_omega_2_1_parameter_constant_1.1.mat;
time_step = 0.01;
time_data = 0:time_step:10;
% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

figure;
% Plot the residuals over time
subplot(3,1,1);
plot(time_data, final_residual(:, 1), 'b-', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

clear; clc; 
load N_10_x_1_nonliear_stiffness_matrix_1_epsilon_fu0.8_alpha_0.5_omega_2_1_parameter_constant_1.1.mat;
time_step = 0.01;
time_data = 0:time_step:10;
% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals over time
subplot(3,1,2);
plot(time_data, final_residual(:, 1), 'b-', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

clear; clc; 
load N_15_x_1_nonliear_stiffness_matrix_1_epsilon_fu0.8_alpha_0.5_omega_2_1_parameter_constant_1.1.mat;
time_step = 0.01;
time_data = 0:time_step:10;
% Calculate final residual
final_residual = calculate_residual(harmonic_coefficients);

% Plot the residuals over time
subplot(3,1,3);
plot(time_data, final_residual(:, 1), 'b-', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
























