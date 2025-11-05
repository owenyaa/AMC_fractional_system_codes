clear; clc; close all;

load N_10_x_1_nonliear_stiffness_matrix_0_epsilon_fu0.8_alpha_0.5_omega_2_1_parameter_constant_1.1.mat;

time_step = 0.01;
time_data = 0:time_step:20;
% Calculate final residual

initial_frequency = harmonic_coefficients(1, 1);
harmonic_parameters = harmonic_coefficients(2:end, :);

% Initialize arrays for displacement, velocity, and acceleration
displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

% Compute displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        frequency_term = (2 * harm_index - 1) * initial_frequency * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
        velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(frequency_term) + initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index) * cos(frequency_term);
        acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
    end
end

% Plot displacement response over time
% subplot(2,1,1);
% plot(time_data, displacement(1, :), 'r-', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% 
% % Plot phase diagrams for displacement and velocity
figure;
plot(displacement(1, :), velocity(1, :), 'r-', 'LineWidth', 1.5);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


clear; clc; 
load N_10_x_1_nonliear_stiffness_matrix_1_epsilon_fu0.8_alpha_0.5_omega_2_1_parameter_constant_1.1.mat;

time_step = 0.01;
time_data = 0:time_step:20;
% Calculate final residual

initial_frequency = harmonic_coefficients(1, 1);
harmonic_parameters = harmonic_coefficients(2:end, :);

% Initialize arrays for displacement, velocity, and acceleration
displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

% Compute displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        frequency_term = (2 * harm_index - 1) * initial_frequency * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
        velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(frequency_term) + initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index) * cos(frequency_term);
        acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(frequency_term) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(frequency_term);
    end
end

% Plot displacement response over time
% subplot(2,1,2);
% plot(time_data, displacement(1, :), 'b-', 'LineWidth', 1.5);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% 
% Plot phase diagrams for displacement and velocity
figure;
plot(displacement(1, :), velocity(1, :), 'b-', 'LineWidth', 1.5);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
















