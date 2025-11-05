clear; clc; close all;
% load bifurcation_epsilon_0.1_omega_1_parameter_constant_1.1_nonliear_stiffness_matrix_1_alpha_0.01_0.01_1.mat;
load bifurcation_epsilon_fu0.8_omega_1_parameter_constant_1.1_nonliear_stiffness_matrix_1_alpha_1_0.01_0.01.mat;
num_degrees_freedom = 1;
num_harmonics = 10;
index_bifurcation=1;
for alpha=1:-0.01:0.01
    %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
    harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
    index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;

    initial_frequency = harmonic_coefficients(1, 1);
    total_time = 2 * pi / initial_frequency;
    step_size = 2 * pi / (initial_frequency * 500);
    time_data = (0:step_size:total_time);
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
    max_displacement(index_bifurcation)=max(displacement);
    alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation).alpha;
    frequency_bif(index_bifurcation)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
end
% Plot displacement response over time
figure;
plot(alpha_bif, max_displacement, 'k-', 'LineWidth', 1);
legend_handle = legend('$$\alpha$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

figure;
plot(alpha_bif, frequency_bif, 'k-', 'LineWidth', 1);
legend_handle = legend('$$\omega$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);