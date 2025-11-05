clear; clc; %close all;
load bifurcation_epsilon_0.1_alpha_0.5_parameter_constant_1.1_nonliear_stiffness_matrix_1_omega_2_1_0.01_fu1.mat;
num_degrees_freedom = 1;
num_harmonics = 10;
index_bifurcation=189;
% for omega_2=1:-0.01:-0.88
for omega_2=-0.88:-0.01:-0.88
    %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
    harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
    index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;

    initial_frequency = harmonic_coefficients(1, 1);
    harmonic_parameters = harmonic_coefficients(2:end, :);

    initial_frequency = harmonic_coefficients(1, 1);
    total_time = 2 * pi / initial_frequency;
    step_size = 0.01;
    time_data = (0:step_size:total_time);

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
    omega_2_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation).omega_2;
    frequency_bif(index_bifurcation)=initial_frequency;
    index_bifurcation=index_bifurcation+1;

    hold on;
    plot(displacement(1, :), velocity(1, :), 'r--', 'LineWidth', 1);
    legend_handle = legend('$$x_1$$');
    set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


end
% Plot displacement response over time
figure;
plot(omega_2_bif, max_displacement, 'r--', 'LineWidth', 1);
legend_handle = legend('$$\dot{x}_{123}$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

figure;
plot(omega_2_bif, frequency_bif, 'k-', 'LineWidth', 1);
legend_handle = legend('$$\omega$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);