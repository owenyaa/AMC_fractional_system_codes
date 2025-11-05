%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 22:13
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 2 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%          \begin{equation}
% 		   \label{eq3.1}
% 		        \ddot{x}+\epsilon(1-x^2)D^\alpha\dot{x}+kx^3=0
% 	       \end{equation}
%          The structure of the program is as follows:%
%          The basis function in this program is
%          \begin{equation}
%          	\label{eq4.17}
%          	x_1\approx x^N_1=\sum_{k=1}^{N}\left[b_k\cos\left((2k-1)\omega t\right)+c_k\sin\left((2k-1)\omega t\right)\right]
%          \end{equation}
%          parameter_a :the coefficients%
%          N_harm: the reserved order%
%          Tdata: the duration%
%          Etol:the convergence error%
%          SSS: response sensitivity matrix%
%          dR:residual vector%
%          lambda_inverse:regularization parameter%
clear; clc; close all;

% Start timing the execution
tic;

% Declare global variables
global total_time step_size time_data harmonic_coefficients initial_frequency num_degrees_freedom num_harmonics
global alpha1 alpha2 mass_matrix damping_matrix stiffness_matrix
global mu k alpha beta epsilon
% Initialize degrees of freedom and number of harmonics
num_degrees_freedom = 2;
num_harmonics = 25;
alpha1=1.5; alpha2=0.9; mu=0.8; k=3; alpha=0.8; beta=0; epsilon=1.2;

mass_matrix=[1,0;0,1];damping_matrix=[-mu 0;-k alpha+k];stiffness_matrix=[1+k -k;0 -1];

initial_frequency = 1.3993731;
% initial_frequency = 0.79;

% Define time parameters
total_time = 2 * pi / initial_frequency;
step_size = 2 * pi / (initial_frequency * 500);
time_data = (0:step_size:total_time);

% Initialize harmonic coefficients matrix
harmonic_coefficients = zeros(num_harmonics + 2, 2 * num_degrees_freedom);
harmonic_coefficients(1, 1) = initial_frequency;

% %x1
% harmonic_coefficients(2, :) = [2.524160005360648e-01    0 3.365630642647408e-01     0]; 
% harmonic_coefficients(3, :) = [-5.198583761705369e-03    -2.905762232650881e-01 7.370734049693495e-02    -1.977397489028559e-01]; % Initial values for the first-order harmonics
% harmonic_coefficients(4, :) = [-2.951697149698812e-03    -9.188683918653191e-03 6.487597579199397e-03    -7.240506183569443e-03]; % Initial values for the first-order harmonics
% harmonic_coefficients(5, :) = [1.411098517448326e-03    -2.004486632877274e-04 6.240699861789847e-04     9.877870057471831e-04]; % Initial values for the first-order harmonics
% harmonic_coefficients(6, :) = [-5.695513797987136e-05    -2.994238851235780e-05 -3.361085423885410e-05    -3.719158493101851e-05]; % Initial values for the first-order harmonics

load matlab_para.mat;

% Store initial harmonic coefficients
ini_harmonic_coefficients = harmonic_coefficients;

% Define the number of iterations based on time data length
iteration = length(time_data);

% Function to check if coefficients are within acceptable range
coefficients_check = @(harmonic_coefficients)(abs(harmonic_coefficients) < 5);

% Algorithm parameters for trust-region algorithm
gamma_trust = 1.414;
rho_trust = 0.5;

% Record values of parameters during iteration
harmonic_coefficients_record = harmonic_coefficients;
trust_region_record = [];

% Maximum number for response sensitivity iteration
max_iterations_rs = 1000;

% Maximum number for trust region iteration
max_iterations_tr = 20;

% Define the relative error tolerance for convergence of the algorithm
Etol = 1e-10;
index_bifurcation=1;
% Main loop for response sensitivity iteration
% for alpha1=1:-0.001:1
    for iteration_index = 1:max_iterations_rs
        % Reload frequency and update time parameters
        initial_frequency = harmonic_coefficients(1, 1);
        total_time = 2 * pi / initial_frequency;
        step_size = 2 * pi / (initial_frequency * 500);
        time_data = (0:step_size:total_time);

        % Compute residual identification
        residual_identification = calculate_residual_even_without_velocity(harmonic_coefficients);
        % Extract response sensitivity matrix

        sensitivity_matrix = reshape(residual_identification(:,num_degrees_freedom+1:2*num_degrees_freedom), num_degrees_freedom*length(time_data), 1);

        for i = 1:2*num_harmonics*num_degrees_freedom+2
            sensitivity_matrix = [sensitivity_matrix, reshape(residual_identification(:,num_degrees_freedom*(i+1)+1:num_degrees_freedom*(i+2)), num_degrees_freedom*length(time_data), 1)];
        end

        % Calculate neg_residual (negative of the first column of residual)
        neg_residual = -reshape(residual_identification(:,1:num_degrees_freedom), num_degrees_freedom*length(time_data), 1);

        % Compute SVD of sensitivity matrix
        [U, s, V] = csvd(sensitivity_matrix);

        % Calculate the optimal lambda using L-curve method
        lambda_inverse = l_curve(U, s, neg_residual);

        % Store a copy of the harmonic coefficients matrix for later use
        temp_harmonic_coefficients = harmonic_coefficients;

        % Trust-region algorithm loop
        for trust_index = 1:max_iterations_tr
            % Calculate real coefficient updates
            % Compute real_coeff_updates (optimal parameter increment)
            real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_inverse);
            ini_real_coeff_updates=real_coeff_updates;            
            % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
            temp_real_frequency = zeros(1, 2 * num_degrees_freedom);
            temp_real_frequency(1, 1) = real_coeff_updates(1, 1);
            real_coeff_updates(3, 1) = 0;

            coeff_updates = reshape(real_coeff_updates(2:end-2,:), 2, num_degrees_freedom*num_harmonics);
            coeff_updates = coeff_updates';

            % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
            sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);

            for dof_index = 1:num_degrees_freedom-1
                sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index*num_harmonics+1:(dof_index+1)*num_harmonics, 1:2)];
            end
            %sensitivity_parameter_updates = [temp_real_frequency; sensitivity_parameter_updates];

            sensitivity_const = zeros(1, 2 * num_degrees_freedom);
            sensitivity_const(1,1) = real_coeff_updates(end-1,:);
            sensitivity_const(1,3) = real_coeff_updates(end,:);
            sensitivity_parameter_updates = [temp_real_frequency;sensitivity_const; sensitivity_parameter_updates];

            % Check if updated parameters are within the acceptable range
            if ~coefficients_check(temp_harmonic_coefficients + sensitivity_parameter_updates)
                lambda_inverse = lambda_inverse * gamma_trust;
                continue;
            end

            % Calculate real coefficient updates
            real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_inverse);
            ini_real_coeff_updates=real_coeff_updates;

            % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
            temp_real_frequency = zeros(1, 2 * num_degrees_freedom);
            temp_real_frequency(1, 1) = real_coeff_updates(1, 1);
            real_coeff_updates(3, 1) = 0;

            coeff_updates = reshape(real_coeff_updates(2:end-2,:), 2, num_degrees_freedom*num_harmonics);
            coeff_updates = coeff_updates';

            % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
            sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);

            for dof_index = 1:num_degrees_freedom-1
                sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index*num_harmonics+1:(dof_index+1)*num_harmonics, 1:2)];
            end
            %sensitivity_parameter_updates = [temp_real_frequency; sensitivity_parameter_updates];

            sensitivity_const = zeros(1, 2 * num_degrees_freedom);
            sensitivity_const(1,1) = real_coeff_updates(end-1,:);
            sensitivity_const(1,3) = real_coeff_updates(end,:);
            sensitivity_parameter_updates = [temp_real_frequency;sensitivity_const; sensitivity_parameter_updates];

            % Update harmonic coefficients and recalculate residuals
            harmonic_coefficients = temp_harmonic_coefficients + sensitivity_parameter_updates;
            initial_frequency = harmonic_coefficients(1, 1);
            total_time = 2 * pi / initial_frequency;
            step_size = 2 * pi / (initial_frequency * 500);
            time_data = (0:step_size:total_time);
            residual_updates = calculate_residual_even_without_velocity(harmonic_coefficients);

            % Calculate agreement indicator and check for convergence
            neg_residual_temp = -reshape(residual_updates(:,1:num_degrees_freedom), num_degrees_freedom*length(time_data), 1);
            %-reshape(residual_identification(:,1:(num_degrees_freedom+1)), (num_degrees_freedom+1)*length(time_data), 1);
            L_neg_residual = sensitivity_matrix * ini_real_coeff_updates - neg_residual;
            rho_agreement = (neg_residual' * neg_residual - neg_residual_temp' * neg_residual_temp) / (neg_residual' * neg_residual - L_neg_residual' * L_neg_residual);
            if rho_agreement >= rho_trust
                break;
            end
            lambda_inverse = lambda_inverse * gamma_trust;
        end

        % Check for overall convergence
        tol_convergence = norm(real_coeff_updates) / norm(harmonic_coefficients);
        harmonic_coefficients_record = [harmonic_coefficients_record, harmonic_coefficients];
        trust_region_record = [trust_region_record; lambda_inverse];
        harmonic_coefficients
        if tol_convergence <= Etol
            break;
        end
        all_coefficients(iteration_index).harmonic_coefficients = harmonic_coefficients;
        iteration_index
    end
%     every_harmonic_coefficients(index_bifurcation).alpha1=alpha1;
%     every_harmonic_coefficients(index_bifurcation).harmonic_coefficients=harmonic_coefficients;
%     index_residual=norm(residual_updates(:, 1));
%     every_harmonic_coefficients(index_bifurcation).index_residual=index_residual;
%     index_bifurcation=index_bifurcation+1;
% end

% End timing the execution
toc;

time_step = 0.01;
time_data = 0:time_step:200;
% Calculate final residual
final_residual = calculate_residual_even_without_velocity(harmonic_coefficients);

% Plot the residuals over time
figure;
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Define new time data for response analysis
% time_data = 0:0.04:200;

initial_frequency = harmonic_coefficients(1, 1);
harmonic_parameters = harmonic_coefficients(3:end, :);

% Initialize arrays for displacement, velocity, and acceleration
displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

% Compute displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        harmonic_term = harm_index  * initial_frequency * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
        velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);
        acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
    end
    displacement(dof_index, :) = displacement(dof_index, :) + harmonic_coefficients(2, 2*dof_index-1);
end

% Plot displacement response over time
figure;
plot(time_data, displacement(1, :), 'k-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Plot displacement response over time
figure;
plot(time_data, velocity(1, :), 'k-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);



% Plot phase diagrams for displacement and velocity
figure;
plot(displacement(1, :), velocity(1, :), 'k-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Store initial harmonic coefficients for later reference
ini_harmonic_coefficients

