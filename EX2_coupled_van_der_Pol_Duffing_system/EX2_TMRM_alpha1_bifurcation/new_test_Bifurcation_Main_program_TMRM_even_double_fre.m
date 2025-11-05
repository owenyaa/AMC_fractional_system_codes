%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 22:13
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : This program combines arc-length increment method and dynamic adjustment of bifurcation parameter step size

clear; %clc; close all;

% Start timing the execution
tic;

% Declare global variables
global total_time step_size time_data harmonic_coefficients initial_frequency num_degrees_freedom num_harmonics
global alpha1 alpha2 mass_matrix damping_matrix stiffness_matrix
global mu k alpha beta epsilon
% Initialize degrees of freedom and number of harmonics
num_degrees_freedom = 2;
num_harmonics = 25;
% alpha1=1.5;
alpha2=0.9;
% alpha1=1; alpha2=1;
mu=0.8; k=3; alpha=0.8; beta=0.05; epsilon=1.2;

mass_matrix=[1,0;0,1]; stiffness_matrix=[1+k -k;0 -1];
damping_matrix=[-mu 0;-k alpha+k];

% Initialize harmonic coefficients matrix
harmonic_coefficients = zeros(num_harmonics + 2, 2 * num_degrees_freedom);
% harmonic_coefficients(1, 1) = initial_frequency;

% Store initial harmonic coefficients
ini_harmonic_coefficients = harmonic_coefficients;
initial_frequency=1;

% Define time parameters
total_time = 2 * pi / initial_frequency;
step_size = 2 * pi / (initial_frequency * 500);
time_data = (0:step_size:total_time);

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

% Load previous results (for testing, use synthetic data)
load matlab_1.503_1.5.mat;

% Initialize I_prev (start with a reasonable default)
iteration_index = 23; % Initial guess for the iteration count

ds = 0.01; % Initial arc length increment
Nd = 20; % Desired number of iterations

index_bifurcation = 5;

% Start loop for bifurcation analysis
if index_bifurcation > 4
    for i = 1:4
        temp_bifurcation_harmonic_coefficients(i).harmonic_coefficients = every_harmonic_coefficients(index_bifurcation - 5 + i).harmonic_coefficients;
    end
    [harmonic_coefficients, ds] = arc_length_new(temp_bifurcation_harmonic_coefficients, Nd, ds, iteration_index);
end

% Define the relative error tolerance for convergence of the algorithm
Etol = 1e-10;

% Main loop for response sensitivity iteration
for alpha1 = 1.499:-0.001:1.48
    for iteration_index = 1:max_iterations_rs
        % Reload frequency and update time parameters
        initial_frequency = harmonic_coefficients(1, 1);
        total_time = 2 * pi / initial_frequency;
        step_size = 2 * pi / (initial_frequency * 500);
        time_data = (0:step_size:total_time);

        % Compute residual identification
        residual_identification = calculate_residual_even(harmonic_coefficients);
        % Extract response sensitivity matrix
        sensitivity_matrix = reshape(residual_identification(:, (num_degrees_freedom + 1) + 1:2 * (num_degrees_freedom + 1)), (num_degrees_freedom + 1) * length(time_data), 1);

        for i = 1:2 * num_harmonics * num_degrees_freedom + 2
            sensitivity_matrix = [sensitivity_matrix, reshape(residual_identification(:, (num_degrees_freedom + 1) * (i + 1) + 1:(num_degrees_freedom + 1) * (i + 2)), (num_degrees_freedom + 1) * length(time_data), 1)];
        end

        % Calculate neg_residual (negative of the first column of residual)
        neg_residual = -reshape(residual_identification(:, 1:(num_degrees_freedom + 1)), (num_degrees_freedom + 1) * length(time_data), 1));

        % Compute SVD of sensitivity matrix
        [U, s, V] = csvd(sensitivity_matrix);

        % Calculate the optimal lambda using L-curve method
        lambda_inverse = l_curve(U, s, neg_residual);

        % Store a copy of the harmonic coefficients matrix for later use
        temp_harmonic_coefficients = harmonic_coefficients;

        % Trust-region algorithm loop
        for trust_index = 1:max_iterations_tr
            % Calculate real coefficient updates
            real_coeff_updates = tikhonov(U, s, V, neg_residual, lambda_inverse);
            ini_real_coeff_updates = real_coeff_updates;

            % Reshape real_coeff_updates to match the harmonic coefficients matrix structure
            temp_real_frequency = zeros(1, 2 * num_degrees_freedom);
            temp_real_frequency(1, 1) = real_coeff_updates(1, 1);
            coeff_updates = reshape(real_coeff_updates(2:end - 2, :), 2, num_degrees_freedom * num_harmonics);
            coeff_updates = coeff_updates';

            % Extract sensitivity_matrix (formerly sensitivity_parameter_da)
            sensitivity_parameter_updates = coeff_updates(1:num_harmonics, 1:2);

            for dof_index = 1:num_degrees_freedom - 1
                sensitivity_parameter_updates = [sensitivity_parameter_updates, coeff_updates(dof_index * num_harmonics + 1:(dof_index + 1) * num_harmonics, 1:2)];
            end

            sensitivity_const = zeros(1, 2 * num_degrees_freedom);
            sensitivity_const(1, 1) = real_coeff_updates(end - 1, :);
            sensitivity_const(1, 3) = real_coeff_updates(end, :);
            sensitivity_parameter_updates = [temp_real_frequency; sensitivity_const; sensitivity_parameter_updates];

            % Check if updated parameters are within the acceptable range
            if ~coefficients_check(temp_harmonic_coefficients + sensitivity_parameter_updates)
                lambda_inverse = lambda_inverse * gamma_trust;
                continue;
            end

            % Update harmonic coefficients and recalculate residuals
            harmonic_coefficients = temp_harmonic_coefficients + sensitivity_parameter_updates;
            initial_frequency = harmonic_coefficients(1, 1);
            total_time = 2 * pi / initial_frequency;
            step_size = 2 * pi / (initial_frequency * 500);
            time_data = (0:step_size:total_time);
            residual_updates = calculate_residual_even(harmonic_coefficients);

            % Calculate agreement indicator and check for convergence
            neg_residual_temp = -reshape(residual_updates(:, 1:(num_degrees_freedom + 1)), (num_degrees_freedom + 1) * length(time_data), 1);
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
    every_harmonic_coefficients(index_bifurcation).alpha1 = alpha1;
    every_harmonic_coefficients(index_bifurcation).harmonic_coefficients = harmonic_coefficients;
    index_residual = norm(residual_updates(:, 1));
    every_harmonic_coefficients(index_bifurcation).index_residual = index_residual;
    index_bifurcation = index_bifurcation + 1;

    % Adjust bifurcation parameter step size dynamically
    % Using arc-length to determine step size
    if index_bifurcation > 4
        for i = 1:4
            temp_bifurcation_harmonic_coefficients(i).harmonic_coefficients = every_harmonic_coefficients(index_bifurcation - 5 + i).harmonic_coefficients;
        end
        [harmonic_coefficients, ds] = arc_length_new(temp_bifurcation_harmonic_coefficients, Nd, ds, iteration_index);
    end
end

% End timing the execution
toc;

% Plot final results and residuals
time_step = 0.01;
time_data = 0:time_step:200;
final_residual = calculate_residual_even(harmonic_coefficients);

% Plot the residuals over time
figure;
plot(time_data, final_residual(:, 1), 'r-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% Calculate displacement, velocity, and acceleration for each degree of freedom
initial_frequency = harmonic_coefficients(1, 1);
harmonic_parameters = harmonic_coefficients(3:end, :);

% Initialize arrays for displacement, velocity, and acceleration
displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

% Compute displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        harmonic_term = harm_index * initial_frequency * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harmonic_term, 2 * dof_index) * sin(harmonic_term);
    end
end
% Plot displacement response over time
figure;
plot(time_data, displacement(1, :), 'k-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


% Arc Length New Method
function [Predictor, ds_new] = arc_length_new(every_a, Nd, ds_old, iteration_index)
    % Arc-length increment method with dynamic step adjustment
    % Inputs:
    % every_a: Structure array containing previous solutions
    % Nd: Desired number of iterations
    % ds_old: Arc length increment from the previous step
    % iteration_index: Iteration count from the main program
    % Outputs:
    % Predictor: Predicted next solution
    % ds_new: Updated arc length increment for the next step

    % Step 1: Calculate arc lengths between known points
    s = zeros(1, 3); % Arc lengths
    for j = 1:3
        s(j) = sqrt(norm(every_a(j + 1).harmonic_coefficients - every_a(j).harmonic_coefficients));
        if s(j) < 1e-15
            s(j) = 1e-15; % Avoid zero arc length
        end
    end

    % Step 2: Define time parameters t based on arc lengths
    t = zeros(5, 1); % Time parameters
    t(1) = 0;
    t(2) = s(1);
    t(3) = t(2) + s(2);
    t(4) = t(3) + s(3);

    % Adaptive arc length increment for t5
    ds_new = (ds_old * Nd) / iteration_index; % Update arc length increment
    t(5) = t(4) + ds_new;

    % Step 3: Perform cubic interpolation to predict the next solution
    Predictor = zeros(size(every_a(1).harmonic_coefficients)); % Initialize Predictor
    for i = 1:4
        L = 1; % Initialize Lagrange basis
        for j = 1:4
            if i ~= j
                L = L * (t(5) - t(j)) / (t(i) - t(j));
            end
        end
        Predictor = Predictor + L * every_a(i).harmonic_coefficients;
    end
end
