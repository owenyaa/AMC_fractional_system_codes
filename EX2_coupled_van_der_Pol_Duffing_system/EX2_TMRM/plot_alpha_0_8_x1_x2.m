%% x2
%% 周期为2的不稳定解
clear;clc;close all;
load 1_2_x1_inter_para_alpha1_1.5_alpha2_0.9_doublecycle.mat;
num_degrees_freedom = 1;
num_harmonics = 25;
alpha=0.8;
index_bifurcation=101;
% for alpha=0.6:0.002:0.748
%alpha=every_harmonic_coefficients(index_bifurcation).alpha;
harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;

initial_frequency = harmonic_coefficients(1, 1);
total_time = pi / initial_frequency;
step_size = 2 * pi / (initial_frequency * 500);
time_data = (0:step_size:total_time);
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

%% 周期为2的不稳定解
subplot(2,1,1);
plot(time_data,displacement,'k-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

figure;
plot(displacement,velocity,'k-', 'LineWidth', 1);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


clear;clc;%close all;
load 1_2_x2_inter_para_alpha1_1.5_alpha2_0.9_doublecycle.mat;
num_degrees_freedom = 1;
num_harmonics = 25;
alpha=0.8;
index_bifurcation=101;
% for alpha=0.6:0.002:0.748
%alpha=every_harmonic_coefficients(index_bifurcation).alpha;
harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;

initial_frequency = harmonic_coefficients(1, 1);
total_time =  2 * pi / initial_frequency;
step_size = 2 * pi / (initial_frequency * 500);
time_data = (0:step_size:total_time);
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

%% 周期为2的不稳定解
% hold on;
% subplot(2,1,2);
% plot(time_data,displacement,'k-', 'LineWidth', 1);
% legend_handle = legend('$$x_2$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

figure;
plot(displacement,velocity,'k-', 'LineWidth', 1);
legend_handle = legend('$$x_2$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

clear;clc;%close all;
load x_2_num_all.mat;
index_bifurcation=101;

x_cal=num_x_cal(index_bifurcation).x_cal;
tf=800;h=0.01;Tdata=(0:h:tf)';


% figure;
% subplot(2,1,1);
% plot(Tdata(1:dd:end),x_cal(1:dd:end-1,1),'k.', 'MarkerSize', 6);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);




dd=20;
figure;
% subplot(2,1,1);
plot(x_cal(79150:dd:end-1,1),x_cal(79150:dd:end-1,3),'k.', 'MarkerSize', 6);
legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

% dd=30;
% figure;
% subplot(2,1,2);
% plot(Tdata(1:dd:end),x_cal(1:dd:end-1,2),'k.', 'MarkerSize', 6);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);



