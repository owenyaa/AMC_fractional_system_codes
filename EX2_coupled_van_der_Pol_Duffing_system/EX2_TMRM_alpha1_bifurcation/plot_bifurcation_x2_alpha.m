clear; clc; close all;

%% x2_double_stable
load x2_alpha1_1_0.001_2_doublecycle.mat;

num_degrees_freedom = 1;
num_harmonics = 25;
index_bifurcation=1;
%% 周期为2的不稳定解
for alpha1=1:0.001:1.906
    %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
    harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
    index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;

    initial_frequency = harmonic_coefficients(1, 1);
    total_time = 4 * pi / initial_frequency;
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

    max_displacement(index_bifurcation)=max(displacement);
    xmax(index_bifurcation).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation)=initial_frequency;
    index_bifurcation=index_bifurcation+1;

end

figure;
% for i=1:2:length(1:0.001:1.906)
%     QQ=alpha1_bif(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'k.','MarkerSize',6);
%     hold on;
% end

for i=1:1:length(1:0.001:1.474)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'b-', 'LineWidth', 1.5);
% hold on;
% plot(alpha1_bif_x2,x_max_low,'b-', 'LineWidth', 1.5);

alpha1_bif_x2=[];x_max_high=[];x_max_low=[];
for i=1:1:length(1.474:0.001:1.533)
    x_max=xmax(474+i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(474+i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'k-', 'LineWidth', 1.5);
hold on;
plot(alpha1_bif_x2,x_max_low,'k-', 'LineWidth', 1.5);


alpha1_bif_x2=[];x_max_high=[];x_max_low=[];
for i=1:1:length(1.534:0.001:1.679)
    x_max=xmax(534+i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(534+i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'k--', 'LineWidth', 1.5);
hold on;
plot(alpha1_bif_x2,x_max_low,'k--', 'LineWidth', 1.5);

alpha1_bif_x2=[];x_max_high=[];x_max_low=[];
for i=1:1:length(1.680:0.001:1.724)
    x_max=xmax(680+i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(680+i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'k-', 'LineWidth', 1.5);
hold on;
plot(alpha1_bif_x2,x_max_low,'k-', 'LineWidth', 1.5);























alpha1_bif_x2=[];x_max_high=[];x_max_low=[];
for i=1:1:length(1.724:0.001:1.906)
    x_max=xmax(724+i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(724+i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'b-', 'LineWidth', 1.5);
% 
% 
% 
clear;
%% x2_single_stable_unstable
load x2_alpha1_1_0.001_2_singlecycle.mat;

num_degrees_freedom = 1;
num_harmonics = 25;
index_bifurcation=1;
%% 周期为2的不稳定解
for alpha1=1:0.001:1.906
    %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
    harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
    index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;

    initial_frequency = harmonic_coefficients(1, 1);
    total_time = 4 * pi / initial_frequency;
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

    max_displacement(index_bifurcation)=max(displacement);
    xmax(index_bifurcation).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation)=initial_frequency;
    index_bifurcation=index_bifurcation+1;

end

hold on;
% for i=1:2:length(1:0.001:1.906)
%     QQ=alpha1_bif(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'k.','MarkerSize',6);
%     hold on;
% end

for i=1:1:length(1:0.001:1.474)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'b-', 'LineWidth', 1.5);
% hold on;
% plot(alpha1_bif_x2,x_max_low,'b-', 'LineWidth', 1.5);

alpha1_bif_x2=[];x_max_high=[];x_max_low=[];
for i=1:1:length(1.474:0.001:1.724)
    x_max=xmax(474+i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(474+i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'b--', 'LineWidth', 1.5);
hold on;
plot(alpha1_bif_x2,x_max_low,'b--', 'LineWidth', 1.5);

% alpha1_bif_x2=[];x_max_high=[];x_max_low=[];
% for i=1:1:length(1.724:0.001:1.906)
%     x_max=xmax(724+i).xmax;
%     alpha1_bif_x2(i)=alpha1_bif(724+i);
%     x_max_high(i)=max(x_max);
%     x_max_low(i)=min(x_max);
% end
% hold on;
% plot(alpha1_bif_x2,x_max_high,'b-', 'LineWidth', 1.5);

legend_handle = legend('$$x_1$$');
set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


%% x1_single 
















