clear;clc;close all;

%% x1_single_stable_fitting
load the_two_parameter_alpha_0.8_alpha1_1_1.033.mat;

%% x1 single 稳定解
num_degrees_freedom = 1;
num_harmonics = 25;
index_bifurcation=1;
index_bifurcation1=1;
%% 周期为2的不稳定解
% for alpha1=1.396:0.001:1.92
for alpha1=1:0.001:1.031
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end

figure;
for i=1:1:length(1:0.001:1.031)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);


%% green x1 single
load x1_alpha1_1.396_0.001_1.86_green.mat;

%% x1 single 稳定解
num_degrees_freedom = 1;
num_harmonics = 25;
index_bifurcation=1;
index_bifurcation1=1;
%% 周期为2的不稳定解
% for alpha1=1.396:0.001:1.92
for alpha1=1.396:0.001:1.86
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end

hold on;
for i=1:1:length(1.396:0.001:1.86)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'g--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

%% green x1 single
clear;
load x1_alpha1_1_0.001_1.031_green.mat;

%% x1 single 稳定解
num_degrees_freedom = 1;
num_harmonics = 25;
index_bifurcation=1;
index_bifurcation1=1;
%% 周期为2的不稳定解
% for alpha1=1.396:0.001:1.92
for alpha1=1:0.001:1.031
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end

hold on;
for i=1:1:length(1:0.001:1.031)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'g-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

clear;
%% x2_single_stable_unstable
load x1_alpha1_1.24_0.001_1.92_singlecycle.mat;

%% x1 single 稳定解
num_degrees_freedom = 1;
num_harmonics = 25;
index_bifurcation=157;
index_bifurcation1=1;
%% 周期为2的不稳定解
% for alpha1=1.396:0.001:1.92
for alpha1=1.396:0.001:1.526
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end

hold on;
for i=1:1:length(1.396:0.001:1.526)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'b-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

%% x1 single 不稳定解
index_bifurcation=287;
index_bifurcation1=1;
%% 周期为2的不稳定解
% for alpha1=1.396:0.001:1.92
for alpha1=1.526:0.001:1.716
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end

for i=1:1:length(1.526:0.001:1.716)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'b--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

%% x1 single 稳定解
index_bifurcation=477;
index_bifurcation1=1;
%% 周期为2的不稳定解
% for alpha1=1.396:0.001:1.92
for alpha1=1.716:0.001:1.874
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end

alpha1_bif_x2=[];x_max_high=[];
for i=1:1:length(1.716:0.001:1.874)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'b-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);




%% x1 double 稳定解
clear;
load parameter_alpha_0.8_alpha1.515_1.730_x1_double.mat;

num_degrees_freedom = 1;
num_harmonics = 25;
index_bifurcation=11;
index_bifurcation1=1;
%% 周期为2的不稳定解
for alpha1=1.525:0.001:1.577
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end

for i=1:1:length(1.525:0.001:1.577)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'k-', 'LineWidth', 1.5);
hold on;
plot(alpha1_bif_x2,x_max_low,'k-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

%% x1 double 不稳定解
index_bifurcation=63;
index_bifurcation1=1;
%% 周期为2的不稳定解
for alpha1=1.577:0.001:1.675
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end

for i=1:1:length(1.577:0.001:1.675)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'k--', 'LineWidth', 1.5);
hold on;
plot(alpha1_bif_x2,x_max_low,'k--', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

%% x1 double 稳定解
index_bifurcation=161;
index_bifurcation1=1;
%% 周期为2的不稳定解
for alpha1=1.675:0.001:1.717
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

    max_displacement(index_bifurcation1)=max(displacement);
    xmax(index_bifurcation1).xmax=getmax(displacement');
    alpha1_bif(index_bifurcation1)=every_harmonic_coefficients(index_bifurcation).alpha1;
    frequency_bif(index_bifurcation1)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation1=index_bifurcation1+1;

end
alpha1_bif_x2=[];x_max_high=[];x_max_low=[];
for i=1:1:length(1.675:0.001:1.717)
    x_max=xmax(i).xmax;
    alpha1_bif_x2(i)=alpha1_bif(i);
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);
end
hold on;
plot(alpha1_bif_x2,x_max_high,'k-', 'LineWidth', 1.5);
hold on;
plot(alpha1_bif_x2,x_max_low,'k-', 'LineWidth', 1.5);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);