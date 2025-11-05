clear; clc; %close all;
% load 1_2_x1_inter_para_alpha1_1.5_alpha2_0.9_doublecycle.mat;
% 
% num_degrees_freedom = 1;
% num_harmonics = 25;
% index_bifurcation=1;
% %% 周期为2的不稳定解
% for alpha=0.6:0.002:0.684
%     %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
%     harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
%     index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;
% 
%     initial_frequency = harmonic_coefficients(1, 1);
%     total_time = 4 * pi / initial_frequency;
%     step_size = 2 * pi / (initial_frequency * 500);
%     time_data = (0:step_size:total_time);
%     harmonic_parameters = harmonic_coefficients(3:end, :);
% 
%     % Initialize arrays for displacement, velocity, and acceleration
%     displacement = zeros(num_degrees_freedom, length(time_data));
%     velocity = zeros(num_degrees_freedom, length(time_data));
%     acceleration = zeros(num_degrees_freedom, length(time_data));
% 
%     % Compute displacement, velocity, and acceleration for each degree of freedom
%     for dof_index = 1:num_degrees_freedom
%         for harm_index = 1:num_harmonics
%             harmonic_term = harm_index  * initial_frequency * time_data;
%             displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%             velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);
%             acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%         end
%         displacement(dof_index, :) = displacement(dof_index, :) + harmonic_coefficients(2, 2*dof_index-1);
%     end
% 
%     max_displacement(index_bifurcation)=max(displacement);
%     xmax(index_bifurcation).xmax=getmax(displacement');
%     alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation).alpha;
%     frequency_bif(index_bifurcation)=initial_frequency;
%     index_bifurcation=index_bifurcation+1;
%     
% end
% 
% figure;
% % for i=1:2:length(0.6:0.002:0.684)
% %     QQ=alpha_bif(i)*ones(1,length(xmax(i).xmax));
% %     plot(QQ,xmax(i).xmax,'k.','MarkerSize',6);
% %     hold on;
% % end
% 
% for i=1:1:length(0.6:0.002:0.684)
%     x_max=xmax(i).xmax;
%     x_max_high(i)=max(x_max);
%     x_max_low(i)=min(x_max);   
% end
% hold on;
% plot(alpha_bif,x_max_high,'k--', 'LineWidth', 1);
% hold on;
% plot(alpha_bif,x_max_low,'k--', 'LineWidth', 1);
% 
% alpha_bif=[];
% index_bifurcation_alpha=44;
% index_bifurcation=1;
% %% 周期为2的稳定解
% for alpha=0.686:0.002:0.758
%     %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
%     harmonic_coefficients=every_harmonic_coefficients(index_bifurcation_alpha).harmonic_coefficients;
%     index_residual=every_harmonic_coefficients(index_bifurcation_alpha).index_residual;
% 
%     initial_frequency = harmonic_coefficients(1, 1);
%     total_time = 4 * pi / initial_frequency;
%     step_size = 2 * pi / (initial_frequency * 500);
%     time_data = (0:step_size:total_time);
%     harmonic_parameters = harmonic_coefficients(3:end, :);
% 
%     % Initialize arrays for displacement, velocity, and acceleration
%     displacement = zeros(num_degrees_freedom, length(time_data));
%     velocity = zeros(num_degrees_freedom, length(time_data));
%     acceleration = zeros(num_degrees_freedom, length(time_data));
% 
%     % Compute displacement, velocity, and acceleration for each degree of freedom
%     for dof_index = 1:num_degrees_freedom
%         for harm_index = 1:num_harmonics
%             harmonic_term = harm_index  * initial_frequency * time_data;
%             displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%             velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);
%             acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%         end
%         displacement(dof_index, :) = displacement(dof_index, :) + harmonic_coefficients(2, 2*dof_index-1);
%     end
% 
%     max_displacement(index_bifurcation)=max(displacement);
%     xmax(index_bifurcation).xmax=getmax(displacement');
%     alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation_alpha).alpha;
%     frequency_bif(index_bifurcation)=initial_frequency;
%     index_bifurcation=index_bifurcation+1;
%     index_bifurcation_alpha=index_bifurcation_alpha+1;
% end
% 
% x_max_high=[];
% x_max_low=[];
% for i=1:1:length(0.686:0.002:0.758)
%     x_max=xmax(i).xmax;
%     x_max_high(i)=max(x_max);
%     x_max_low(i)=min(x_max);   
% end
% hold on;
% plot(alpha_bif,x_max_high,'k-', 'LineWidth', 1);
% hold on;
% plot(alpha_bif,x_max_low,'k-', 'LineWidth', 1);
% % hold on;
% % for i=1:2:length(0.686:0.002:0.758)
% %     QQ=alpha_bif(i)*ones(1,length(xmax(i).xmax));
% %     plot(QQ,xmax(i).xmax,'b.','MarkerSize',6);
% %     hold on;
% % end
% 
% %% 周期为1的稳定解
% index_bifurcation_alpha=80;
% xmax=[];alpha_bif=[];frequency_bif=[];
% index_bifurcation=1;
% for alpha=0.758:0.002:0.9
%     %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
%     harmonic_coefficients=every_harmonic_coefficients(index_bifurcation_alpha).harmonic_coefficients;
%     index_residual=every_harmonic_coefficients(index_bifurcation_alpha).index_residual;
% 
%     initial_frequency = harmonic_coefficients(1, 1);
%     total_time = 2 * pi / initial_frequency;
%     step_size = 2 * pi / (initial_frequency * 500);
%     time_data = (0:step_size:total_time);
%     harmonic_parameters = harmonic_coefficients(3:end, :);
% 
%     % Initialize arrays for displacement, velocity, and acceleration
%     displacement = zeros(num_degrees_freedom, length(time_data));
%     velocity = zeros(num_degrees_freedom, length(time_data));
%     acceleration = zeros(num_degrees_freedom, length(time_data));
% 
%     % Compute displacement, velocity, and acceleration for each degree of freedom
%     for dof_index = 1:num_degrees_freedom
%         for harm_index = 1:num_harmonics
%             harmonic_term = harm_index  * initial_frequency * time_data;
%             displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%             velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);
%             acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%         end
%         displacement(dof_index, :) = displacement(dof_index, :) + harmonic_coefficients(2, 2*dof_index-1);
%     end
% 
%     max_displacement(index_bifurcation)=max(displacement);
%     xmax(index_bifurcation)=getmax(displacement');
%     alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation_alpha).alpha;
%     frequency_bif(index_bifurcation)=initial_frequency;
%     index_bifurcation=index_bifurcation+1;
%     index_bifurcation_alpha=index_bifurcation_alpha+1;
% end
% 
% hold on;
% for i=1:1:length(0.758:0.002:0.9)
%     %QQ=alpha_bif(i)*ones(1,length(xmax(i).xmax));
%     plot(alpha_bif,xmax,'b-', 'LineWidth', 1);
%     hold on;
% end
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% 
% 
% %% 周期为1的不稳定解
% clear;
% load 1_2_x1_inter_para_alpha1_1.5_alpha2_0.9_singlecycle.mat;
% 
% num_degrees_freedom = 1;
% num_harmonics = 25;
% index_bifurcation=1;
% for alpha=0.6:0.002:0.758
%     %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
%     harmonic_coefficients=every_harmonic_coefficients(index_bifurcation).harmonic_coefficients;
%     index_residual=every_harmonic_coefficients(index_bifurcation).index_residual;
% 
%     initial_frequency = harmonic_coefficients(1, 1);
%     total_time = 4 * pi / initial_frequency;
%     step_size = 2 * pi / (initial_frequency * 500);
%     time_data = (0:step_size:total_time);
%     harmonic_parameters = harmonic_coefficients(3:end, :);
% 
%     % Initialize arrays for displacement, velocity, and acceleration
%     displacement = zeros(num_degrees_freedom, length(time_data));
%     velocity = zeros(num_degrees_freedom, length(time_data));
%     acceleration = zeros(num_degrees_freedom, length(time_data));
% 
%     % Compute displacement, velocity, and acceleration for each degree of freedom
%     for dof_index = 1:num_degrees_freedom
%         for harm_index = 1:num_harmonics
%             harmonic_term = harm_index  * initial_frequency * time_data;
%             displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%             velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);
%             acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%         end
%         displacement(dof_index, :) = displacement(dof_index, :) + harmonic_coefficients(2, 2*dof_index-1);
%     end
% 
%     max_displacement(index_bifurcation)=max(displacement);
%     xmax(index_bifurcation)=getmax(displacement');
%     alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation).alpha;
%     frequency_bif(index_bifurcation)=initial_frequency;
%     index_bifurcation=index_bifurcation+1;
% end
% 
% hold on;
% for i=1:1:length(0.6:0.002:0.758)
%     plot(alpha_bif,xmax,'b--', 'LineWidth', 1);
%     hold on;
% end

%% x2 
%% 周期为2的不稳定解
clear;
load 1_2_x2_inter_para_alpha1_1.5_alpha2_0.9_doublecycle.mat;
num_degrees_freedom = 1;
num_harmonics = 25;

index_bifurcation=1;
for alpha=0.6:0.002:0.748
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
    alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation).alpha;
    frequency_bif(index_bifurcation)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
end

x_max_high=[];
x_max_low=[];
for i=1:1:length(0.6:0.002:0.748)
    x_max=xmax(i).xmax;
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);   
end
%% 周期为2的不稳定解
figure;
plot(alpha_bif,x_max_high,'k--', 'LineWidth', 1);
hold on;
plot(alpha_bif,x_max_low,'k--', 'LineWidth', 1);

% for i=1:2:length(0.6:0.002:0.748)
%     QQ=alpha_bif(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'k.','MarkerSize',6);
%     hold on;
% end

%% 周期为2的稳定解
index_bifurcation_alpha=76;
index_bifurcation=1;alpha_bif=[];
for alpha=0.75:0.002:0.842
    %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
    harmonic_coefficients=every_harmonic_coefficients(index_bifurcation_alpha).harmonic_coefficients;
    index_residual=every_harmonic_coefficients(index_bifurcation_alpha).index_residual;

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
    alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation_alpha).alpha;
    frequency_bif(index_bifurcation)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation_alpha=index_bifurcation_alpha+1;
end

x_max_high=[];
x_max_low=[];
for i=1:1:length(0.75:0.002:0.842)
    x_max=xmax(i).xmax;
    x_max_high(i)=max(x_max);
    x_max_low(i)=min(x_max);   
end
%% 周期为2的稳定解
hold on;
plot(alpha_bif,x_max_high,'k-', 'LineWidth', 1);
hold on;
plot(alpha_bif,x_max_low,'k-', 'LineWidth', 1);

% figure;
% for i=1:2:length(0.75:0.002:0.842)
%     QQ=alpha_bif(i)*ones(1,length(xmax(i).xmax));
%     plot(QQ,xmax(i).xmax,'k.','MarkerSize',6);
%     hold on;
% end

%% 周期为1的稳定解
index_bifurcation_alpha=122;
xmax=[];alpha_bif=[];frequency_bif=[];
index_bifurcation=1;
for alpha=0.842:0.002:0.9
    %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
    harmonic_coefficients=every_harmonic_coefficients(index_bifurcation_alpha).harmonic_coefficients;
    index_residual=every_harmonic_coefficients(index_bifurcation_alpha).index_residual;

    initial_frequency = harmonic_coefficients(1, 1);
    total_time = 2 * pi / initial_frequency;
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
    xmax(index_bifurcation)=getmax(displacement');
    alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation_alpha).alpha;
    frequency_bif(index_bifurcation)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
    index_bifurcation_alpha=index_bifurcation_alpha+1;
end

hold on;
for i=1:1:length(0.842:0.002:0.9)
    %QQ=alpha_bif(i)*ones(1,length(xmax(i).xmax));
    plot(alpha_bif,xmax,'b-', 'LineWidth', 1);
    hold on;
end
set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);

%% 周期为1的不稳定解
clear;
load 1_2_x2_inter_para_alpha1_1.5_alpha2_0.9_singlecycle.mat;

num_degrees_freedom = 1;
num_harmonics = 25;
index_bifurcation=1;
for alpha=0.6:0.002:0.842
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
    xmax(index_bifurcation)=getmax(displacement');
    alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation).alpha;
    frequency_bif(index_bifurcation)=initial_frequency;
    index_bifurcation=index_bifurcation+1;
end

hold on;
for i=1:1:length(0.6:0.002:0.842)
    plot(alpha_bif,xmax,'b--', 'LineWidth', 1);
    hold on;
end



% signgle;

% index_bifurcation_alpha=80;
% xmax=[];alpha_bif=[];frequency_bif=[];
% index_bifurcation=1;
% for alpha=0.758:0.002:0.9
%     %alpha=every_harmonic_coefficients(index_bifurcation).alpha;
%     harmonic_coefficients=every_harmonic_coefficients(index_bifurcation_alpha).harmonic_coefficients;
%     index_residual=every_harmonic_coefficients(index_bifurcation_alpha).index_residual;
% 
%     initial_frequency = harmonic_coefficients(1, 1);
%     total_time = 4 * pi / initial_frequency;
%     step_size = 2 * pi / (initial_frequency * 500);
%     time_data = (0:step_size:total_time);
%     harmonic_parameters = harmonic_coefficients(3:end, :);
% 
%     % Initialize arrays for displacement, velocity, and acceleration
%     displacement = zeros(num_degrees_freedom, length(time_data));
%     velocity = zeros(num_degrees_freedom, length(time_data));
%     acceleration = zeros(num_degrees_freedom, length(time_data));
% 
%     % Compute displacement, velocity, and acceleration for each degree of freedom
%     for dof_index = 1:num_degrees_freedom
%         for harm_index = 1:num_harmonics
%             harmonic_term = harm_index  * initial_frequency * time_data;
%             displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%             velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);
%             acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
%         end
%         displacement(dof_index, :) = displacement(dof_index, :) + harmonic_coefficients(2, 2*dof_index-1);
%     end
% 
%     max_displacement(index_bifurcation)=max(displacement);
%     xmax(index_bifurcation)=getmax(displacement');
%     alpha_bif(index_bifurcation)=every_harmonic_coefficients(index_bifurcation_alpha).alpha;
%     frequency_bif(index_bifurcation)=initial_frequency;
%     index_bifurcation=index_bifurcation+1;
%     index_bifurcation_alpha=index_bifurcation_alpha+1;
% end
% 
% hold on;
% for i=1:1:length(0.758:0.002:0.9)
%     %QQ=alpha_bif(i)*ones(1,length(xmax(i).xmax));
%     plot(alpha_bif,xmax,'b-', 'LineWidth', 1);
%     hold on;
% end
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);










% Plot displacement response over time
% figure;
% plot(alpha_bif, max_displacement, 'k-', 'LineWidth', 1);
% legend_handle = legend('$$x_1$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);
% 
% figure;
% plot(alpha_bif, frequency_bif, 'k-', 'LineWidth', 1);
% legend_handle = legend('$$\omega$$');
% set(legend_handle, 'Interpreter', 'latex', 'FontSize', 15);
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 15, 'LineWidth', 1.5);