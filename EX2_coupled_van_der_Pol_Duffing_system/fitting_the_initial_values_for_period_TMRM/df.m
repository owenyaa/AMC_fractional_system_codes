function y=df(harmonic_params)
global num_harmonics time_data %N_harm index_global
num_degrees_freedom1=1;
initial_frequency=harmonic_params(1,1);
harmonic_parameters=harmonic_params(3:end,:);

%% 计算频率的灵敏度  只需要计算基频w2的就可以
% Initialize arrays for displacement, velocity, and acceleration sensitivity to frequency
displacement_freq_sensitivity = zeros(num_degrees_freedom1, length(time_data));
% velocity_freq_sensitivity = zeros(num_degrees_freedom1, length(time_data));
% Dx_alpha_w0 = zeros(num_degrees_freedom1, length(time_data));
% acceleration_freq_sensitivity = zeros(num_degrees_freedom1, length(time_data));

% Compute the frequency sensitivity of displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom1
    for harm_index = 1:num_harmonics
        % Compute the frequency derivative term for each harmonic
        frequency_derivative = harm_index  * time_data;
        % Calculate the frequency sensitivity of displacement
        displacement_freq_sensitivity(dof_index, :) = displacement_freq_sensitivity(dof_index, :) - harmonic_parameters(harm_index, 2 * dof_index - 1) * frequency_derivative .* sin(harm_index  * initial_frequency * time_data) + harmonic_parameters(harm_index, 2 * dof_index) * frequency_derivative .* cos(harm_index  * initial_frequency * time_data);

        % Calculate frequency sensitivity of velocity
        % velocity_freq_sensitivity(dof_index, :) = velocity_freq_sensitivity(dof_index, :) - (harm_index  * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harm_index  * initial_frequency * time_data) + initial_frequency * time_data * harm_index ^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) .* cos(harm_index  * initial_frequency * time_data)) +(harm_index  * harmonic_parameters(harm_index, 2 * dof_index) * cos(harm_index  * initial_frequency * time_data) - initial_frequency * time_data * harm_index ^2 * harmonic_parameters(harm_index, 2 * dof_index) .* sin(harm_index  * initial_frequency * time_data));
        %
        % Y1_N_harm=(2*harm_index-1)^alpha*cos(alpha*pi/2);Y2_N_harm=(2*harm_index-1)^alpha*sin(alpha*pi/2);
        % Dx_alpha_w0(dof_index,:)=Dx_alpha_w0(dof_index,:)+initial_frequency^alpha*(2*harm_index-1)*time_data*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm*harmonic_parameters(harm_index,2*dof_index-1)).*cos((2*harm_index-1)*initial_frequency*time_data)+...
        %     alpha*initial_frequency^(alpha-1)*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm*harmonic_parameters(harm_index,2*dof_index-1)).*sin((2*harm_index-1)*initial_frequency*time_data)+...
        %     -initial_frequency^alpha*(2*harm_index-1)*time_data*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_parameters(harm_index,2*dof_index)).*sin((2*harm_index-1)*initial_frequency*time_data)+...
        %     alpha*initial_frequency^(alpha-1)*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_parameters(harm_index,2*dof_index)).*cos((2*harm_index-1)*initial_frequency*time_data);
        %
        % % Calculate frequency sensitivity of acceleration
        % acceleration_freq_sensitivity(dof_index, :) = acceleration_freq_sensitivity(dof_index, :) - (2 * initial_frequency * harm_index ^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harm_index  * initial_frequency * time_data) - initial_frequency^2 * time_data * harm_index ^3 * harmonic_parameters(harm_index, 2 * dof_index - 1) .* sin(harm_index  * initial_frequency * time_data)) - (2 * initial_frequency * harm_index ^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harm_index  * initial_frequency * time_data) + initial_frequency^2 * time_data * harm_index ^3 * harmonic_parameters(harm_index, 2 * dof_index) .* cos(harm_index  * initial_frequency * time_data));

    end
end
y_frequency=displacement_freq_sensitivity';

%% 计算谐波系数的灵敏度
for i = 1:2 * num_harmonics * num_degrees_freedom1
    sensitivity_params = zeros(2 * num_harmonics * num_degrees_freedom1, 1);
    sensitivity_params(i, 1) = 1;
    reshaped_params = reshape(sensitivity_params, 2, num_harmonics * num_degrees_freedom1);
    reshaped_params = reshaped_params';

    harmonic_sensitivity = reshaped_params(1:num_harmonics, 1:2);

    for num_dof = 1:num_degrees_freedom1 - 1
        harmonic_sensitivity = [harmonic_sensitivity, reshaped_params(num_dof * num_harmonics + 1:(num_dof + 1) * num_harmonics, 1:2)];
    end

    displacement_harm_sensitivity = zeros(num_degrees_freedom1, length(time_data));
    % velocity_harm_sensitivity = zeros(num_degrees_freedom1, length(time_data));
    % Dx_alpha_a = zeros(num_degrees_freedom1, length(time_data));
    % acceleration_harm_sensitivity = zeros(num_degrees_freedom1, length(time_data));

    for dof_index = 1:num_degrees_freedom1
        for harm_index = 1:num_harmonics
            harmonic_term = harm_index  * initial_frequency * time_data;
            displacement_harm_sensitivity(dof_index, :) = displacement_harm_sensitivity(dof_index, :) + harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);

            % Y1_N_harm=(2*harm_index-1)^alpha*cos(alpha*pi/2);Y2_N_harm=(2*harm_index-1)^alpha*sin(alpha*pi/2);
            % Dx_alpha_a(dof_index,:)=Dx_alpha_a(dof_index,:)+initial_frequency^alpha*(Y1_N_harm*harmonic_sensitivity(harm_index,2*dof_index)-Y2_N_harm*harmonic_sensitivity(harm_index,2*dof_index-1))*sin((2*harm_index-1)*initial_frequency*time_data)+...
            %     initial_frequency^alpha*(Y1_N_harm*harmonic_sensitivity(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_sensitivity(harm_index,2*dof_index))*cos((2*harm_index-1)*initial_frequency*time_data);
            %
            % velocity_harm_sensitivity(dof_index, :) = velocity_harm_sensitivity(dof_index, :) - initial_frequency * harm_index  * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_sensitivity(harm_index, 2 * dof_index) * cos(harmonic_term);
            % acceleration_harm_sensitivity(dof_index, :) = acceleration_harm_sensitivity(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);
        end
    end
    y_harmonic(:,i)=displacement_harm_sensitivity';
end


%% 计算谐波系数的灵敏度
% for i = 1:num_degrees_freedom1
%     sensitivity_cons = zeros(2 * num_harmonics * num_degrees_freedom1, 1);
%     %sensitivity_cons(i, 1) = 1;
%     reshaped_cons = reshape(sensitivity_cons, 2, num_harmonics * num_degrees_freedom1);
%     reshaped_cons = reshaped_cons';
%
%     cons_sensitivity = reshaped_cons(1:num_harmonics, 1:2);
%
%     for num_dof = 1:num_degrees_freedom1 - 1
%         reshaped_cons = [reshaped_cons, reshaped_cons(num_dof * num_harmonics + 1:(num_dof + 1) * num_harmonics, 1:2)];
%     end
%
displacement_cons_sensitivity = zeros(num_degrees_freedom1, length(time_data));
%     % velocity_harm_sensitivity = zeros(num_degrees_freedom1, length(time_data));
%     % Dx_alpha_a = zeros(num_degrees_freedom1, length(time_data));
%     % acceleration_harm_sensitivity = zeros(num_degrees_freedom1, length(time_data));
%
%     for dof_index = 1:num_degrees_freedom1
%         for harm_index = 1:num_harmonics
%             harmonic_term = harm_index  * initial_frequency * time_data;
%             displacement_cons_sensitivity(dof_index, :) = displacement_cons_sensitivity(dof_index, :) + reshaped_cons(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + reshaped_cons(harm_index, 2 * dof_index) * sin(harmonic_term);
%
%             % Y1_N_harm=(2*harm_index-1)^alpha*cos(alpha*pi/2);Y2_N_harm=(2*harm_index-1)^alpha*sin(alpha*pi/2);
%             % Dx_alpha_a(dof_index,:)=Dx_alpha_a(dof_index,:)+initial_frequency^alpha*(Y1_N_harm*harmonic_sensitivity(harm_index,2*dof_index)-Y2_N_harm*harmonic_sensitivity(harm_index,2*dof_index-1))*sin((2*harm_index-1)*initial_frequency*time_data)+...
%             %     initial_frequency^alpha*(Y1_N_harm*harmonic_sensitivity(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_sensitivity(harm_index,2*dof_index))*cos((2*harm_index-1)*initial_frequency*time_data);
%             %
%             % velocity_harm_sensitivity(dof_index, :) = velocity_harm_sensitivity(dof_index, :) - initial_frequency * harm_index  * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * harm_index  * harmonic_sensitivity(harm_index, 2 * dof_index) * cos(harmonic_term);
%             % acceleration_harm_sensitivity(dof_index, :) = acceleration_harm_sensitivity(dof_index, :) - (initial_frequency * harm_index )^2 * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * harm_index )^2 * harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);
%         end
%     end
displacement_cons_sensitivity(1, :) = ones(1, length(time_data));
y_cons(:,1)=displacement_cons_sensitivity';
% end

y=[y_frequency,y_harmonic,y_cons];
end