%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 22:13
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 2 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%          \begin{equation}%
%          	\label{eq4.15}%
%          	\begin{cases}%
%          		\begin{aligned}%
%          		  & \ddot{\xi}+0.25\ddot{\alpha}+0.1\dot{\xi}+0.2\xi+0.1Q\alpha=0           \\%
%          		  & 0.25\ddot{\xi}+0.5\ddot{\alpha}+0.1\dot{\alpha}-0.04Q\alpha+f(\alpha)=0 %
%          		\end{aligned}%
%          	\end{cases}%
%          \end{equation}%
%          The structure of the program is as follows:%
%          The basis function in this program is
%          \begin{equation}
%          	\label{eq4.17}
%          	\begin{cases}
%         		\begin{aligned}
%         		  & x_1\approx x^N_1=\sum_{k=1}^{N}\left[b_k\cos\left((2k-1)\omega t\right)+c_k\sin\left((2k-1)\omega t\right)\right] \\
%         		  & x_2\approx x^N_2=\sum_{k=1}^{N}\left[d_k\cos\left((2k-1)\omega t\right)+e_k\sin\left((2k-1)\omega t\right)\right]
%         		\end{aligned}
%         	\end{cases}
%          \end{equation}
%          The function of this program is to calculate the residuals, including the%
%          residuals of the control equations, the residuals caused by the initial value
%          conditions, and the sensitivity response of the residuals with %
%          respect to the coefficients.%
function residual = calculate_residual(harmonic_params)
global num_degrees_freedom num_harmonics time_data
global alpha epsilon omega_2 nonliear_stiffness_matrix mass_matrix parameter_constant

initial_frequency = harmonic_params(1, 1);
harmonic_parameters = harmonic_params(2:end, :);

displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
Dx_alpha=zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        harmonic_term = (2 * harm_index - 1) * initial_frequency * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
        velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);

        Y1_N_harm=(2*harm_index-1)^alpha*cos(alpha*pi/2);Y2_N_harm=(2*harm_index-1)^alpha*sin(alpha*pi/2);
        Dx_alpha(dof_index,:)=Dx_alpha(dof_index,:)+initial_frequency^alpha*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm*harmonic_parameters(harm_index,2*dof_index-1))*sin(harmonic_term)+...
            initial_frequency^alpha*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_parameters(harm_index,2*dof_index))*cos(harmonic_term);

        acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
    end
end

residual(1:num_degrees_freedom, :) = mass_matrix * acceleration + epsilon*(1-parameter_constant*displacement.^2).*Dx_alpha + omega_2*displacement + nonliear_stiffness_matrix * displacement.^3;

% Initialize arrays for displacement, velocity, and acceleration sensitivity to frequency
displacement_freq_sensitivity = zeros(num_degrees_freedom, length(time_data));
velocity_freq_sensitivity = zeros(num_degrees_freedom, length(time_data));
Dx_alpha_w0 = zeros(num_degrees_freedom, length(time_data));
acceleration_freq_sensitivity = zeros(num_degrees_freedom, length(time_data));

% Compute the frequency sensitivity of displacement, velocity, and acceleration for each degree of freedom
for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        % Compute the frequency derivative term for each harmonic
        frequency_derivative = (2 * harm_index - 1) * time_data;
        % Calculate the frequency sensitivity of displacement
        displacement_freq_sensitivity(dof_index, :) = displacement_freq_sensitivity(dof_index, :) - harmonic_parameters(harm_index, 2 * dof_index - 1) * frequency_derivative .* sin((2 * harm_index - 1) * initial_frequency * time_data) + harmonic_parameters(harm_index, 2 * dof_index) * frequency_derivative .* cos((2 * harm_index - 1) * initial_frequency * time_data);

        % Calculate frequency sensitivity of velocity
        velocity_freq_sensitivity(dof_index, :) = velocity_freq_sensitivity(dof_index, :) - ((2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin((2 * harm_index - 1) * initial_frequency * time_data) + initial_frequency * time_data * (2 * harm_index - 1)^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) .* cos((2 * harm_index - 1) * initial_frequency * time_data)) +((2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index) * cos((2 * harm_index - 1) * initial_frequency * time_data) - initial_frequency * time_data * (2 * harm_index - 1)^2 * harmonic_parameters(harm_index, 2 * dof_index) .* sin((2 * harm_index - 1) * initial_frequency * time_data));

        Y1_N_harm=(2*harm_index-1)^alpha*cos(alpha*pi/2);Y2_N_harm=(2*harm_index-1)^alpha*sin(alpha*pi/2);
        Dx_alpha_w0(dof_index,:)=Dx_alpha_w0(dof_index,:)+initial_frequency^alpha*frequency_derivative*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm*harmonic_parameters(harm_index,2*dof_index-1)).*cos((2*harm_index-1)*initial_frequency*time_data)+...
            alpha*initial_frequency^(alpha-1)*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm*harmonic_parameters(harm_index,2*dof_index-1)).*sin((2*harm_index-1)*initial_frequency*time_data)+...
            -initial_frequency^alpha*frequency_derivative*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_parameters(harm_index,2*dof_index)).*sin((2*harm_index-1)*initial_frequency*time_data)+...
            alpha*initial_frequency^(alpha-1)*(Y1_N_harm*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_parameters(harm_index,2*dof_index)).*cos((2*harm_index-1)*initial_frequency*time_data);

        % Calculate frequency sensitivity of acceleration
        acceleration_freq_sensitivity(dof_index, :) = acceleration_freq_sensitivity(dof_index, :) - (2 * initial_frequency * (2 * harm_index - 1)^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos((2 * harm_index - 1) * initial_frequency * time_data) - initial_frequency^2 * time_data * (2 * harm_index - 1)^3 * harmonic_parameters(harm_index, 2 * dof_index - 1) .* sin((2 * harm_index - 1) * initial_frequency * time_data)) - (2 * initial_frequency * (2 * harm_index - 1)^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin((2 * harm_index - 1) * initial_frequency * time_data) + initial_frequency^2 * time_data * (2 * harm_index - 1)^3 * harmonic_parameters(harm_index, 2 * dof_index) .* cos((2 * harm_index - 1) * initial_frequency * time_data));

    end
end

% Compute the residual considering the frequency sensitivity
residual(num_degrees_freedom+1:2*num_degrees_freedom,:) = mass_matrix * acceleration_freq_sensitivity + epsilon.*Dx_alpha_w0 -...
    epsilon*parameter_constant*(2*displacement.*Dx_alpha.*displacement_freq_sensitivity+Dx_alpha_w0.*displacement.^2) + omega_2*displacement_freq_sensitivity + 3*nonliear_stiffness_matrix * displacement.^2.*displacement_freq_sensitivity;

for i = 1:2 * num_harmonics * num_degrees_freedom
    sensitivity_params = zeros(2 * num_harmonics * num_degrees_freedom, 1);
    sensitivity_params(i, 1) = 1;
    reshaped_params = reshape(sensitivity_params, 2, num_harmonics * num_degrees_freedom);
    reshaped_params = reshaped_params';

    harmonic_sensitivity = reshaped_params(1:num_harmonics, 1:2);

    for num_dof = 1:num_degrees_freedom - 1
        harmonic_sensitivity = [harmonic_sensitivity, reshaped_params(num_dof * num_harmonics + 1:(num_dof + 1) * num_harmonics, 1:2)];
    end

    displacement_harm_sensitivity = zeros(num_degrees_freedom, length(time_data));
    velocity_harm_sensitivity = zeros(num_degrees_freedom, length(time_data));
    Dx_alpha_a = zeros(num_degrees_freedom, length(time_data));
    acceleration_harm_sensitivity = zeros(num_degrees_freedom, length(time_data));

    for dof_index = 1:num_degrees_freedom
        for harm_index = 1:num_harmonics
            harmonic_term = (2 * harm_index - 1) * initial_frequency * time_data;
            displacement_harm_sensitivity(dof_index, :) = displacement_harm_sensitivity(dof_index, :) + harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);

            Y1_N_harm=(2*harm_index-1)^alpha*cos(alpha*pi/2);Y2_N_harm=(2*harm_index-1)^alpha*sin(alpha*pi/2);
            Dx_alpha_a(dof_index,:)=Dx_alpha_a(dof_index,:)+initial_frequency^alpha*(Y1_N_harm*harmonic_sensitivity(harm_index,2*dof_index)-Y2_N_harm*harmonic_sensitivity(harm_index,2*dof_index-1))*sin(harmonic_term)+...
                initial_frequency^alpha*(Y1_N_harm*harmonic_sensitivity(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_sensitivity(harm_index,2*dof_index))*cos(harmonic_term);

            velocity_harm_sensitivity(dof_index, :) = velocity_harm_sensitivity(dof_index, :) - initial_frequency * (2 * harm_index - 1) * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * (2 * harm_index - 1) * harmonic_sensitivity(harm_index, 2 * dof_index) * cos(harmonic_term);
            acceleration_harm_sensitivity(dof_index, :) = acceleration_harm_sensitivity(dof_index, :) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);
        end
    end

    residual(num_degrees_freedom*(i+1)+1:num_degrees_freedom*(i+2),:) = mass_matrix * acceleration_harm_sensitivity + epsilon.*Dx_alpha_a -...
        epsilon*parameter_constant*(2*displacement.*Dx_alpha.*displacement_harm_sensitivity+Dx_alpha_a.*displacement.^2) + omega_2*displacement_harm_sensitivity + 3*nonliear_stiffness_matrix * displacement.^2.*displacement_harm_sensitivity;

end
residual = residual';
end


