%% Author : GUANG_LIU  * owenyaa@gmail.com *
% Created Time : 2020-11-12 22:13
% Last Revised : GUANG_LIU ,2020-11-12
% Remark : The structure, function and variables of this program are explained:%
%          This program is to calculate the example 2 in the paper "A New Semi-%
%          analytical technique for Nonlinear Systems Based on Response Sensitivity%
%           Analysis", that is, the semi-analytical solution of equation %
%          \begin{equation}
% 		   \label{eq3.1}
% 		        D^\alpha1+\epsilon(1-x^2)D^\alpha2\dot{x}+kx^3=0
% 	       \end{equation}
%          The structure of the program is as follows:%
%          The basis function in this program is
%          \begin{equation}
%          	\label{eq4.17}
%          	x_1\approx x^N_1=\sum_{k=1}^{N}\left[b_k\cos\left((2k-1)\omega t\right)+c_k\sin\left((2k-1)\omega t\right)\right] 
%          \end{equation}
%          The function of this program is to calculate the residuals, including the%
%          residuals of the control equations, the residuals caused by the initial value
%          conditions, and the sensitivity response of the residuals with %
%          respect to the coefficients.%
function residual = calculate_residual(harmonic_params)
global num_degrees_freedom num_harmonics time_data
global alpha1 alpha2 epsilon omega_2 nonliear_stiffness_matrix mass_matrix parameter_constant
beta=alpha1-1;
initial_frequency = harmonic_params(1, 1);
harmonic_parameters = harmonic_params(2:end, :);

displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
Dx_alpha1=zeros(num_degrees_freedom, length(time_data));
Dx_alpha2=zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));

for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        harmonic_term = (2 * harm_index - 1) * initial_frequency * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
        velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * (2 * harm_index - 1) * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);

        Y1_N_harm_alpha1=(2*harm_index-1)^alpha1*cos(beta*pi/2);Y2_N_harm_alpha1=(2*harm_index-1)^alpha1*sin(beta*pi/2);
        Dx_alpha1(dof_index,:)=Dx_alpha1(dof_index,:)+initial_frequency^alpha1*(Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1))*cos(harmonic_term)+...
        initial_frequency^alpha1*(-Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index))*sin(harmonic_term);

        Y1_N_harm_alpha2=(2*harm_index-1)^alpha2*cos(alpha2*pi/2);Y2_N_harm_alpha2=(2*harm_index-1)^alpha2*sin(alpha2*pi/2);
        Dx_alpha2(dof_index,:)=Dx_alpha2(dof_index,:)+initial_frequency^alpha2*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1))*sin(harmonic_term)+...
            initial_frequency^alpha2*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index))*cos(harmonic_term);

        acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
    end
end

residual(1:num_degrees_freedom, :) = mass_matrix * Dx_alpha1 + epsilon*(1-parameter_constant*displacement.^2).*Dx_alpha2 + omega_2*displacement + nonliear_stiffness_matrix * displacement.^3;
                                         
% Initialize arrays for displacement, velocity, and acceleration sensitivity to frequency
displacement_freq_sensitivity = zeros(num_degrees_freedom, length(time_data));
velocity_freq_sensitivity = zeros(num_degrees_freedom, length(time_data));
Dx_alpha1_w0 = zeros(num_degrees_freedom, length(time_data));
Dx_alpha2_w0 = zeros(num_degrees_freedom, length(time_data));
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

        Y1_N_harm_alpha1=(2*harm_index-1)^alpha1*cos(beta*pi/2);Y2_N_harm_alpha1=(2*harm_index-1)^alpha1*sin(beta*pi/2);
        Dx_alpha1_w0(dof_index,:)=Dx_alpha1_w0(dof_index,:)+initial_frequency^alpha1*frequency_derivative*(-Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)+Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)).*sin((2*harm_index-1)*initial_frequency*time_data)+...
        alpha1*initial_frequency^(alpha1-1)*(Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)).*cos((2*harm_index-1)*initial_frequency*time_data)+...
        +initial_frequency^alpha1*frequency_derivative*(-Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)).*cos((2*harm_index-1)*initial_frequency*time_data)+...
        alpha1*initial_frequency^(alpha1-1)*(-Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)).*sin((2*harm_index-1)*initial_frequency*time_data);

        Y1_N_harm_alpha2=(2*harm_index-1)^alpha2*cos(alpha2*pi/2);Y2_N_harm_alpha2=(2*harm_index-1)^alpha2*sin(alpha2*pi/2);
        Dx_alpha2_w0(dof_index,:)=Dx_alpha2_w0(dof_index,:)+initial_frequency^alpha2*frequency_derivative*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1)).*cos((2*harm_index-1)*initial_frequency*time_data)+...
            alpha2*initial_frequency^(alpha2-1)*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1)).*sin((2*harm_index-1)*initial_frequency*time_data)+...
            -initial_frequency^alpha2*frequency_derivative*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index)).*sin((2*harm_index-1)*initial_frequency*time_data)+...
            alpha2*initial_frequency^(alpha2-1)*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index)).*cos((2*harm_index-1)*initial_frequency*time_data);

        % Calculate frequency sensitivity of acceleration
        acceleration_freq_sensitivity(dof_index, :) = acceleration_freq_sensitivity(dof_index, :) - (2 * initial_frequency * (2 * harm_index - 1)^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos((2 * harm_index - 1) * initial_frequency * time_data) - initial_frequency^2 * time_data * (2 * harm_index - 1)^3 * harmonic_parameters(harm_index, 2 * dof_index - 1) .* sin((2 * harm_index - 1) * initial_frequency * time_data)) - (2 * initial_frequency * (2 * harm_index - 1)^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin((2 * harm_index - 1) * initial_frequency * time_data) + initial_frequency^2 * time_data * (2 * harm_index - 1)^3 * harmonic_parameters(harm_index, 2 * dof_index) .* cos((2 * harm_index - 1) * initial_frequency * time_data));

    end
end

% Compute the residual considering the frequency sensitivity
residual(num_degrees_freedom+1:2*num_degrees_freedom,:) = mass_matrix * Dx_alpha1_w0 + epsilon.*Dx_alpha2_w0 -...
    epsilon*parameter_constant*(2*displacement.*Dx_alpha2.*displacement_freq_sensitivity+Dx_alpha2_w0.*displacement.^2) + omega_2*displacement_freq_sensitivity + 3*nonliear_stiffness_matrix * displacement.^2.*displacement_freq_sensitivity;
                  
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
    Dx_alpha1_a = zeros(num_degrees_freedom, length(time_data));
    Dx_alpha2_a = zeros(num_degrees_freedom, length(time_data));
    acceleration_harm_sensitivity = zeros(num_degrees_freedom, length(time_data));

    for dof_index = 1:num_degrees_freedom
        for harm_index = 1:num_harmonics
            harmonic_term = (2 * harm_index - 1) * initial_frequency * time_data;
            displacement_harm_sensitivity(dof_index, :) = displacement_harm_sensitivity(dof_index, :) + harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);

            Y1_N_harm_alpha1=(2*harm_index-1)^alpha1*cos(beta*pi/2);Y2_N_harm_alpha1=(2*harm_index-1)^alpha1*sin(beta*pi/2);
            Dx_alpha1_a(dof_index,:)=Dx_alpha1_a(dof_index,:)+initial_frequency^alpha1*(Y1_N_harm_alpha1*harmonic_sensitivity(harm_index,2*dof_index)-Y2_N_harm_alpha1*harmonic_sensitivity(harm_index,2*dof_index-1))*cos(harmonic_term)+...
                       initial_frequency^alpha1*(-Y1_N_harm_alpha1*harmonic_sensitivity(harm_index,2*dof_index-1)-Y2_N_harm_alpha1*harmonic_sensitivity(harm_index,2*dof_index))*sin(harmonic_term);        

            Y1_N_harm_alpha2=(2*harm_index-1)^alpha2*cos(alpha2*pi/2);Y2_N_harm_alpha2=(2*harm_index-1)^alpha2*sin(alpha2*pi/2);
            Dx_alpha2_a(dof_index,:)=Dx_alpha2_a(dof_index,:)+initial_frequency^alpha2*(Y1_N_harm_alpha2*harmonic_sensitivity(harm_index,2*dof_index)-Y2_N_harm_alpha2*harmonic_sensitivity(harm_index,2*dof_index-1))*sin(harmonic_term)+...
                initial_frequency^alpha2*(Y1_N_harm_alpha2*harmonic_sensitivity(harm_index,2*dof_index-1)+Y2_N_harm_alpha2*harmonic_sensitivity(harm_index,2*dof_index))*cos(harmonic_term);

            velocity_harm_sensitivity(dof_index, :) = velocity_harm_sensitivity(dof_index, :) - initial_frequency * (2 * harm_index - 1) * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency * (2 * harm_index - 1) * harmonic_sensitivity(harm_index, 2 * dof_index) * cos(harmonic_term);
            acceleration_harm_sensitivity(dof_index, :) = acceleration_harm_sensitivity(dof_index, :) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency * (2 * harm_index - 1))^2 * harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);
        end
    end

    residual(num_degrees_freedom*(i+1)+1:num_degrees_freedom*(i+2),:) = mass_matrix * Dx_alpha1_a + epsilon.*Dx_alpha2_a -...
    epsilon*parameter_constant*(2*displacement.*Dx_alpha2.*displacement_harm_sensitivity+Dx_alpha2_a.*displacement.^2) + omega_2*displacement_harm_sensitivity + 3*nonliear_stiffness_matrix * displacement.^2.*displacement_harm_sensitivity;
                 
end
residual = residual';
end


