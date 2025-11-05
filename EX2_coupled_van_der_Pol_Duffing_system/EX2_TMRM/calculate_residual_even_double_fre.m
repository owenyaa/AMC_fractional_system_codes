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
function residual = calculate_residual_even_double_fre(harmonic_params)
global time_data initial_frequency num_degrees_freedom num_harmonics
global alpha1 alpha2 mass_matrix  stiffness_matrix
global mu k alpha beta epsilon
beta2=alpha1-1;
% beta2=alpha1;

initial_frequency = harmonic_params(1, 1);
initial_frequency_new=initial_frequency/2;
harmonic_parameters = harmonic_params(3:end, :);

displacement = zeros(num_degrees_freedom, length(time_data));
velocity = zeros(num_degrees_freedom, length(time_data));
Dx_alpha1=zeros(num_degrees_freedom, length(time_data));
Dx_alpha2=zeros(num_degrees_freedom, length(time_data));
acceleration = zeros(num_degrees_freedom, length(time_data));


for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        harmonic_term = harm_index * initial_frequency_new * time_data;
        displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
        velocity(dof_index, :) = velocity(dof_index, :) - initial_frequency_new * harm_index * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency_new * harm_index * harmonic_parameters(harm_index, 2 * dof_index) * cos(harmonic_term);

        Y1_N_harm_alpha1=harm_index^alpha1*cos(beta2*pi/2);Y2_N_harm_alpha1=harm_index^alpha1*sin(beta2*pi/2);
        Dx_alpha1(dof_index,:)=Dx_alpha1(dof_index,:)+initial_frequency_new^alpha1*(Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1))*cos(harmonic_term)+...
            initial_frequency_new^alpha1*(-Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index))*sin(harmonic_term);

        Y1_N_harm_alpha2=harm_index^alpha2*cos(alpha2*pi/2);Y2_N_harm=harm_index^alpha2*sin(alpha2*pi/2);
        Dx_alpha2(dof_index,:)=Dx_alpha2(dof_index,:)+initial_frequency_new^alpha2*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm*harmonic_parameters(harm_index,2*dof_index-1))*sin(harmonic_term)+...
            initial_frequency_new^alpha2*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_parameters(harm_index,2*dof_index))*cos(harmonic_term);

        acceleration(dof_index, :) = acceleration(dof_index, :) - (initial_frequency_new * harm_index).^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency_new * harm_index).^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
    end
    displacement(dof_index, :) = displacement(dof_index, :) + harmonic_params(2, 2*dof_index-1);
end

frac_damping_matrix=[0 0;-k k];
nonlinear_part=[-mu*(1-displacement(1,:).^2).*Dx_alpha1(1,:);alpha*Dx_alpha2(2,:)+beta*displacement(2,:).^2+epsilon*displacement(2,:).^3];
residual(1:num_degrees_freedom, :) = mass_matrix * acceleration + frac_damping_matrix * velocity + stiffness_matrix * displacement + nonlinear_part;

velocity_dof1 = zeros(1, length(time_data));
for dof_index = 1
    for harm_index = 1:num_harmonics
        velocity_dof1(dof_index, :) = velocity_dof1(dof_index, :) + initial_frequency_new * harm_index  * harmonic_parameters(harm_index, 2 * dof_index);
    end
end
residual(num_degrees_freedom+1, :) = velocity_dof1;

% Initialize arrays for displacement, velocity, and acceleration sensitivity to frequency
displacement_freq_sensitivity = zeros(num_degrees_freedom, length(time_data));
velocity_freq_sensitivity = zeros(num_degrees_freedom, length(time_data));
Dx_alpha1_w0 = zeros(num_degrees_freedom, length(time_data));
Dx_alpha2_w0 = zeros(num_degrees_freedom, length(time_data));
acceleration_freq_sensitivity = zeros(num_degrees_freedom, length(time_data));

for dof_index = 1:num_degrees_freedom
    for harm_index = 1:num_harmonics
        % Compute the frequency derivative term for each harmonic
        frequency_derivative = harm_index * time_data;
        % Calculate the frequency sensitivity of displacement
        displacement_freq_sensitivity(dof_index, :) = displacement_freq_sensitivity(dof_index, :) - harmonic_parameters(harm_index, 2 * dof_index - 1) * frequency_derivative .* sin(harm_index * initial_frequency_new * time_data) + harmonic_parameters(harm_index, 2 * dof_index) * frequency_derivative .* cos(harm_index * initial_frequency_new * time_data);

        % Calculate frequency sensitivity of velocity
        velocity_freq_sensitivity(dof_index, :) = velocity_freq_sensitivity(dof_index, :) - (harm_index * harmonic_parameters(harm_index, 2 * dof_index - 1) * sin(harm_index * initial_frequency_new * time_data) + initial_frequency_new * time_data * harm_index^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) .* cos(harm_index * initial_frequency_new * time_data)) +(harm_index * harmonic_parameters(harm_index, 2 * dof_index) * cos(harm_index * initial_frequency_new * time_data) - initial_frequency_new * time_data * harm_index^2 * harmonic_parameters(harm_index, 2 * dof_index) .* sin(harm_index * initial_frequency_new * time_data));

        Y1_N_harm_alpha1=harm_index^alpha1*cos(beta2*pi/2);Y2_N_harm_alpha1=harm_index^alpha1*sin(beta2*pi/2);
        Dx_alpha1_w0(dof_index,:)=Dx_alpha1_w0(dof_index,:)+initial_frequency_new^alpha1*frequency_derivative*(-Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)+Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)).*sin(harm_index*initial_frequency_new*time_data)+...
            alpha1*initial_frequency_new^(alpha1-1)*(Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)).*cos(harm_index*initial_frequency_new*time_data)+...
            +initial_frequency_new^alpha1*frequency_derivative*(-Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)).*cos(harm_index*initial_frequency_new*time_data)+...
            alpha1*initial_frequency_new^(alpha1-1)*(-Y1_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index-1)-Y2_N_harm_alpha1*harmonic_parameters(harm_index,2*dof_index)).*sin(harm_index*initial_frequency_new*time_data);

        Y1_N_harm_alpha2=harm_index^alpha2*cos(alpha2*pi/2);Y2_N_harm=harm_index^alpha2*sin(alpha2*pi/2);
        Dx_alpha2_w0(dof_index,:)=Dx_alpha2_w0(dof_index,:)+initial_frequency_new^alpha2*frequency_derivative*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm*harmonic_parameters(harm_index,2*dof_index-1)).*cos(harm_index*initial_frequency_new*time_data)+...
            alpha2*initial_frequency_new^(alpha2-1)*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index)-Y2_N_harm*harmonic_parameters(harm_index,2*dof_index-1)).*sin(harm_index*initial_frequency_new*time_data)+...
            -initial_frequency_new^alpha2*frequency_derivative*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_parameters(harm_index,2*dof_index)).*sin(harm_index*initial_frequency_new*time_data)+...
            alpha2*initial_frequency_new^(alpha2-1)*(Y1_N_harm_alpha2*harmonic_parameters(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_parameters(harm_index,2*dof_index)).*cos(harm_index*initial_frequency_new*time_data);

        % Calculate frequency sensitivity of acceleration
        acceleration_freq_sensitivity(dof_index, :) = acceleration_freq_sensitivity(dof_index, :) - (2 * initial_frequency_new * harm_index^2 * harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harm_index * initial_frequency_new * time_data) - initial_frequency_new^2 * time_data * harm_index^3 * harmonic_parameters(harm_index, 2 * dof_index - 1) .* sin(harm_index * initial_frequency_new * time_data)) - (2 * initial_frequency_new * harm_index^2 * harmonic_parameters(harm_index, 2 * dof_index) * sin(harm_index * initial_frequency_new * time_data) + initial_frequency_new^2 * time_data * harm_index^3 * harmonic_parameters(harm_index, 2 * dof_index) .* cos(harm_index * initial_frequency_new * time_data));

    end
end

% Compute the residual considering the frequency sensitivity
residual(num_degrees_freedom + 2:2 * num_degrees_freedom+1,:) = mass_matrix * acceleration_freq_sensitivity + frac_damping_matrix * velocity_freq_sensitivity + stiffness_matrix * displacement_freq_sensitivity...
    + [-mu*Dx_alpha1_w0(1,:)+2*mu*displacement(1,:).*displacement_freq_sensitivity(1,:).*Dx_alpha1(1,:)+mu*displacement(1,:).^2.*Dx_alpha1_w0(1,:);...
    alpha*Dx_alpha2_w0(2,:)+2*beta*displacement(2,:).*displacement_freq_sensitivity(2,:)+3*epsilon*displacement(2,:).^2.*displacement_freq_sensitivity(2,:)];

velocity_freq_sensitivity_dof1 = zeros(1, length(time_data));
for dof_index = 1
    for harm_index = 1:num_harmonics
        % Calculate frequency sensitivity of velocity
        velocity_freq_sensitivity_dof1(dof_index, :) = velocity_freq_sensitivity_dof1(dof_index, :) +  harm_index  * harmonic_parameters(harm_index, 2 * dof_index);
    end
end
residual(2 * num_degrees_freedom+2, :) = velocity_freq_sensitivity_dof1;

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
            harmonic_term = harm_index * initial_frequency_new * time_data;
            displacement_harm_sensitivity(dof_index, :) = displacement_harm_sensitivity(dof_index, :) + harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);

            Y1_N_harm_alpha1=harm_index^alpha1*cos(beta2*pi/2);Y2_N_harm_alpha1=harm_index^alpha1*sin(beta2*pi/2);
            Dx_alpha1_a(dof_index,:)=Dx_alpha1_a(dof_index,:)+initial_frequency_new^alpha1*(Y1_N_harm_alpha1*harmonic_sensitivity(harm_index,2*dof_index)-Y2_N_harm_alpha1*harmonic_sensitivity(harm_index,2*dof_index-1))*cos(harmonic_term)+...
                initial_frequency_new^alpha1*(-Y1_N_harm_alpha1*harmonic_sensitivity(harm_index,2*dof_index-1)-Y2_N_harm_alpha1*harmonic_sensitivity(harm_index,2*dof_index))*sin(harmonic_term);

            Y1_N_harm_alpha2=harm_index^alpha2*cos(alpha2*pi/2);Y2_N_harm=harm_index^alpha2*sin(alpha2*pi/2);
            Dx_alpha2_a(dof_index,:)=Dx_alpha2_a(dof_index,:)+initial_frequency_new^alpha2*(Y1_N_harm_alpha2*harmonic_sensitivity(harm_index,2*dof_index)-Y2_N_harm*harmonic_sensitivity(harm_index,2*dof_index-1))*sin(harmonic_term)+...
                initial_frequency_new^alpha2*(Y1_N_harm_alpha2*harmonic_sensitivity(harm_index,2*dof_index-1)+Y2_N_harm*harmonic_sensitivity(harm_index,2*dof_index))*cos(harmonic_term);

            velocity_harm_sensitivity(dof_index, :) = velocity_harm_sensitivity(dof_index, :) - initial_frequency_new * harm_index * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * sin(harmonic_term) + initial_frequency_new * harm_index * harmonic_sensitivity(harm_index, 2 * dof_index) * cos(harmonic_term);
            acceleration_harm_sensitivity(dof_index, :) = acceleration_harm_sensitivity(dof_index, :) - (initial_frequency_new * harm_index).^2 * harmonic_sensitivity(harm_index, 2 * dof_index - 1) * cos(harmonic_term) - (initial_frequency_new * harm_index).^2 * harmonic_sensitivity(harm_index, 2 * dof_index) * sin(harmonic_term);
        end
    end

    residual((num_degrees_freedom+1) * (i + 1) + 1:(num_degrees_freedom+1) * (i + 1)+ 2, :) = mass_matrix * acceleration_harm_sensitivity + frac_damping_matrix * velocity_harm_sensitivity + stiffness_matrix * displacement_harm_sensitivity...
        + [-mu*Dx_alpha1_a(1,:)+2*mu*displacement(1,:).*displacement_harm_sensitivity(1,:).*Dx_alpha1(1,:)+mu*displacement(1,:).^2.*Dx_alpha1_a(1,:);...
        alpha*Dx_alpha2_a(2,:)+2*beta*displacement(2,:).*displacement_harm_sensitivity(2,:)+3*epsilon*displacement(2,:).^2.*displacement_harm_sensitivity(2,:)];

    velocity_harm_sensitivity_dof1 = zeros(1, length(time_data));
    for dof_index = 1
        for harm_index = 1:num_harmonics
            velocity_harm_sensitivity_dof1(dof_index, :) = velocity_harm_sensitivity_dof1(dof_index, :) + initial_frequency_new * harm_index * harmonic_sensitivity(harm_index, 2 * dof_index);
        end
    end

    residual((num_degrees_freedom+1) * (i + 2), :) = velocity_harm_sensitivity_dof1;
end

%% Compute sensitivities for constant terms
for dof = 1:num_degrees_freedom
    displacement_sensitivity_constant = zeros(num_degrees_freedom, length(time_data));
    velocity_sensitivity_constant = zeros(num_degrees_freedom, length(time_data));
    acceleration_sensitivity_constant = zeros(num_degrees_freedom, length(time_data));
    Dx_alpha1_sensitivity_constant = zeros(num_degrees_freedom, length(time_data));
    Dx_alpha2_sensitivity_constant = zeros(num_degrees_freedom, length(time_data));
    displacement_sensitivity_constant(dof, :) = ones(1, length(time_data));

    % Calculate residuals for constant sensitivities
    residual((num_degrees_freedom+1)*(2 * num_harmonics * num_degrees_freedom+1+dof)+1:(num_degrees_freedom+1)*(2 * num_harmonics * num_degrees_freedom+1+dof)+2,:) = mass_matrix * acceleration_sensitivity_constant + frac_damping_matrix * velocity_sensitivity_constant + stiffness_matrix * displacement_sensitivity_constant...
        + [-mu*Dx_alpha1_sensitivity_constant(1,:)+2*mu*displacement(1,:).*displacement_sensitivity_constant(1,:).*Dx_alpha1(1,:)+mu*displacement(1,:).^2.*Dx_alpha1_sensitivity_constant(1,:);...
        alpha*Dx_alpha2_sensitivity_constant(2,:)+2*beta*displacement(2,:).*displacement_sensitivity_constant(2,:)+3*epsilon*displacement(2,:).^2.*displacement_sensitivity_constant(2,:)];


    residual((num_degrees_freedom+1)*(2 * num_harmonics * num_degrees_freedom+2+dof), :) = velocity_sensitivity_constant(1,:);
end

residual = residual';
end


