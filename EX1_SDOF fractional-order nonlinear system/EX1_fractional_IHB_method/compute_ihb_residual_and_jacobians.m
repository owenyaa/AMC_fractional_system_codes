function ihb_results = compute_ihb_residual_and_jacobians(coefficient_matrix, angular_frequency)
%% Compute residual R, Ja=dR/da, and Jw=dR/dω for the IHB linearization

global number_of_harmonics mass_matrix vdp_coefficient stiffness_matrix ...
       num_dofs fractional_order

H = number_of_harmonics;

% Zero increment used by the nonlinear builder for R term (interface stable)
zero_increment_coeffs = zeros(H, 2*num_dofs);

% Spectral derivatives (series-coefficient layout)
first_deriv_coeffs      = derive_series_coeffs(coefficient_matrix, 1); %#ok<NASGU>
second_deriv_coeffs     = derive_series_coeffs(coefficient_matrix, 2);
fractional_deriv_coeffs = fractional_derive_series_coeffs(coefficient_matrix, fractional_order);

% Nonlinear contributions (R / Ja / Jw)
nonlinear_terms = build_nonlinear_terms(angular_frequency, coefficient_matrix, zero_increment_coeffs);

%% ------------------------------- Residual R --------------------------------
% residual_series = -[ M*(ω^2 x'') + K*x + vdp*(ω^α D^α x) ] - nonlinear_terms(1)
rblk(1).matrix = apply_linear_operator_to_series(mass_matrix,     -angular_frequency^2          * second_deriv_coeffs);
rblk(2).matrix = apply_linear_operator_to_series(vdp_coefficient, -angular_frequency^fractional_order * fractional_deriv_coeffs);
rblk(3).matrix = apply_linear_operator_to_series(stiffness_matrix, -coefficient_matrix);
rblk(4).matrix = -nonlinear_terms(1).part;

residual_series = sum_series_coeffs(rblk);
residual_series = residual_series(1:H, :);  % trim to H rows

residual_packed = pack_series_coeffs_per_dof(residual_series);
residual_vector = stack_columns_vector(residual_packed);

%% ---------------------------- Ja: dR/d(coeffs) -----------------------------
jacobian_wrt_coeffs = assemble_coefficients_jacobian(coefficient_matrix, angular_frequency);

%% ----------------------------- Jw: dR/dω -----------------------------------
% d/dω[ M*(ω^2 x'') ] = 2ω M x''
% d/dω[ vdp*(ω^α D^α x) ]  →  vdp * D^α x  (kept consistent with original pipeline)
wblk(1).matrix = apply_linear_operator_to_series(mass_matrix,     2*angular_frequency * second_deriv_coeffs);
wblk(2).matrix = apply_linear_operator_to_series(vdp_coefficient, fractional_deriv_coeffs);
wblk(3).matrix = nonlinear_terms(3).part;

freq_jac_series = sum_series_coeffs(wblk);
freq_jac_series = freq_jac_series(1:H, :);  % trim

freq_packed            = pack_series_coeffs_per_dof(freq_jac_series);
jacobian_wrt_frequency = stack_columns_vector(freq_packed);

%% ------------------------------ Pack outputs -------------------------------
ihb_results(1).vector = residual_vector;          % R
ihb_results(2).vector = jacobian_wrt_coeffs;      % Ja
ihb_results(3).vector = jacobian_wrt_frequency;   % Jw
end
