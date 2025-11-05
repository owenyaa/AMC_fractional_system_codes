function nonlinear_terms = build_nonlinear_terms(angular_frequency, coefficient_matrix, increment_matrix)
%% Build nonlinear contributions for IHB (Residual / Ja / Jw) in MDOF form
% Outputs (contract kept compatible with original `nonterm`):
%   nonlinear_terms(1).part : nonlinear contribution to Residual   (R)
%   nonlinear_terms(2).part : nonlinear contribution to dR/d(a)    (Ja)
%   nonlinear_terms(3).part : nonlinear contribution to dR/d(ω)    (Jw)
%
% Model (van der Pol–Duffing type, per-DOF uncoupled nonlinearity):
%   x'' + vdp*(1 - param*x^2)*D^{alpha}(x) + k*x + k3*x^3 = 0
% Nonlinear pieces assembled here (consistent with previous pipeline):
%   R_nl  =  k3*x^3  -  vdp*param*ω^α * x^2 * D^{α}x
%   Ja_nl =  3*k3*x^2*Δx  -  vdp*param*ω^α * x^2 * D^{α}Δx  -  2*vdp*param*ω^α * x * D^{α}x * Δx
%   Jw_nl =  (-vdp) * ( x * x * D^{α}x )

global num_dofs cubic_stiffness_matrix fractional_order ...
       vdp_coefficient parameter_coefficient number_of_harmonics

H = number_of_harmonics;

% Split per DOF
for dof_index = 1:num_dofs
    disp_series_by_dof(dof_index).series = coefficient_matrix(:, 2*dof_index-1 : 2*dof_index); %#ok<AGROW>
    incr_series_by_dof(dof_index).series = increment_matrix(:,   2*dof_index-1 : 2*dof_index); %#ok<AGROW>
end

% Precompute fractional derivatives (series space)
for dof_index = 1:num_dofs
    frac_deriv_disp_by_dof(dof_index).series = ...
        fractional_derive_series_coeffs(disp_series_by_dof(dof_index).series, fractional_order); %#ok<AGROW>
    frac_deriv_incr_by_dof(dof_index).series = ...
        fractional_derive_series_coeffs(incr_series_by_dof(dof_index).series, fractional_order); %#ok<AGROW>
end

% Allocate accumulators (final size must be H x (2*num_dofs))
residual_nonlinear_series            = zeros(H, 2*num_dofs);
jacobian_wrt_coeffs_nonlinear_series = zeros(H, 2*num_dofs);
jacobian_wrt_frequency_nonlinear     = zeros(H, 2*num_dofs);

%% ------------------------ (1) Residual nonlinearity ------------------------
for dof_index = 1:num_dofs
    x_series   = disp_series_by_dof(dof_index).series;
    Dx_series  = frac_deriv_disp_by_dof(dof_index).series;

    % k3 * x^3
    blk1(1).matrix = cubic_stiffness_matrix * x_series;
    blk1(2).matrix = x_series;
    blk1(3).matrix = x_series;
    x_cubed = multiply_series_coeffs(blk1);
    x_cubed = x_cubed(1:H, :);  % trim

    % - vdp*param*ω^α * x^2 * D^{α}x
    blk2(1).matrix = -vdp_coefficient * parameter_coefficient * angular_frequency^fractional_order * x_series;
    blk2(2).matrix =  x_series;
    blk2(3).matrix =  Dx_series;
    x2_Dx = multiply_series_coeffs(blk2);
    x2_Dx = x2_Dx(1:H, :);      % trim

    sum_list(1).matrix = x_cubed;
    sum_list(2).matrix = x2_Dx;
    residual_for_dof   = sum_series_coeffs(sum_list);
    residual_for_dof   = residual_for_dof(1:H, :); % safety trim

    residual_nonlinear_series(:, 2*dof_index-1 : 2*dof_index) = residual_for_dof;
end

%% ------------------ (2) Ja nonlinearity (w.r.t. coefficients) --------------
for dof_index = 1:num_dofs
    x_series    = disp_series_by_dof(dof_index).series;
    Dx_series   = frac_deriv_disp_by_dof(dof_index).series;
    da_series   = incr_series_by_dof(dof_index).series;
    Dda_series  = frac_deriv_incr_by_dof(dof_index).series;

    % 3*k3*x^2*Δx
    jb1(1).matrix = 3 * cubic_stiffness_matrix * x_series;
    jb1(2).matrix =     x_series;
    jb1(3).matrix =     da_series;
    t1 = multiply_series_coeffs(jb1); t1 = t1(1:H, :);

    % - vdp*param*ω^α * x^2 * D^{α}Δx
    jb2(1).matrix = -vdp_coefficient * parameter_coefficient * angular_frequency^fractional_order * x_series;
    jb2(2).matrix =  x_series;
    jb2(3).matrix =  Dda_series;
    t2 = multiply_series_coeffs(jb2); t2 = t2(1:H, :);

    % - 2*vdp*param*ω^α * x * D^{α}x * Δx
    jb3(1).matrix = -2 * vdp_coefficient * parameter_coefficient * angular_frequency^fractional_order * x_series;
    jb3(2).matrix =      Dx_series;
    jb3(3).matrix =      da_series;
    t3 = multiply_series_coeffs(jb3); t3 = t3(1:H, :);

    s(1).matrix = t1; s(2).matrix = t2; s(3).matrix = t3;
    ja_for_dof  = sum_series_coeffs(s);
    ja_for_dof  = ja_for_dof(1:H, :);

    jacobian_wrt_coeffs_nonlinear_series(:, 2*dof_index-1 : 2*dof_index) = ja_for_dof;
end

%% --------------------- (3) Jw nonlinearity (w.r.t. ω) ----------------------
for dof_index = 1:num_dofs
    x_series   = disp_series_by_dof(dof_index).series;
    Dx_series  = frac_deriv_disp_by_dof(dof_index).series;

    fb(1).matrix = -vdp_coefficient * x_series;
    fb(2).matrix =  x_series;
    fb(3).matrix =  Dx_series;
    jw_for_dof   = multiply_series_coeffs(fb);
    jw_for_dof   = jw_for_dof(1:H, :);

    jacobian_wrt_frequency_nonlinear(:, 2*dof_index-1 : 2*dof_index) = jw_for_dof;
end

%% ----------------------------- Pack outputs --------------------------------
nonlinear_terms(1).part = residual_nonlinear_series;            % R (nonlinear)
nonlinear_terms(2).part = jacobian_wrt_coeffs_nonlinear_series; % Ja (nonlinear)
nonlinear_terms(3).part = jacobian_wrt_frequency_nonlinear;     % Jw (nonlinear)
end
