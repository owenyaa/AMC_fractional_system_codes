function big_jacobian = assemble_coefficients_jacobian(coefficient_matrix, angular_frequency)
%% Assemble Jacobian dR/da for self-excited periodic solution (IHB)
% Contract matches original `jacobi`:
%   Input:
%     coefficient_matrix : [H x (2*num_dofs)] current Fourier coeffs
%     angular_frequency  : ω
%   Output:
%     big_jacobian       : [(2H-1)*num_dofs  x  (2H-1)*num_dofs]
%
% Method
%   For each unit increment basis δa(i,j):
%     - Build derivatives δa_2 and δa_alpha in series space
%     - Linear operator pieces: M*ω^2*δa_2,  ε*ω^α*δa_alpha,  K*δa
%     - Nonlinear Ja contribution via build_nonlinear_terms(…, δa)
%     - Sum → trim to first H rows → pack per DOF → stack columns → place into S
%   Finally, rearrange S columns to form the square Jacobian (same as original).

global number_of_harmonics mass_matrix vdp_coefficient stiffness_matrix ...
       num_dofs fractional_order

H = number_of_harmonics;

% S will store column-wise assembled blocks before final rearrangement
S = zeros( (2*H-1)*num_dofs, 2*H*num_dofs );

for dof_col = 1:(2*num_dofs)            % j-index over columns of a
    for row_h = 1:H                      % i-index over harmonic rows
        % Unit increment basis δa
        delta_coeffs = zeros(H, 2*num_dofs);
        delta_coeffs(row_h, dof_col) = 1;

        % Derivatives in series space
        delta_ddot   = derive_series_coeffs(delta_coeffs, 2);
        delta_frac   = fractional_derive_series_coeffs(delta_coeffs, fractional_order);

        % Linear operator contributions: M*ω^2*δa_2 + ε*ω^α*δa_α + K*δa
        Ja_terms(1).matrix = apply_linear_operator_to_series(mass_matrix,      angular_frequency^2       * delta_ddot);
        Ja_terms(2).matrix = apply_linear_operator_to_series(vdp_coefficient,  angular_frequency^fractional_order * delta_frac);
        Ja_terms(3).matrix = apply_linear_operator_to_series(stiffness_matrix,                                 delta_coeffs);

        % Nonlinear contribution (Ja part) via build_nonlinear_terms
        nl = build_nonlinear_terms(angular_frequency, coefficient_matrix, delta_coeffs);
        Ja_terms(4).matrix = nl(2).part;

        % Sum, trim to first H rows, pack per DOF, then stack into a vector column
        summed = sum_series_coeffs(Ja_terms);
        summed = summed(1:H, :);
        packed = pack_series_coeffs_per_dof(summed);
        colvec = stack_columns_vector(packed);

        S(:, row_h + (dof_col-1)*H) = colvec;
    end
end

% Rearrangement to final square Jacobian (same as original indexing scheme)
big_jacobian = zeros( (2*H-1)*num_dofs, (2*H-1)*num_dofs );
for dof = 1:num_dofs
    % cos columns block
    big_jacobian(:, (dof-1)*(2*H-1)+1 : (dof-1)*(2*H-1)+H) = ...
        S(:, (dof-1)*(2*H) + 1 : (dof-1)*(2*H) + H);
    % sin columns block (skip the DC sine column inside S)
    big_jacobian(:, (dof-1)*(2*H-1) + H + 1 : (dof-1)*(2*H-1) + 2*H - 1) = ...
        S(:, (dof-1)*(2*H) + H + 2 : (dof-1)*(2*H) + 2*H);
end
end
