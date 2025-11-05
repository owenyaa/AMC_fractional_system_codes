function derivative_matrix = derive_series_coeffs(series_coeff_matrix, derivative_order)
%% Derive Fourier-series coefficient matrix for trigonometric terms
% Purpose:
%   Compute the n-th derivative of a Fourier-series coefficient matrix
%   with respect to time, for a representation:
%       x(t) = Σ [a_i*cos(iωt) + b_i*sin(iωt)]
%
% Inputs:
%   series_coeff_matrix : [H x (2*num_dofs)] matrix of Fourier coefficients
%                         Each DOF has two columns [cos_coeffs, sin_coeffs].
%                         Row 1 is DC term; rows 2..H correspond to harmonics.
%   derivative_order    : order of differentiation (integer ≥ 0)
%
% Outputs:
%   derivative_matrix   : coefficient matrix after differentiation
%                         Same dimension as input (H x 2*num_dofs)
%
% Method:
%   - Differentiation in time introduces factors of (iω)^n.
%   - For odd derivatives: cos→±sin, sin→∓cos.
%   - For even derivatives: cos→±cos, sin→±sin.
%   - DC term (row 1) has zero derivative for n ≥ 1.

% -------------------------------------------------------------------------

[num_harmonics, num_columns] = size(series_coeff_matrix);
derivative_matrix = zeros(num_harmonics, num_columns);

is_odd_derivative = mod(derivative_order, 2);  % 1=odd, 0=even

for harmonic_index = 1:num_harmonics
    for col_index = 1:num_columns
        % DC term has zero derivative when order >= 1
        if harmonic_index == 1 && derivative_order >= 1
            derivative_matrix(harmonic_index, col_index) = 0;
            continue;
        end

        harmonic_number = harmonic_index - 1; % actual harmonic order

        if is_odd_derivative
            % Odd derivative → cos↔sin with sign alternation
            % (-1) term switches between + and - depending on parity of column
            sign_factor = (-1)^((derivative_order - 3 + 2*(2 - mod(col_index,2))) / 2);
            paired_col = col_index + (-1)^(mod(col_index,2) + 1);
            derivative_matrix(harmonic_index, col_index) = ...
                (harmonic_number)^derivative_order * sign_factor * ...
                series_coeff_matrix(harmonic_index, paired_col);
        else
            % Even derivative → same function type, possible sign inversion
            sign_factor = (-1)^(derivative_order / 2);
            derivative_matrix(harmonic_index, col_index) = ...
                (harmonic_number)^derivative_order * sign_factor * ...
                series_coeff_matrix(harmonic_index, col_index);
        end
    end
end

% Ensure DC component of sine coefficients is zero (cos derivative vanishes)
for dof = 1:(num_columns/2)
    derivative_matrix(1, 2*dof) = 0;
end

end
