function fractional_derivative_matrix = fractional_derive_series_coeffs(series_coeff_matrix, fractional_order)
%% Fractional derivative of Fourier-series coefficient matrix (series layout)
% Purpose
%   Compute D^{alpha}{x(t)} in the *coefficient space* for a trigonometric
%   Fourier expansion:
%       x(t) = sum_{k>=0} [ a_k * cos(k*.) + b_k * sin(k*.) ]
%   Here we operate directly on the coefficient "series matrix":
%       size(series_coeff_matrix) = [H  x  (2*num_dofs)]
%       per DOF: [cos_coeffs, sin_coeffs], row 1 is DC (k=0), rows 2..H are k=1..H-1.
%
% Inputs
%   series_coeff_matrix : [H x (2*num_dofs)] coefficients (DC row + harmonics)
%   fractional_order    : alpha in [0, 2]
%
% Output
%   fractional_derivative_matrix : same size as input, coefficients of D^{alpha}x
%
% Spectral rule used (kept consistent with your original codebase):
%   For harmonic k>=1:
%     c1 = k^alpha * cos(alpha*pi/2)
%     c2 = k^alpha * sin(alpha*pi/2)
%   Then, per DOF:
%     D^{alpha}{cos} ↔  c1 * cos  ±  c2 * sin
%     D^{alpha}{sin} ↔  c1 * sin  ∓  c2 * cos
%   DC (k=0) sine column remains 0.

% ----------------------------- argument checks -----------------------------
if fractional_order < 0 || fractional_order > 2
    warning('fractional_derive_series_coeffs:alphaOutOfRange', ...
        'fractional_order (alpha) should be within [0, 2].');
end

[num_harmonics, num_columns] = size(series_coeff_matrix);
fractional_derivative_matrix = zeros(num_harmonics, num_columns);

% ---------------------------- main (0 < alpha <= 1) ------------------------
% NOTE: Your original code used "if alpha>0 || alpha<=1" which is always true
% (for alpha in [0,2]) and made the else-branch unreachable. We keep the *same*
% behavior here (i.e., always execute this branch) to avoid any functional change.
for harmonic_idx = 1:num_harmonics
    k = harmonic_idx - 1;                 % harmonic order (k=0 is DC)
    c1 = k^fractional_order * cos(fractional_order*pi/2);
    c2 = k^fractional_order * sin(fractional_order*pi/2);

    for col_idx = 1:num_columns
        if mod(col_idx, 2) == 0
            % even column → sine coefficients
            % D^{alpha}{sin} = c1 * sin  -  c2 * cos
            fractional_derivative_matrix(harmonic_idx, col_idx) = ...
                c1 * series_coeff_matrix(harmonic_idx, col_idx) ...
              - c2 * series_coeff_matrix(harmonic_idx, col_idx - 1);
        else
            % odd column  → cosine coefficients
            % D^{alpha}{cos} = c1 * cos  +  c2 * sin
            fractional_derivative_matrix(harmonic_idx, col_idx) = ...
                c1 * series_coeff_matrix(harmonic_idx, col_idx) ...
              + c2 * series_coeff_matrix(harmonic_idx, col_idx + 1);
        end
    end
end

% ------------------------------ DC sine row --------------------------------
% Ensure DC (k=0) sine columns are exactly zero
for dof = 1:(num_columns/2)
    fractional_derivative_matrix(1, 2*dof) = 0;
end

% ----------------------------- historical note -----------------------------
% The original function contained an "else" branch with different algebra and
% references to undefined symbols (alpha1, Y1_N_harm_alpha1, etc.). Because the
% original condition was `if alpha>0 || alpha<=1`, that else-block was never
% reachable for alpha in [0,2]. To preserve behavior and avoid introducing bugs,
% we have *not* re-enabled that branch here. If you later need a distinct rule
% for alpha in (1,2], we can add it explicitly and update all callers together.

end
