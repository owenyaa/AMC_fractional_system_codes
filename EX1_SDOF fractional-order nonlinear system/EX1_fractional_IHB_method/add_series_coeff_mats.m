function out_series = add_series_coeff_mats(series_A, series_B)
%% Add two Fourier series coefficient matrices with row alignment
% Inputs
%   series_A, series_B : [H_A x (2*num_dofs)], [H_B x (2*num_dofs)]
% Output
%   out_series         : [max(H_A,H_B) x (2*num_dofs)] element-wise sum
%
% Behavior
%   - Pads the shorter matrix in rows.
%   - Keeps per-DOF [cos, sin] layout.
%   - Enforces DC sine to be zero in the output (row 1, even columns).

rows_A  = size(series_A, 1);
rows_B  = size(series_B, 1);
num_dofs = size(series_A, 2) / 2;

rows_out   = max(rows_A, rows_B);
out_series = zeros(rows_out, 2*num_dofs);

% ensure A is the taller (swap if needed)
if rows_A < rows_B
    tmp      = series_B;
    series_B = series_A;
    series_A = tmp;
    rows_B   = rows_A;
    rows_A   = rows_out;
end

% rows where both have data
if rows_B > 0
    out_series(1:rows_B, :) = series_A(1:rows_B, :) + series_B(1:rows_B, :);
end

% remaining rows from the taller A
if rows_A > rows_B
    out_series(rows_B+1:rows_A, :) = series_A(rows_B+1:rows_A, :);
end

% enforce DC sine (row 1, even columns) = 0
for d = 1:num_dofs
    out_series(1, 2*d) = 0;
end
end
