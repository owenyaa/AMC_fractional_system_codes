function series_out = apply_linear_operator_to_series(operator_matrix, series_in)
%% Apply a linear operator (scalar or square matrix) to a series coefficient matrix
% Inputs
%   operator_matrix : scalar or [num_dofs x num_dofs] matrix (e.g., mass, stiffness)
%   series_in       : [H x (2*num_dofs)] series coeffs, columns per DOF = [cos, sin]
% Output
%   series_out      : same size as series_in; operator applied DOF-wise
%
% Notes
%   For each DOF block (two columns: cos/sin), we form a linear combination
%   with operator rows. This preserves the [cos, sin] pairing per DOF.

[num_harmonics, total_cols] = size(series_in);
num_dofs = total_cols / 2;

% Split per DOF into 2-column blocks
for d = 1:num_dofs
    dof_block(d).part = series_in(:, 2*d-1:2*d); %#ok<AGROW>
end

% First column combination
for i = 1:num_dofs
    series_blocks(i).part = operator_matrix(i,1) * dof_block(1).part; %#ok<AGROW>
end

% Accumulate remaining columns
for i = 1:num_dofs
    for j = 1:(num_dofs - 1)
        series_blocks(i).part = add_series_coeff_mats( ...
            series_blocks(i).part, operator_matrix(i, j+1) * dof_block(j+1).part);
    end
end

% Stitch back to full series matrix
series_out = zeros(num_harmonics, 2*num_dofs);
for d = 1:num_dofs
    series_out(:, 2*d-1:2*d) = series_blocks(d).part;
end
end
