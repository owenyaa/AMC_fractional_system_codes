function series_matrix = columns_to_series_matrix(stacked_by_dof)
%% Rebuild series coefficient matrix from DOF-stacked column layout
% Purpose
%   Convert the stacked-by-DOF layout (each DOF is one column made by
%   stacking its cosine series followed by sine series-without-DC)
%   back to the standard series matrix layout:
%       [H x (2*num_dofs)] with per-DOF columns [cos, sin].
%
% Input
%   stacked_by_dof : [ (ceil(H/?) * 2 - 1) x num_dofs ]  layout produced by your solver:
%                    For each DOF column:
%                      - Top L rows are cosine coefficients (including DC).
%                      - Next (L-1) rows are sine coefficients (starting at k=1).
% Output
%   series_matrix  : [L x (2*num_dofs)] where
%                    - column 2*d-1 : cosine coefficients for DOF d
%                    - column 2*d   : sine   coefficients for DOF d (row 1 = 0)
%
% Notes
%   - This is the exact inverse of your packing convention where sine DC is omitted
%     and implicitly considered zero.

[num_rows_stacked, num_dofs] = size(stacked_by_dof);

% The input format uses L cosine rows and (L-1) sine rows stacked,
% so total rows per DOF column = L + (L-1) = 2L - 1  =>  L = ceil(num_rows_stacked/2)
num_harmonics = ceil(num_rows_stacked / 2);

series_matrix = zeros(num_harmonics, 2*num_dofs);

for d = 1:num_dofs
    % cosine coefficients (including DC at row 1)
    series_matrix(:, 2*d-1) = stacked_by_dof(1:num_harmonics, d);

    % sine coefficients start from harmonic 1 (row 2..end)
    if num_harmonics > 1
        series_matrix(2:end, 2*d) = stacked_by_dof(num_harmonics+1:end, d);
    end
    % row 1, sine column stays 0 by construction
end
end
