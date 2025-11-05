function packed_cols = pack_series_coeffs_per_dof(series_matrix)
%% Pack per-DOF Fourier series coefficients into stacked columns
% Purpose
%   Convert [H x (2*num_dofs)] series matrix (per DOF: [cos, sin])
%   into [2H-1 x num_dofs] stacked columns:
%     column d = [cos_coeffs(:,d); sin_coeffs(2:end,d)]
%
% Input
%   series_matrix : [H x (2*num_dofs)] with [cos, sin] for each DOF
%
% Output
%   packed_cols   : [2*H-1 x num_dofs], per-DOF stacked column

[num_harmonics, two_ndofs] = size(series_matrix);
num_dofs = two_ndofs / 2;

packed_cols = zeros(2*num_harmonics - 1, num_dofs);
for dof = 1:num_dofs
    cos_col = series_matrix(:, 2*dof - 1);
    sin_col = series_matrix(:, 2*dof    );  % sin DC = row1 (usually 0 by convention)
    packed_cols(1:num_harmonics, dof)               = cos_col;
    packed_cols(num_harmonics+1:end, dof)           = sin_col(2:end);
end
end
