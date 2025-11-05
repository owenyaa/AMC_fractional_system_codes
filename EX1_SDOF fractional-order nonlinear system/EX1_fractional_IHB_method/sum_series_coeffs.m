function summed_series = sum_series_coeffs(series_list)
%% Sum multiple Fourier series coefficient matrices (cos/sin paired by DOF)
% Input
%   series_list : struct array with field .matrix
%                 each .matrix is [H x (2*num_dofs)] with per-DOF columns [cos, sin]
% Output
%   summed_series : the sum over all .matrix, row-aligned (zero-padded as needed)
%
% Notes
%   - Handles different row counts by expanding to the max row length.
%   - Preserves the convention: DC sine (row 1, even columns) set to 0.

num_terms = numel(series_list);
max_rows  = 0;
num_cols  = size(series_list(1).matrix, 2);

% find max number of rows among inputs
for k = 1:num_terms
    max_rows = max(max_rows, size(series_list(k).matrix, 1));
end

summed_series = zeros(max_rows, num_cols);

% accumulate using the two-matrix adder (handles different heights)
for k = 1:num_terms
    summed_series = add_series_coeff_mats(summed_series, series_list(k).matrix);
end
end
