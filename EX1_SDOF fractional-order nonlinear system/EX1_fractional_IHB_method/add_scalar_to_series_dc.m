function series_out = add_scalar_to_series_dc(scalar_value, series_in)
%% Add a scalar to the DC cosine entry of a series coefficient matrix
% Inputs
%   scalar_value : scalar to add
%   series_in    : [H x (2*num_dofs)] series coeffs (DC row = row 1)
% Output
%   series_out   : same as series_in, with (1,1) increased by scalar_value

series_in(1,1) = series_in(1,1) + scalar_value;
series_out = series_in;
end
