function product_series = multiply_series_coeffs(series_operands)
%% Multiply multiple Fourier series (in coefficient-matrix form) via convolution
% Input
%   series_operands : struct array with field .matrix
%                     each .matrix is a series coefficient matrix:
%                       size = [H x (2*num_dofs)], columns per DOF = [cos, sin]
% Output
%   product_series  : resulting series coefficient matrix after sequential multiplications
%
% Behavior
%   ((...((A ⊗ B) ⊗ C) ⊗ D) ... ), where ⊗ denotes spectral convolution
%   implemented by harmonic_product_convolution.

num_operands = numel(series_operands);
product_series = series_operands(1).matrix;
for k = 2:num_operands
    product_series = harmonic_product_convolution(product_series, series_operands(k).matrix);
end
end
