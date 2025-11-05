function colvec = stack_columns_vector(matrix_in)
%% Stack a 2D matrix column-wise into a single column vector
% Input : matrix_in [m x n]
% Output: colvec    [m*n x 1]  (column-major stacking)
[m, n] = size(matrix_in);
colvec = zeros(m*n, 1);
for j = 1:n
    colvec( (j-1)*m + 1 : j*m, 1 ) = matrix_in(:, j);
end
end
