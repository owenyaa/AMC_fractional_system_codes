function matrix_out = unstack_columns_by_dof(vector_in)
%% Unstack a vector back to [n_rows x num_dofs] by DOF columns
% Purpose
%   Inverse of stack_columns_vector (for DOF-wise column matrices).
% Inputs
%   vector_in  : [(n_rows*num_dofs) x 1]
% Globals
%   num_dofs   : number of DOFs
% Output
%   matrix_out : [n_rows x num_dofs]

global num_dofs
L  = numel(vector_in);
nr = L / num_dofs;

matrix_out = zeros(nr, num_dofs);
for dof = 1:num_dofs
    matrix_out(:, dof) = vector_in( 1 + (dof-1)*nr : dof*nr );
end
end
