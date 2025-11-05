function y=f(num,harmonic_params)
global num_harmonics time_data temp1 temp2 %N_harm index_global

initial_frequency = harmonic_params(1, 1);
harmonic_parameters = harmonic_params(3:end, :);

displacement = zeros(1, length(time_data));

% for 
dof_index = 1;
for harm_index = 1:num_harmonics
    harmonic_term = harm_index * initial_frequency * time_data;
    displacement(dof_index, :) = displacement(dof_index, :) + harmonic_parameters(harm_index, 2 * dof_index - 1) * cos(harmonic_term) + harmonic_parameters(harm_index, 2 * dof_index) * sin(harmonic_term);
end
displacement(dof_index, :) = displacement(dof_index, :) + harmonic_params(2, 2 * dof_index - 1);
% end

f1=displacement'-num(temp1:temp2,1);
% f2=velocity'-num(temp1:temp2,3:4);
% y=[f1(:,1);f1(:,2);f2(:,1);f2(:,2)];
y=f1;
end