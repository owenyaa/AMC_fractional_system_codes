function [Predictor, ds_new] = arc_length_new(every_a, Nd, ds_old, iteration_index)
% ARC_LENGTH_NEW Predicts the next solution in a nonlinear system using the arc-length increment method with dynamic step adjustment.
%
% This function implements the arc-length increment method, which is used to predict the next solution point in a nonlinear system
% based on previous solution points. It dynamically adjusts the arc-length increment to ensure stability and accuracy in the prediction.
%
% Inputs:
%   every_a        - Structure array containing previous solutions (x0, x1, x2, x3).
%                    Each structure should have a field 'harmonic_coefficients' representing the solution at that point.
%   Nd             - Desired number of iterations or a control parameter related to the response curve (commonly 4 or 5).
%   ds_old         - Arc length increment from the previous step, representing the step size used in the last prediction.
%   iteration_index - Current iteration count from the main program, indicating how many iterations have been performed so far.
%
% Outputs:
%   Predictor      - Predicted next solution (x4), obtained through cubic interpolation based on previous solutions.
%   ds_new         - Updated arc length increment for the next step, dynamically adjusted based on the iteration count.

%% Step 1: Calculate arc lengths between consecutive known solution points
% Initialize an array to store the arc lengths between the four known points.
s = zeros(1, 3); % Arc lengths between x0-x1, x1-x2, and x2-x3

for j = 1:3
    % Calculate the Euclidean distance (norm) between consecutive harmonic coefficients.
    % This represents the arc length between each pair of consecutive solution points.
    s(j) = sqrt(norm(every_a(j+1).harmonic_coefficients - every_a(j).harmonic_coefficients));

    % To prevent numerical issues, ensure that the arc length is not effectively zero.
    if s(j) < 1e-15
        s(j) = 1e-15; % Set a minimal arc length threshold
    end
end

%% Step 2: Define time parameters 't' based on the calculated arc lengths
% The 't' array represents the parameterization of the solution path based on arc length.
% t(1) to t(4) correspond to the cumulative arc lengths up to each known solution point.
t = zeros(5, 1); % Initialize the time parameter array for four known points and one predicted point
t(1) = 0;         % Start parameter at 0 for the first solution point (x0)
t(2) = s(1);      % Parameter for the second point (x1) is the first arc length
t(3) = t(2) + s(2); % Parameter for the third point (x2) is cumulative arc lengths s1 + s2
t(4) = t(3) + s(3); % Parameter for the fourth point (x3) is cumulative arc lengths s1 + s2 + s3

%% Step 3: Adaptive arc length increment for predicting the next solution (x4)
% The arc length increment 'ds_new' is dynamically adjusted based on the desired number of iterations (Nd)
% and the current iteration count. This adjustment ensures that the step size adapts as the algorithm progresses,
% potentially improving convergence and stability.
ds_new = (ds_old * Nd) / iteration_index; % Update the arc length increment based on iteration progress
t(5) = t(4) + ds_new; % Parameter for the predicted point (x4) is the last parameter plus the new increment

%% Step 4: Perform cubic (Lagrange) interpolation to predict the next solution
% Initialize the Predictor with zeros, matching the size of the harmonic coefficients.
Predictor = zeros(size(every_a(1).harmonic_coefficients)); % Initialize Predictor as a zero vector

% Loop over the four known solution points to construct the Lagrange interpolating polynomial.
for i = 1:4
    L = 1; % Initialize the Lagrange basis polynomial for the i-th point

    for j = 1:4
        if i ~= j
            % Construct the Lagrange basis polynomial by multiplying terms of the form:
            % (t5 - tj) / (ti - tj), where ti and tj are the parameters of the i-th and j-th points respectively.
            L = L * (t(5) - t(j)) / (t(i) - t(j));
        end
    end

    % Accumulate the contribution of each known point's harmonic coefficients weighted by the Lagrange basis.
    Predictor = Predictor + L * every_a(i).harmonic_coefficients;
end

% The Predictor now holds the predicted harmonic coefficients for the next solution point (x4).
end
