function s = sign_nonzero(x)
%% Sign function with sign_nonzero(0) = 0
% Returns:
%   -1 if x < 0,  0 if x == 0,  +1 if x > 0

if x == 0
    s = 0;
else
    s = x / abs(x);
end
end
