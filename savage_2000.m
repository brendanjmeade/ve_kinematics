function v = savage_2000(x, y, t, H, mu, eta, T)
% This calculates the velocity at a point (x, y) in time from Savage's
% viscoelastic earthquake cycle result from 2000
% (v. 105, no. B11, p. 25,525-25-532)
%
% Arguments:
%   x           : x position (horizontal distance from fault) [km]
%   y           : y position (depth) [km]
%   t           : time since last earthquake [years]
%   H           : thickness of elastic layer [km]
%   mu          : shear modulus [Pa]
%   eta         : dynamic viscocity [Pa s]
%   T           : length of earthquake cycle [years]
%
% Returned variables:
%   v           : velocity [dimensionless]


% Convert both t and T from years to seconds
seconds_in_a_year = 60 * 60 * 24 * 365.25;
t = seconds_in_a_year * t;
T = seconds_in_a_year * T;

% Calculate relaxation time
tau_0 = (mu * T) / (2 * eta);
tau = (mu * t) / (2 * eta);

% Set resolution, time and counting parameters
max_m = 50;
max_k = 250 / tau_0;

% Calculate the sum just over S_m (vectorized in x, not y)
if (y < H)
   sum_S_m = 1 / (2 * pi) * (-atan((H + y) ./ x) - atan((H - y) ./ x) + pi .* sign(x));
else
   sum_S_m = 1 / (2 * pi) * (-atan((H + y) ./ x) - atan((y - H) ./ x) + pi .* sign(x));
end

% Calculate S_m as a function of x, y, and m
for m = 1 : max_m
   if (y < H)
      S_m(m, :) = 1 / (2 * pi) * (atan(((2 * m + 1) * H + y) ./ x) - atan(((2 * m - 1) * H + y) ./ x) + atan(((2 * m + 1) * H - y) ./ x) - atan(((2 * m - 1) * H - y) ./ x));
   else
      S_m(m, :) = 1 / (2 * pi) * (atan(((2 * m + 1) * H + y) ./ x) - atan(((2 * m - 3) * H + y) ./ x));
   end
end

% Calculate C_m as a function of tau, tau_0, and m
for m = 1 : max_m
   sum_k = 0;
   for k = 0 : max_k
      sum_k = sum_k + exp(-(tau + k * tau_0)) * (tau + k * tau_0)^(m - 1);
   end
   C_m(m) = tau_0 / factorial(m - 1) * sum_k + gammainc(m, tau + (max_k + 1) * tau_0) / factorial(m - 1);
end

% Calculate velocities current time
v = zeros(size(x));
for m = 1 : max_m
   v = v + (C_m(m) - 1) .* S_m(m, :);
end
v = v + sum_S_m;



