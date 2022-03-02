% Rishav Mallick, EOS, 2021
close all;
clearvars;
addpath functions/

% Specify parameters
M2KM = 1e-3;
H = 20e3; % Transition depth
T = 100; % Duration of earthquake cycle
eta = 3e18; % Viscosity-like parameter
n = 3; % Power-law exponent
x = linspace(-200e3, 200e3, 400);
t = 1; % Observation year

% Call to Mallick wrapper function
velocity_profile = mallick_power_law(x, t, H, eta, n, T);

% Plot single velocity profile
figure;
plot(M2KM * x, velocity_profile, "-r"); hold on;
xlabel("x (km)");
ylabel('v / v_0');
set(gca, "LineWidth", 1, "tickdir", "out")
