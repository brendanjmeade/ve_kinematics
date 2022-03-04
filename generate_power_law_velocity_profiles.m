% Generate and save velocity profiles for later use in
% basal slip estimation experiments
%
% n values = [1, 3, 6]
% A values = 1e18
% t values = [0.001, 5, 10, 25, 50, 100, 150, 200]

close all;
clearvars;
addpath functions/

% Specify parameters
M2KM = 1e-3;
H = 20e3; % Transition depth
T = 200; % Duration of earthquake cycle
eta = 3e18; % Viscosity-like parameter
x = linspace(-200e3, 200e3, 400);

n_array = [1, 3, 6];
t_array = [0.001, 1, 5, 10, 25, 50, 100, 150, 200];
n_store = {};
t_store = {};
velocity_profile_store = {};
tic;
for i = 1:numel(n_array)
    for j = 1:numel(t_array)
        disp([i, j]);
        n_store{i, j} = n_array(i);
        t_store{i, j} = n_array(i);
        velocity_profile_store{i, j} = mallick_power_law(x, t_array(i), H, eta, n_array(i), T);
    end
end
toc;

tic;
save("generate_power_law_velocity_profiles.mat");
toc;

% % Call to Mallick wrapper function
% velocity_profile = mallick_power_law(x, t, H, eta, n, T);
% 
% % Plot single velocity profile
% figure;
% plot(M2KM * x, velocity_profile, "-r"); hold on;
% xlabel("x (km)");
% ylabel('v / v_0');
% set(gca, "LineWidth", 1, "tickdir", "out")
