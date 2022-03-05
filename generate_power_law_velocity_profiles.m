% Generate and save velocity profiles for later use in
% basal slip estimation experiments
%
% n values = [1, 3, 6]
% A values = 1e18
% t values = [0.001, 5, 10, 25, 50, 100, 150, 200]

warning off;
close all;
clearvars;
addpath functions/

% Specify parameters
SECONDS_IN_A_YEAR = 60 * 60 * 24 * 365;
M2KM = 1e-3;
H = 20e3; % Transition depth
T = 200; % Duration of earthquake cycle
eta = 2e19; % Viscosity-like parameter
x = linspace(-200e3, 200e3, 40);

n_array = [1, 3, 6];
eta_array = [3e18, 3e19, 3e23];
t_array = [0.001, 1, 5, 10, 25, 50, 100, 150, 200] * SECONDS_IN_A_YEAR;

% n_array = [6];
% eta_array = [3e23];
% t_array = [0.001] * SECONDS_IN_A_YEAR;


n_store = {};
eta_store = {};
t_store = {};
velocity_profile_store = {};
for i = 1:numel(n_array)
    figure;
    hold on;
    eta = eta_array(i);
    for j = 1:numel(t_array)
        disp([i, j, n_array(i), t_array(j)]);
        n_store{i, j} = n_array(i);
        eta_store{i, j} = eta_array(i);
        t_store{i, j} = t_array(j);
        velocity_profile_store{i, j} = mallick_power_law(x, t_array(j), H, eta, n_array(i), T);

        plot(M2KM * x, velocity_profile_store{i, j}, "-r"); hold on;
        xlabel("x (km)");
        ylabel('v / v_0');
        set(gca, "LineWidth", 1, "tickdir", "out");
        drawnow;
    end
end

save("generate_power_law_velocity_profiles.mat");

% % Call to Mallick wrapper function
% velocity_profile = mallick_power_law(x, t, H, eta, n, T);
% 
% % Plot single velocity profile
% figure;
% plot(M2KM * x, velocity_profile, "-r"); hold on;
% xlabel("x (km)");
% ylabel('v / v_0');
% set(gca, "LineWidth", 1, "tickdir", "out")
