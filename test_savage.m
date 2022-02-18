close all;
clearvars;

% Parameters for Savage (2000) model
mu = 3e10;
eta = 5e18;
n_pts = 100;
x = linspace(-200, 200, n_pts);
y = zeros(size(x));
T = 100;
D = 20;

% Velocities early in the earthquake cycle
t = 1;
v_early = savage_2000(x, y, t, D, mu, eta, T);

% Velocities late in the earthquake cycle
t = 99;
v_late = savage_2000(x, y, t, D, mu, eta, T);

% Plot velocities relative to steady state
figure("color", "w");
hold on;
plot(x/D, 1/pi * atan(x./D), "-k");
plot(x/D, v_early, "-r");
plot(x/D, v_late, "--r");
ylabel("v / v_0");
xlabel("x / D");
legend_handle = legend("steady state", "early", "late");
set(legend_handle, "Location", "northwest");
box on;
set(gca, "TickDir", "out");
set(gca, "fontsize", 18);
ytickformat("%0.1f")