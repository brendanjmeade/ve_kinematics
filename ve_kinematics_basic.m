close all;
clear variables;

% Reference parameters
siay = 1;
n_pts = 1000;
x = linspace(-200e3, 200e3, n_pts);
x_buffer_fraction = 0.1;
y = zeros(size(x));
t = 50 * siay;
H = 15e3;
mu = 30e10;
eta = 5e19;
T = 100 * siay;
fontsize = 18;

% Calculate and plot viscoelastic velocity profile at one time
n_ve_profiles = 5;
t_vec = linspace(0, T, n_ve_profiles);
v_viscoelastic = zeros(n_ve_profiles, n_pts);
for i = 1:length(t_vec)
    v_viscoelastic(i, :) = savage_2000(x, y, t_vec(i), H, mu, eta, T);
end
v_steady_state = 1/pi .* atan(x./H);
v_steady_state = v_steady_state(:);

figure;
hold on;
ve_handle = plot(x / 1e3, v_viscoelastic, '-r', LineWidth=2.0);
steady_state_handle = plot(x / 1e3, v_steady_state, '-b', LineWidth=2.0);
xlabel("x (km)");
ylabel("v (mm/yr)");
set(gca, "fontsize", fontsize);
set(gca, "Tickdir", "out");
set(gcf, "Color", "white")
box on;
set(gca, "XTick", [-200, -100, 0, 100, 200]);
set(gca, "YTick", [-2.0, -1.0, 0.0, 1.0, 2.0]);
set(gca, "YTickLabel", {"-2.0", "-1.0", "0.0", "1.0", "2.0"});
legend([ve_handle(1), steady_state_handle], "viscoelastic", "steady state")
title("steady state and viscoelastic velocity profiles", "fontsize", fontsize, FontWeight="normal")

% Define patches for basal detachment surface
F.xmin = -220;
F.xmax = 220;
F.nfaults = 20;
F.fault_width = (F.xmax - F.xmin) / F.nfaults;
F.xf = F.xmin+F.fault_width:F.fault_width:F.xmax;
F.length = 20000; % along strike length
F.yf = -F.length / 2; % anchor point y
F.strike = deg2rad(90);
F.depth = 15; % anchor point depth
F.dip = 0; % dip
F.poissonsratio = 0.25;


% Calculate partials for each of the basal fault patches
G = zeros(n_pts, F.nfaults);
for i = 1:F.nfaults
   [~, uy, ~] = okada_1985(F.xf(i), F.yf, F.strike, F.depth, F.dip, F.length, ...
      F.fault_width, 1, 0, 0, x/1e3, y, F.poissonsratio);
   G(:, i) = uy;
end


% Plot partials for each of the basal fault patches
figure(Position=[0, 0 1000, 600]);
hold on;
tiledlayout(5, 5, TileSpacing="none");
for i = 1:F.nfaults
    nexttile;
    hold on;

    % Draw patch width
    fh = fill([F.xf(i) - F.fault_width, F.xf(i), F.xf(i), F.xf(i) - F.fault_width], [0.0, 0.0, 0.5, 0.5], "y");
    fh.FaceColor = 0.85 * [1, 1, 1];
    fh.EdgeColor = "None";

    % Draw velocity profile
    plot(x / 1e3, G(:, i), '-r', LineWidth=2);
    text(-180, 0.43, string(i), "fontsize", fontsize)
    set(gca, "xlim", [-200, 200])
    set(gca, "ylim", [0, 0.5])
    box on;

    set(gca, "Tickdir", "out");
    set(gca, "XTick", []);
    set(gca, "YTick", []);
    if i == 16
        set(gca, "XTick", [-200, 0, 200]);
        set(gca, "YTick", [0.0, 0.25, 0.5]);
        set(gca, "YTickLabel", {"0.00", "0.25", "0.50"});
        xlabel("x (km)");
        ylabel("v (mm/yr)");
        set(gca, "fontsize", fontsize);
    end
end
% set(gca, "fontsize", fontsize);
set(gcf, "Color", "white")
sgtitle("velocities from unit slip on 20 basal dislocation patches", "fontsize", fontsize);


% Solve for equivalent slip disribution
v = v_viscoelastic(1, :);
v = v(:);
effectiveslip = G \ v;
if size(G, 1) > size(G, 2)
   effectiveslip = inv(G' * G) * G' * v;
elseif size(G, 2) > size(G, 1)
   effectiveslip = G' * inv(G * G') * v;
end
basal_model_velocities = G * effectiveslip;

% Plot comparision of direct VE calculation and basal approximation
markersize = 15;
figure;
set(gcf, "Color", "w");
hold on;
plot(x / 1e3, v, '-b', 'markersize', markersize, LineWidth=2.0);
plot(x(1:20:end) / 1e3, basal_model_velocities(1:20:end), '.r', 'markersize', markersize);
set(gca, "XTick", [-200, -100, 0, 100, 200]);
set(gca, "YTick", [-2.0, -1.0, 0.0, 1.0, 2.0]);
set(gca, "YTickLabel", {"-2.0", "-1.0", "0.0", "1.0", "2.0"});
xlabel("x (km)");
ylabel("v (mm/yr)");
set(gca, "fontsize", fontsize);
set(gca, "Tickdir", "out")
box on;
legend("Maxwell viscoelastic model", "basal dislocations representation", Location="northwest");

return;

figure;
plot(effectiveslip, '-k+');



