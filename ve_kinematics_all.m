close all;
clear variables;

addpath(genpath(pwd));

% Load precomputed velocity profiles
load generate_power_law_velocity_profiles.mat

% Reference parameters
fontsize = 18;
SAVE_FIGURES = false;

% Calculate and plot viscoelastic velocity profile at one time
n_ve_profiles = numel(t_array);
v_steady_state = 1/pi .* atan(x./H);
velocity_perturbation_profile_store = {};
for i = 1:3
    for j = 1:n_ve_profiles
        velocity_perturbation_profile_store{i, j} = velocity_profile_store{i, j} - v_steady_state;
    end
end

% Single idealized plot of 
% 1 - Total
% 2 - Block motion
% 3 - Coseismic slip deficit
% 4 - postseismic perturbation
v_block = zeros(size(x));
v_block(x < 0) = -0.5;
v_block(x >= 0) = 0.5;
v_csd = -1/pi .* atan(H ./ x);
ve_index = 6;
figure("Position", [0, 0, 600, 400]);
hold on;
line_width = 2.0;
ve_handle = plot(x / 1e3, velocity_profile_store{1, ve_index}, "-k", LineWidth=line_width);
block_handle = plot(x / 1e3, v_block, ".-r", LineWidth=line_width);
csd_handle = plot(x(x<0) / 1e3, v_csd(x<0), ":r", LineWidth=line_width);
csd_handle = plot(x(x>0) / 1e3, v_csd(x>0), ":r", LineWidth=line_width);
perturbation_handle = plot(x / 1e3, velocity_perturbation_profile_store{1, ve_index}, "--r", LineWidth=line_width);
xlabel("$x \; \mathrm{(km)}$", "interpreter", "latex");
ylabel("$v \; \mathrm{(mm/yr)}$",  "interpreter", "latex");
set(gca, "fontsize", fontsize, "TickLabelInterpreter", "latex");
set(gca, "fontsize", fontsize);
set(gca, "Tickdir", "out");
set(gcf, "Color", "white")
box on;
set(gca, "YLim", [-0.7, 0.7]);
set(gca, "XTick", [-200, -100, 0, 100, 200]);
set(gca, "YTick", [-0.50, -0.25, 0.00, 0.25, 0.50]);
set(gca, "YTickLabel", ["-0.50", "-0.25", "0.00", "0.25", "0.50"]);
legend([ve_handle, block_handle, csd_handle, perturbation_handle], "$v_\mathrm{t}$", "$v_\mathrm{b}$", "$v_\mathrm{sd}$", "$v_\mathrm{p}$", "location", "northwest", "interpreter", "latex", fontsize=18);
title("$v_\mathrm{t} = v_\mathrm{b} + v_\mathrm{sd} + v_\mathrm{p}$", "interpreter", "latex");
if SAVE_FIGURES
    export_fig("vt_vb_vsd_vp.pdf");
end

% Plot VE profiles and perturbations from steady state
figure("Position", [0, 0, 1400, 800]);
cspec = flipud(cool(numel(t_array)));
for i = 1:3
    subplot(2, 3, i);
    hold on;
    for j = 1:n_ve_profiles
        ve_handle = plot(x / 1e3, velocity_profile_store{i, j}, "-r", "LineWidth", 1.0, "Color", cspec(j, :));
    end
    steady_state_handle = plot(x / 1e3, v_steady_state, "-k", LineWidth=1.0);
    xlabel("$x \; \mathrm{(km)}$", "interpreter", "latex");
    ylabel("$v_\mathrm{t} \; \mathrm{(mm/yr)}$",  "interpreter", "latex");
    set(gca, "fontsize", fontsize, "TickLabelInterpreter", "latex");
    set(gca, "fontsize", fontsize);
    set(gca, "Tickdir", "out");
    set(gcf, "Color", "white")
    box on;
    set(gca, "YLim", [-10.0, 10.0]);
    set(gca, "XTick", [-200, -100, 0, 100, 200]);
    set(gca, "YTick", [-10.0, -5.0, 0.0, 5.0, 10.0]);
    set(gca, "YTickLabel", ["-10.0", "-5.0", "0.0", "5.0", "10.0"]);
    legend([ve_handle(1), steady_state_handle], "viscoelastic", "steady state", "location", "northwest")
    legend([ve_handle(1), steady_state_handle],"$v_\mathrm{t}$", "$v_\mathrm{sd}$", "Location", "northwest", "interpreter", "latex", "Fontsize", 18);
    title(strcat("$v_\mathrm{t} \; (n = ", string(n_array(i)), ")$"), "interpreter", "latex", "fontsize", fontsize, FontWeight="normal")

    subplot(2, 3, i + 3);
    hold on;
    for j = 1:n_ve_profiles
        ve_handle = plot(x / 1e3, sign(velocity_perturbation_profile_store{i, j}) .* (abs(velocity_perturbation_profile_store{i, j})).^(1.0/3.0), "-r", "LineWidth", 1.0, "Color", cspec(j, :));
    end
    steady_state_handle = plot(x / 1e3, zeros(size(v_steady_state)), "-k", LineWidth=1.0);
    set(gca, "fontsize", fontsize);
    set(gca, "Tickdir", "out");
    set(gcf, "Color", "white")
    box on;
    xlabel("$x \; \mathrm{(km)}$", "interpreter", "latex");
    ylabel("$v_\mathrm{p} \; \mathrm{(mm/yr) ^ {1/3}}$",  "interpreter", "latex");
    set(gca, "fontsize", fontsize, "TickLabelInterpreter", "latex");
    set(gca, "YLim", [-3.0, 3.0]);
    set(gca, "XTick", [-200, -100, 0, 100, 200]);
    set(gca, "YTick", [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]);
    set(gca, "YTickLabel", ["-3.0", "-2.0", "-1.0", "0.0", "1.0", "2.0", "3.0"]);
    legend([ve_handle(1), steady_state_handle],"$v_\mathrm{p}$", "$v_\mathrm{sd}$", "Location", "northwest", "interpreter", "latex", "Fontsize", 18);
    title(strcat("$v_\mathrm{p} \; (n = ", string(n_array(i)), ")$"), "interpreter", "latex", "fontsize", fontsize, FontWeight="normal")
end
if SAVE_FIGURES
    export_fig("velocity_profiles_and_perturbations.pdf");
end


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
G = zeros(numel(x), F.nfaults);
for i = 1:F.nfaults
   [~, uy, ~] = okada_1985(F.xf(i), F.yf, F.strike, F.depth, F.dip, F.length, ...
      F.fault_width, 1, 0, 0, x/1e3, zeros(size(x)), F.poissonsratio);
   G(:, i) = uy;
end

% Plot partials for each of the basal fault patches
figure(Position=[0, 0 1000, 600]);
set(gcf, "Color", "white");
hold on;
tiledlayout(5, 5, TileSpacing="none");
for i = 1:F.nfaults
    nexttile;
    hold on;

    % Draw patch width
    fh = fill([F.xf(i) - F.fault_width, F.xf(i), F.xf(i), F.xf(i) - F.fault_width], [0.0, 0.0, 0.5, 0.5], "y");
    fh.FaceColor = 0.70 * [1, 1, 1];
    fh.EdgeColor = "k";

    % Draw velocity profile
    plot(x / 1e3, G(:, i), "-b", LineWidth=2);
    text(-180, 0.43, string(i), "fontsize", fontsize, "interpreter", "latex");
    set(gca, "xlim", [-200, 200])
    set(gca, "ylim", [0, 0.5])
    box on;

    set(gca, "Tickdir", "out");
    set(gca, "XTick", []);
    set(gca, "YTick", []);
    if i == 16
        set(gca, "XTick", [-200, 0, 200]);
        set(gca, "YTick", [0.0, 0.25, 0.5]);
        set(gca, "YTickLabel", ["0.00", "0.25", "0.50"]);
        xlabel("$x \; \mathrm{(km)}$", "interpreter", "latex");
        ylabel("$v \; \mathrm{(mm/yr)}$",  "interpreter", "latex");
        set(gca, "fontsize", fontsize, "TickLabelInterpreter", "latex");
    end
end
if SAVE_FIGURES
    export_fig("basal_slip_partials.pdf");
end

% Solve for equivalent slip disribution
% for j = 1:n_ve_profiles
for j = 4
    markersize = 18;
    figure("Position", [0, 0, 1400, 800]);
    set(gcf, "Color", "w");

    for i = 1:3        
        v = velocity_perturbation_profile_store{i, j};
        v = v(:);
        if size(G, 1) > size(G, 2)
           effectiveslip = inv(G' * G) * G' * v; %#ok<*MINV> 
        elseif size(G, 2) > size(G, 1)
           effectiveslip = G' * inv(G * G') * v;
        end
        basal_model_velocities = G * effectiveslip;
        
        % Plot comparision of direct VE calculation and basal approximation
        subplot(2, 3, i)
        hold on;
        plot(x / 1e3, v, '--r', 'markersize', markersize, LineWidth=2.0);
        plot(x(1:20:end) / 1e3, basal_model_velocities(1:20:end), '.b', 'markersize', markersize);
        set(gca, "XLim", [-200, 200]);
        set(gca, "XTick", [-200, -100, 0, 100, 200]);
        set(gca, "YLim", [-2, 2]);
        set(gca, "YTick", [-2.0, -1.0, 0.0, 1.0, 2.0]);
        set(gca, "YTickLabel", ["-2.0", "-1.0", "0.0", "1.0", "2.0"]);
        xlabel("$x \; \mathrm{(km)}$", "interpreter", "latex");
        ylabel("$v \; \mathrm{(mm/yr)}$",  "interpreter", "latex");
        set(gca, "fontsize", fontsize);
        set(gca, "Tickdir", "out")
        box on;
        legend("$v_\mathrm{p}$", "$v_\mathrm{p^*}$", "Location", "northwest", "interpreter", "latex", "Fontsize", 18);
        title(strcat("$v_\mathrm{p^*} \approx v_\mathrm{p} \; (n = ", string(n_array(i)), ", \;", sprintf("t/T = %04.2f", t_array(j) / T / SECONDS_IN_A_YEAR), ")$"), "interpreter", "latex", "fontsize", fontsize, FontWeight="normal")
        set(gca, "TickLabelInterpreter", "latex");
        
        subplot(2, 3, i + 3)
        bar_x = F.xf - F.fault_width / 2;
        bh = bar(bar_x, effectiveslip, 1.0, "blue");
        bh.EdgeColor = [1, 1, 1];
        set(gca, "XLim", [-200, 200]);
        set(gca, "XTick", [-200, -100, 0, 100, 200]);
        set(gca, "YLim", [-4, 4]);
        set(gca, "YTick", [-4.0, -2.0, 0.0, 2.0, 4.0]);
        set(gca, "YTickLabel", ["-4.0", "-2.0", "0.0", "2.0", "4.0"]);
        xlabel("$x \; \mathrm{(km)}$", "interpreter", "latex");
        ylabel("$s_\mathrm{p^*} \; \mathrm{(mm/yr)}$",  "interpreter", "latex");
        set(gca, "fontsize", fontsize);
        set(gca, "Tickdir", "out")
        box on;
        title(strcat("$s_\mathrm{p^*} \; (n = ", string(n_array(i)), ", \; ", sprintf("t/T = %04.2f", t_array(j) / T / SECONDS_IN_A_YEAR), ")$"), "interpreter", "latex", "fontsize", fontsize, FontWeight="normal")
        set(gca, "TickLabelInterpreter", "latex");
    end
end
if SAVE_FIGURES
    export_fig("vp_sp.pdf");
end

