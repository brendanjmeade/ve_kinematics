close all;
clear all;

siay = 1;
n_pts = 100;
x = linspace(-200e3, 200e3, n_pts);
y = zeros(size(x));
t = 5 * siay;
H = 15e3;
mu = 30e10;
eta = 1e19;
T = 100 * siay;

% Calculate partials for a single fault
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
G = zeros(n_pts, F.nfaults);

for i = 1:F.nfaults
   [~, G(:, i), ~] = okada_1985(F.xf(i), F.yf, F.strike, F.depth, F.dip, F.length, ...
      F.fault_width, 1, 0, 0, x/1e3, y, F.poissonsratio);
end

% Loop over times in earthquake cycle
t_vec = 1:5:50;
effective_slip_mat = zeros(F.nfaults, numel(t_vec));
for i = 1:numel(t_vec)
    % Forward velocities
    v = savage_2000(x, y, t_vec(i), H, mu, eta, T);
    v = v - 1 / pi .* atan(x ./ H);
    v = v(:);

    % Solve for equivalent slip disribution
    if size(G, 1) > size(G, 2)
        effective_slip = inv(G' * G) * G' * v;
    elseif size(G, 2) > size(G, 1)
        effective_slip = G' * inv(G * G') * v;
    end
    effective_slip_mat(:, i) = effective_slip;
end

fontsize = 18;
markersize = 15;

figure;
hold on
plot(x, v, '-r');
plot(x, v, '-b');

figure;
set(gcf, "Color", "w");
hold on;
plot(x / 1e3, v, '+k', 'markersize', markersize);
plot(x / 1e3, G * effective_slip, 'or', 'markersize', markersize);
xlabel("x (km)");
ylabel("v (mm/yr)");
set(gca, "fontsize", fontsize);
set(gca, "Tickdir", "out")
box on;
legend("Savage and Prescott", "effective basal dislocations")

figure;
set(gcf, "Color", "w");
hold on;
for i = 1:numel(t_vec)
    plot(effective_slip_mat(:, i), '-k+');
end
box on;
xlabel("dislocation patch index");
ylabel("v / s_0", "Fontsize", fontsize);
set(gca, "fontsize", fontsize);
set(gca, "Tickdir", "out")

