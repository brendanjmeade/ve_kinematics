close all;
clear all;

siay = 1;
n_pts = 100;
x = linspace(-200e3, 200e3, n_pts);
y = zeros(size(x));
t = 50 * siay;
H = 10e3;
mu = 30e10;
eta = 1e19;
T = 100 * siay;
v = savage_2000(x, y, t, H, mu, eta, T);
figure;
hold on
plot(x, v, '-r');
v = v - 1/pi .* atan(x./H);
plot(x, v, '-b');
v = v(:);

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

figure;
hold on;
for i = 1:F.nfaults
   [~, uy, ~] = okada_1985(F.xf(i), F.yf, F.strike, F.depth, F.dip, F.length, ...
      F.fault_width, 1, 0, 0, x/1e3, y, F.poissonsratio);
   G(:, i) = uy;
   plot(x / 1e3, uy);
end

% Solve for equivalent slip disribution
effectiveslip = G \ v;

fontsize = 18;
markersize = 15;
figure;
set(gcf, "Color", "w");
hold on;
plot(x / 1e3, v, '+k', 'markersize', markersize);
plot(x / 1e3, G * effectiveslip, 'or', 'markersize', markersize);
xlabel("x (km)");
ylabel("v (mm/yr)");
set(gca, "fontsize", fontsize);
set(gca, "Tickdir", "out")
box on;
legend("Savage and Prescott", "effective basal dislocations")



