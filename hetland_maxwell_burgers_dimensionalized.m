clear all;
close all;

siay = 60 * 60 * 24 * 365.25;
cycles = 40; % number of cycles
period = 10 * siay; % earthquake every 10 years
t = [0:period / 10:period*cycles]; % s
x = linspace(-100e3, 100e3, 100);  % m
D = 15e3; % m
muC = 3.0e10; % Pa, reference shear modulus???
Maxwell.muM = 3.0e10; % Pa
Maxwell.etaM = 1e20;  % Pa s ???
Burgers.muM = 3.0e10;  % Pa
Burgers.muV = 3.0e10;  % Pa
Burgers.etaM = 1e19; % Pa s ???
Burgers.etaV = 1e18; % Pa s ???
b = struct('type', 'periodic', 'magnitude', 1.0, 'cycles', cycles, 'period', period, 'phase', 0.0);

[phi, psi] = GIDVE_MaterialConstruct('Maxwell', 'muM', Maxwell.muM, 'etaM', Maxwell.etaM);
[~, ~, ~, Maxwell.v] = IseisDispVE(muC, phi, psi, 40, x, t, b, D);
% [phi, psi] = GIDVE_MaterialConstruct('Burgers', 'muM', Burgers.muM, 'muV', Burgers.muV, 'etaM', Burgers.etaM, 'etaV', Burgers.etaV);
% [~, ~, ~, Burgers.v] = IseisDispVE(muC, phi, psi, 25, x, t, b, D);

% Plot the cycle invariance velocities
figure("color", "w");
hold on;
tpos = (cycles*period-period+1):1:cycles*period-1;
for i = 1:numel(tpos)
    plot(x, Maxwell.v(:, tpos(i)), 'b-');
    % plot(x, Burgers.v(:, tpos(i)), 'r-');
end
xlabel('distance / locking depth', 'FontSize', 14);
ylabel('velocity', 'FontSize', 14);
legend('Maxwell model', 'Burgers model');
box on;
set(gca, "TickDir", "out");
set(gca, "FontSize", 14);
