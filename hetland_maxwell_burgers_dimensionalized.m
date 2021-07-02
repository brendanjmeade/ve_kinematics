clear all;
close all;

cycles = 100;
period = 10;
delt = period / 10;
x = linspace(-100e3, 100e3, 100);
D = 15e3;
t  = [0:delt:period*cycles];
muC = 1.0;
tau = 0.5;
Maxwell.muM = 3.0e10; % Pa
Maxwell.etaM = tau;  % not Pa s
Burgers.muM = 3.0e10;  % Pa
Burgers.muV = 3.0e10;  % Pa
Burgers.etaM = tau * 10; % not Pa s
Burgers.etaV = tau; % not Pa s
b = struct('type', 'periodic', 'magnitude', 1.0, 'cycles', cycles, 'period', period, 'phase', 0.0);

[phi, psi] = GIDVE_MaterialConstruct('Maxwell', 'muM', Maxwell.muM, 'etaM', Maxwell.etaM);
[~, ~, ~, Maxwell.v] = IseisDispVE(muC, phi, psi, 40, x, t, b, D);
[phi, psi] = GIDVE_MaterialConstruct('Burgers', 'muM', Burgers.muM, 'muV', Burgers.muV, 'etaM', Burgers.etaM, 'etaV', Burgers.etaV);
[~, ~, ~, Burgers.v] = IseisDispVE(muC, phi, psi, 25, x, t, b, D);

% Plot the cycle invariance velocities
figure("color", "w");
hold on;
tpos = (cycles*period-period+1):1:cycles*period-1;
for i = 1:numel(tpos)
    plot(x, Maxwell.v(:, tpos(i)), 'b-');
    plot(x, Burgers.v(:, tpos(i)), 'r-');
end
xlabel('distance / locking depth', 'FontSize', 14);
ylabel('velocity', 'FontSize', 14);
legend('Maxwell model', 'Burgers model');
box on;
set(gca, "TickDir", "out");
set(gca, "FontSize", 14);
