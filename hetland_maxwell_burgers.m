clear all;
close all;

Cycles = 40;
Period = 10;
delt = Period / 100;
x = linspace(-10, 10, 100);

D = 1;
t  = [0:delt:Period*Cycles];
vt = t(1:end-1) + diff(t);
muC = 1.0;
Tau = 0.5;
b = struct('type', 'periodic', 'magnitude', 1.0, 'cycles', Cycles, 'period', Period, 'phase', 0.0);

disp(sprintf('\t\n\nMaxwell:'));
Maxwell.muM = 1.0;
Maxwell.etaM = Tau;
[phi, psi] = GIDVE_MaterialConstruct('Maxwell', 'muM', Maxwell.muM, 'etaM', Maxwell.etaM);
% [Maxwell.UM Maxwell.N_mxl Maxwell.vt Maxwell.VM] = IseisDispVE(muC, phi, psi, 40, x, t, b, D);
[~, ~, ~, Maxwell.v] = IseisDispVE(muC, phi, psi, 40, x, t, b, D);

disp(sprintf('\t\n\nBurgers:'));
Burgers.muM = 1.0;
Burgers.muV = 1.0;
Burgers.etaM = Tau * 10;
Burgers.etaV = Tau;
[phi, psi] = GIDVE_MaterialConstruct('Burgers', 'muM', Burgers.muM, 'muV', Burgers.muV, 'etaM', Burgers.etaM, 'etaV', Burgers.etaV);
[~, ~, ~, Burgers.v] = IseisDispVE(muC, phi, psi, 25, x, t, b, D);

% Plot the cycle invariance velocities
figure("color", "w");
tpos = find(vt>Period*(Cycles-1) & vt<Period*Cycles);
hold on;
plot(x, Maxwell.v(:, tpos(1)), 'b-', x, Burgers.v(:, tpos(1)), 'r-');
plot(x, Maxwell.v(:, tpos(9:8:end)), 'b-', x, Burgers.v(:, tpos(9:8:end)), 'r-');
xlabel('distance / locking depth', 'FontSize', 14)
ylabel('velocity', 'FontSize', 14)
title('cycle invariant interseismic velocities', 'FontSize', 14)
legend('Maxwell model', 'Burgers model')
box on;
set(gca, "TickDir", "out");
