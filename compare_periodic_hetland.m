clear

Cycles = 40;
Period = 10;
Print = 1;

delt = Period/100;
delx = 0.2;

x = [0:delx:10];


D = 1;

t  = [0:delt:Period*Cycles];
vt = t(1:end-1)+diff(t);

b = struct('type','periodic',...
	   'magnitude',1.0,...
	   'cycles',Cycles,...
	   'period',Period,...
	   'phase',0.0);

muC = 1.0;

Tau = 0.5;

disp(sprintf('\t\n\nMaxwell:'));
muM = 1.0;
etaM = Tau;
[phi,psi] = GIDVE_MaterialConstruct('Maxwell',...
				    'muM',muM,'etaM',etaM);
[UM N_mxl vt VM] = IseisDispVE(muC,phi,psi,40,x,t,b,D);

disp(sprintf('\t\n\nBurgers:'));
muM = 1.0;
muV = 1.0;
etaM = Tau*10;
etaV = Tau;
[phi,psi] = GIDVE_MaterialConstruct('Burgers',...
				    'muM',muM,'muV',muV,...
				    'etaM',etaM,'etaV',etaV);
[UB N_bb vt VB] = IseisDispVE(muC,phi,psi,25,x,t,b,D);

% plot displacements for first and last cycle at y=D and 2D
figure(1)
clf
tpos1 = find(t>=0&t<Period);
tpos2 = find(t>=Period*(Cycles-1)&t<Period*Cycles);
xpos1 = find(x==1);
xpos2 = find(x==2);
subplot(1,2,1)
plot(t(tpos1)./Period,UM(xpos1,tpos1),'b-',...
     t(tpos1)./Period,UM(xpos2,tpos1),'b--',...
     t(tpos1)./Period,UB(xpos1,tpos1),'r-',...
     t(tpos1)./Period,UB(xpos2,tpos1),'r--')
xlabel('time / Period','FontSize',14)
ylabel('displacement / slip','FontSize',14)
title('first cycle','FontSize',14)
legend('Maxwell (x=D)','Maxwell (x=2D)','Burgers (x=D)','Burgers (x=2D)')
subplot(1,2,2)
plot(t(tpos2)./Period,UM(xpos1,tpos2),'b-',...
     t(tpos2)./Period,UM(xpos2,tpos2),'b--',...
     t(tpos2)./Period,UB(xpos1,tpos2),'r-',...
     t(tpos2)./Period,UB(xpos2,tpos2),'r--')
xlabel('time / Period','FontSize',14)
ylabel('displacement','FontSize',14)
title('last cycle','FontSize',14)
legend('Maxwell (x=D)','Maxwell (x=2D)','Burgers (x=D)','Burgers (x=2D)')
     
% plot velocities for first and last cycle at y=D
figure(2)
clf
tpos1 = find(vt>0&vt<Period);
tpos2 = find(vt>Period*(Cycles-1)&vt<Period*Cycles);
plot(vt(tpos1)./Period,VM(xpos1,tpos1),'b--',...
     vt(tpos1)./Period,VB(xpos1,tpos1),'r--',...
     vt(tpos2)./Period-Cycles+1,VM(xpos1,tpos2),'b-',...
     vt(tpos2)./Period-Cycles+1,VB(xpos1,tpos2),'r-')
hold off
xlabel('time / Period','FontSize',14)
ylabel('velocity','FontSize',14)
title('last cycle','FontSize',14)
legend('Maxwell (first cycle)','Burgers (first cycle)',...
       'Maxwell (last cycle)','Bugers (last cycle)')

% plot the cycle invariance velocities
figure(3)
clf
% plot just the first two velocity profiles so that legend will key
% correctly
plot(x,VM(:,tpos1(1)),'b-',...
     x,VB(:,tpos1(1)),'r-');
hold on
plot(x,VM(:,tpos1(9:8:end)),'b-',...
     x,VB(:,tpos1(9:8:end)),'r-');
hold off
xlabel('distance / locking depth','FontSize',14)
ylabel('velocity','FontSize',14)
title('cycle invariant interseismic velocities','FontSize',14)
legend('Maxwell model','Burgers model')
