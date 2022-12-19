clear; clc;
load results.mat;

scaleDots = 5000.0;

figure('Position', [100 100 800 600]);
subplot(1, 2, 1);
scatter(cP, cR, stationary * scaleDots + 0.0001, stationary, 'd', 'filled'); hold on;
scatter( [ cDiscretion, xBar ], [ xBar, cDiscretion ] , 'ro', 'filled' );
plot( [ xBar, Y-xBar-xBar, xBar, xBar ], [ xBar, xBar, Y-xBar-xBar, xBar ], '--k', 'LineWidth', 1);
xlabel('Consumption of Poor');
ylabel('Consumption of Rich');
legend({'Stationary Distribution', 'Discretion'}, 'Location', 'Best');
legend('boxoff');

% figure;
subplot(1, 2, 2);
scatter(tau, e, stationary * scaleDots + 0.0001, stationary, 'd', 'filled'); hold on;
scatter( [ tauDiscretionR, tauDiscretionP ], [ eDiscretionR, eDiscretionP ] , 'ro', 'filled' );
plot( [ xBar+xBar-yP, yR-xBar, yR-xBar, xBar+xBar-yP ], [ xBar-yP, yR-xBar-xBar, xBar-yP, xBar-yP ], '--k', 'LineWidth', 1);
xlabel('Tax on Rich (\tau)');
ylabel('Entitlements to Poor (e)');
% legend({'Stationary Distribution', 'Discretion'}, 'Location', 'Best');
% legend('boxoff');

figure;
subplot(1, 2, 1);
scatter3(tau, e, Vr)
hold on;
scatter3(tau, e, Wr)
title('Rich: Value Functions');
legend('In power', 'Out of power', 'Location', 'best');
xlabel('$\bar \tau$', 'Interpreter', 'latex');
ylabel('$\bar e$', 'Interpreter', 'latex');
subplot(1, 2, 2);
scatter3(tau, e, Vp)
hold on;
scatter3(tau, e, Wp)
title('Poor: Value Functions');
legend('In power', 'Out of power', 'Location', 'best');
xlabel('$\bar \tau$', 'Interpreter', 'latex');
ylabel('$\bar e$', 'Interpreter', 'latex');

fixElow = 0.05;
fixEhigh = 0.15;
legLow = sprintf('$ \\bar e = %4.2f $', fixElow);
legHigh = sprintf('$ \\bar e = %4.2f $', fixEhigh);

figure;
subplot(1, 3, 1);
plot(tau, iP_cR(tau, fixElow * ones([nn, 1]))); hold on;
plot(tau, iP_cR(tau, fixEhigh * ones([nn, 1]))); 
xlabel('$\bar \tau$', 'Interpreter', 'latex');
legend(legLow, legHigh, 'Location', 'best', 'Interpreter', 'latex');
title('Poor in power: C_R');
ylim([ xBar, Y-xBar-xBar ]);

subplot(1, 3, 2);
plot(tau, iP_cP(tau, fixElow * ones([nn, 1]))); hold on;
plot(tau, iP_cP(tau, fixEhigh * ones([nn, 1]))); 
xlabel('$\bar \tau$', 'Interpreter', 'latex');
title('Poor in power: C_P');
ylim([ xBar, Y-xBar-xBar ]);

subplot(1, 3, 3);
plot(tau, iP_g(tau, fixElow * ones([nn, 1]))); hold on;
plot(tau, iP_g(tau, fixEhigh * ones([nn, 1]))); 
xlabel('$\bar \tau$', 'Interpreter', 'latex');
title('Poor in power: g');
ylim([ xBar, Y-2*xBar ]);
