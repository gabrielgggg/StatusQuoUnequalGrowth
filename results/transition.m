clear; clc;
load results.mat;

initTau = 49.573 / 100 * 0.6;
initE = 17.094 / 100 * 0.6;
initG = initTau - initE;

dataYears = 1962:2019;

dataInPowerHouse = [
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  1
  1
  1
  1
  1
  1
  1
  1
  1
  1
  1
  0
  0
  0
  0
  1
  1
  1
  1
  1
  1
  1
  1
  0
  ];

dataInPowerSenate = [
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  1
  1
  1
  1
  1
  0
  0
  0
  0
  0
  0
  0
  0
  1
  1
  1
  1
  1
  1
  0
  0
  1
  1
  1
  1
  0
  0
  0
  0
  0
  0
  0
  0
  1
  1
  1
  1
  1
  ];

dataInPowerPresident = [
  1
  0
  0
  0
  0
  0
  0
  0
  0
  1
  1
  1
  1
  1
  1
  1
  1
  0
  0
  0
  0
  1
  1
  1
  1
  1
  1
  1
  1
  1
  1
  1
  1
  0
  0
  0
  0
  0
  0
  0
  0
  1
  1
  1
  1
  1
  1
  1
  1
  0
  0
  0
  0
  0
  0
  0
  0
  0
  1
  1
  ];

dataRatio = [
  % initE/(initE+initG) %
  % initE/(initE+initG) % imputed?
  0.279
  0.273
  0.283
  0.290
  0.280
  0.277
  0.294
  0.314
  0.337
  0.373
  0.403
  0.429
  0.443
  0.489
  0.491
  0.480
  0.483
  0.480
  0.487
  0.495
  0.507
  0.508
  0.488
  0.491
  0.487
  0.487
  0.491
  0.499
  0.532
  0.528
  0.548
  0.554
  0.570
  0.576
  0.596
  0.597
  0.609
  0.611
  0.608
  0.608
  0.601
  0.589
  0.580
  0.577
  0.581
  0.582
  0.584
  0.628
  0.587
  0.601
  0.614
  0.628
  0.640
  0.662
  0.672
  0.677
  0.667
  0.671
  ];

T = 70;
pathTau = ones([T, 2]) * initTau;
pathE = ones([T, 2]) * initE;
pathG = ones([T, 2]) * initG;

for tIx = 1:T-1
  pathTau(tIx+1, 1) = iP_tau(pathTau(tIx, 1), pathE(tIx, 1));
  pathE(tIx+1, 1) = iP_e(pathTau(tIx, 1), pathE(tIx, 1));
  pathG(tIx+1, 1) = pathTau(tIx+1, 1) - pathE(tIx+1, 1);

  pathTau(tIx+1, 2) = iR_tau(pathTau(tIx, 2), pathE(tIx, 2));
  pathE(tIx+1, 2) = iR_e(pathTau(tIx, 2), pathE(tIx, 2));
  pathG(tIx+1, 2) = pathTau(tIx+1, 2) - pathE(tIx+1, 2);
end
pathEG = pathE ./ (pathE + pathG);

sPathEG = pathEG;
% sPathEG(2:end, 1) = smooth(sPathEG(2:end, 1));
% sPathEG(2:end, 2) = smooth(sPathEG(2:end, 2));

figure('Name', 'In Power: Extremes');
% subplot(2, 2, 1);
% scatter( 1960, sPathEG(1, 1), 'ok'); hold on; %, 'filled'); hold on;
% scatter( 2025, 0.439, 'ok', 'filled');
scatter( dataYears, dataRatio, 'ok'); hold on; %, 'filled'); hold on;
plot( 1960:1960+T-1, sPathEG(:, 1), '-b', 'LineWidth', 2);
plot( 1960:1960+T-1, sPathEG(:, 2), '-r', 'LineWidth', 2 );
yline(sPathEG(1, 1), '--k', 'Optimum $\Delta_{low}$', 'Interpreter', 'latex', 'LineWidth', 1);
yline(0.439, '--k', 'Optimum $\Delta_{high}$', 'Interpreter', 'latex', 'LineWidth', 1);
xlim([1955 2025]);
% ylim([0 0.8]);
box on; xlabel('Time', 'Interpreter', 'latex');
ylabel('Entitlement Share $\frac{e}{e+g}$', 'Interpreter','latex');


timeVals = (dataYears(1)-1:dataYears(end))';
T = size(timeVals, 1);
pathTau = ones([T, 1]) * initTau;
pathE = ones([T, 1]) * initE;
pathG = ones([T, 1]) * initG;
inPower = zeros([T, 1]);

%
% House
%

pathTauLikeData = ones([T, 1]) * initTau;
pathElikeData = ones([T, 1]) * initE;
pathGlikeData = ones([T, 1]) * initG;
inPowerLikeData = zeros([T, 1]);

whichData = dataInPowerHouse;

for tIx = 2:T
  ppTau = iP_tau(pathTau(tIx-1), pathE(tIx-1));
  ppE = iP_e(pathTau(tIx-1), pathE(tIx-1));
  ppG = ppTau - ppE;
  ratioP = ppE / (ppE + ppG);

  rrTau = iR_tau(pathTau(tIx-1), pathE(tIx-1));
  rrE = iR_e(pathTau(tIx-1), pathE(tIx-1));
  rrG = rrTau - rrE;
  ratioR = rrE / (rrE + rrG);

  if tIx < T
    ratioData = dataRatio(tIx-1);

    if abs(ratioR - ratioData) < abs(ratioP -ratioData) % || tIx == 2
      inPower(tIx-1) = 1;
      pathTau(tIx) = rrTau;
      pathE(tIx) = rrE;
      pathG(tIx) = rrG;
    else
      inPower(tIx-1) = 0;
      pathTau(tIx) = ppTau;
      pathE(tIx) = ppE;
      pathG(tIx) = ppG;
    end

    if whichData(tIx-1) == 1 % Repubs
      inPowerLikeData(tIx-1) = 1;
      pathTauLikeData(tIx) = rrTau;
      pathElikeData(tIx) = rrE;
      pathGlikeData(tIx) = rrG;
    else % Dems
      inPowerLikeData(tIx-1) = 0;
      pathTauLikeData(tIx) = ppTau;
      pathElikeData(tIx) = ppE;
      pathGlikeData(tIx) = ppG;
    end
  else
    pathTau(tIx) = ppTau;
    pathE(tIx) =ppE;
    pathG(tIx) = ppG;

    pathTauLikeData(tIx) = ppTau;
    pathElikeData(tIx) =ppE;
    pathGlikeData(tIx) = ppG;
  end
end
plotRatio = pathE ./ (pathE + pathG);
plotRatioLikeData = pathElikeData ./ (pathElikeData + pathGlikeData);
sPlotRatio = plotRatio; % smooth?
sPlotRatioLikeData = plotRatioLikeData; % smooth?

figure('Name', 'In Power: Best Period-by-Period Fit');
% subplot(2, 2, 2);
scatter( dataYears, dataRatio, 'ok'); hold on; %, 'filled'); hold on;
for tIx = 1:T-1
  if inPower(tIx) == 1
    style = '-r';
  else
    style = '-b';
  end
  plot(timeVals(tIx:tIx+1), sPlotRatio(tIx:tIx+1), style, 'LineWidth', 2);
end
yline(sPathEG(1, 1), '--k', 'Optimum $\Delta_{low}$', 'Interpreter', 'latex', 'LineWidth', 1);
yline(0.439, '--k', 'Optimum $\Delta_{high}$', 'Interpreter', 'latex', 'LineWidth', 1);
xlim([1955 2025]);
% ylim([0 0.8]);
box on; xlabel('Time', 'Interpreter', 'latex');
ylabel('Entitlement Share $\frac{e}{e+g}$', 'Interpreter','latex');

fprintf('MSE Best-fit %f \n', sqrt( mean( (dataRatio - sPlotRatio(2:end)).^2 ) ) );

figure('Name', 'In Power: House Data');
% subplot(2, 2, 3);
scatter( dataYears, dataRatio, 'ok'); hold on; %, 'filled'); hold on;
for tIx = 1:T-1
  if inPowerLikeData(tIx) == 1
    style = '-r';
  else
    style = '-b';
  end
  plot(timeVals(tIx:tIx+1), sPlotRatioLikeData(tIx:tIx+1), style, 'LineWidth', 2);
end
yline(sPathEG(1, 1), '--k', 'Optimum $\Delta_{low}$', 'Interpreter', 'latex', 'LineWidth', 1);
yline(0.439, '--k', 'Optimum $\Delta_{high}$', 'Interpreter', 'latex', 'LineWidth', 1);
xlim([1955 2025]);
% ylim([0 0.8]);
box on; xlabel('Time', 'Interpreter', 'latex');
ylabel('Entitlement Share $\frac{e}{e+g}$', 'Interpreter','latex');

fprintf('MSE House %f \n', sqrt( mean( (dataRatio - sPlotRatioLikeData(2:end)).^2 ) ) );















%
% Senate
%

pathTauLikeData = ones([T, 1]) * initTau;
pathElikeData = ones([T, 1]) * initE;
pathGlikeData = ones([T, 1]) * initG;
inPowerLikeData = zeros([T, 1]);

whichData = dataInPowerSenate;

for tIx = 2:T
  ppTau = iP_tau(pathTau(tIx-1), pathE(tIx-1));
  ppE = iP_e(pathTau(tIx-1), pathE(tIx-1));
  ppG = ppTau - ppE;
  ratioP = ppE / (ppE + ppG);

  rrTau = iR_tau(pathTau(tIx-1), pathE(tIx-1));
  rrE = iR_e(pathTau(tIx-1), pathE(tIx-1));
  rrG = rrTau - rrE;
  ratioR = rrE / (rrE + rrG);

  if tIx < T
    if whichData(tIx-1) == 1 % Repubs
      inPowerLikeData(tIx-1) = 1;
      pathTauLikeData(tIx) = rrTau;
      pathElikeData(tIx) = rrE;
      pathGlikeData(tIx) = rrG;
    else % Dems
      inPowerLikeData(tIx-1) = 0;
      pathTauLikeData(tIx) = ppTau;
      pathElikeData(tIx) = ppE;
      pathGlikeData(tIx) = ppG;
    end
  else
    pathTauLikeData(tIx) = ppTau;
    pathElikeData(tIx) =ppE;
    pathGlikeData(tIx) = ppG;
  end
end
plotRatioLikeData = pathElikeData ./ (pathElikeData + pathGlikeData);
sPlotRatioLikeData = plotRatioLikeData; % smooth?

figure('Name', 'In Power: Senate Data');
% subplot(2, 2, 4);
scatter( dataYears, dataRatio, 'ok'); hold on; %, 'filled'); hold on;
for tIx = 1:T-1
  if inPowerLikeData(tIx) == 1
    style = '-r';
  else
    style = '-b';
  end
  plot(timeVals(tIx:tIx+1), sPlotRatioLikeData(tIx:tIx+1), style, 'LineWidth', 2);
end
yline(sPathEG(1, 1), '--k', 'Optimum $\Delta_{low}$', 'Interpreter', 'latex', 'LineWidth', 1);
yline(0.439, '--k', 'Optimum $\Delta_{high}$', 'Interpreter', 'latex', 'LineWidth', 1);
xlim([1955 2025]);
% ylim([0 0.8]);
box on; xlabel('Time', 'Interpreter', 'latex');
ylabel('Entitlement Share $\frac{e}{e+g}$', 'Interpreter','latex');

fprintf('MSE Senate %f \n', sqrt( mean( (dataRatio - sPlotRatioLikeData(2:end)).^2 ) ) );











%
% President
%

pathTauLikeData = ones([T, 1]) * initTau;
pathElikeData = ones([T, 1]) * initE;
pathGlikeData = ones([T, 1]) * initG;
inPowerLikeData = zeros([T, 1]);

whichData = dataInPowerPresident;

for tIx = 2:T
  ppTau = iP_tau(pathTau(tIx-1), pathE(tIx-1));
  ppE = iP_e(pathTau(tIx-1), pathE(tIx-1));
  ppG = ppTau - ppE;
  ratioP = ppE / (ppE + ppG);

  rrTau = iR_tau(pathTau(tIx-1), pathE(tIx-1));
  rrE = iR_e(pathTau(tIx-1), pathE(tIx-1));
  rrG = rrTau - rrE;
  ratioR = rrE / (rrE + rrG);

  if tIx < T
    if whichData(tIx-1) == 1 % Repubs
      inPowerLikeData(tIx-1) = 1;
      pathTauLikeData(tIx) = rrTau;
      pathElikeData(tIx) = rrE;
      pathGlikeData(tIx) = rrG;
    else % Dems
      inPowerLikeData(tIx-1) = 0;
      pathTauLikeData(tIx) = ppTau;
      pathElikeData(tIx) = ppE;
      pathGlikeData(tIx) = ppG;
    end
  else
    pathTauLikeData(tIx) = ppTau;
    pathElikeData(tIx) =ppE;
    pathGlikeData(tIx) = ppG;
  end
end
plotRatioLikeData = pathElikeData ./ (pathElikeData + pathGlikeData);
sPlotRatioLikeData = plotRatioLikeData; % smooth?

figure('Name', 'In Power: President Data');
% subplot(2, 2, 4);
scatter( dataYears, dataRatio, 'ok'); hold on; %, 'filled'); hold on;
for tIx = 1:T-1
  if inPowerLikeData(tIx) == 1
    style = '-r';
  else
    style = '-b';
  end
  plot(timeVals(tIx:tIx+1), sPlotRatioLikeData(tIx:tIx+1), style, 'LineWidth', 2);
end
yline(sPathEG(1, 1), '--k', 'Optimum $\Delta_{low}$', 'Interpreter', 'latex', 'LineWidth', 1);
yline(0.439, '--k', 'Optimum $\Delta_{high}$', 'Interpreter', 'latex', 'LineWidth', 1);
xlim([1955 2025]);
% ylim([0 0.8]);
box on; xlabel('Time', 'Interpreter', 'latex');
ylabel('Entitlement Share $\frac{e}{e+g}$', 'Interpreter','latex');

fprintf('MSE President %f \n', sqrt( mean( (dataRatio - sPlotRatioLikeData(2:end)).^2 ) ) );

























%
% House + Discretion
%

pathTauLikeData = ones([T, 1]) * initTau;
pathElikeData = ones([T, 1]) * initE;
pathGlikeData = ones([T, 1]) * initG;
inPowerLikeData = zeros([T, 1]);

whichData = dataInPowerHouse;

for tIx = 2:T
  ppTau = tauDiscretionP; % iP_tau(pathTau(tIx-1), pathE(tIx-1));
  ppE = eDiscretionP; % iP_e(pathTau(tIx-1), pathE(tIx-1));
  ppG = ppTau - ppE;
  ratioP = ppE / (ppE + ppG);

  rrTau = tauDiscretionR; % iR_tau(pathTau(tIx-1), pathE(tIx-1));
  rrE = eDiscretionR; % iR_e(pathTau(tIx-1), pathE(tIx-1));
  rrG = rrTau - rrE;
  ratioR = rrE / (rrE + rrG);

  if tIx < T
    ratioData = dataRatio(tIx-1);

    if whichData(tIx-1) == 1 % Repubs
      inPowerLikeData(tIx-1) = 1;
      pathTauLikeData(tIx) = rrTau;
      pathElikeData(tIx) = rrE;
      pathGlikeData(tIx) = rrG;
    else % Dems
      inPowerLikeData(tIx-1) = 0;
      pathTauLikeData(tIx) = ppTau;
      pathElikeData(tIx) = ppE;
      pathGlikeData(tIx) = ppG;
    end
  else
    pathTauLikeData(tIx) = ppTau;
    pathElikeData(tIx) =ppE;
    pathGlikeData(tIx) = ppG;
  end
end
plotRatioLikeData = pathElikeData ./ (pathElikeData + pathGlikeData);
sPlotRatioLikeData = plotRatioLikeData; % smooth?

figure('Name', 'Discretion: House Data');
% subplot(2, 2, 3);
scatter( dataYears, dataRatio, 'ok'); hold on; %, 'filled'); hold on;
for tIx = 1:T-1
  if inPowerLikeData(tIx) == 1
    style = '-r';
  else
    style = '-b';
  end
  plot(timeVals(tIx:tIx+1), sPlotRatioLikeData(tIx:tIx+1), style, 'LineWidth', 2);
end
yline(sPathEG(1, 1), '--k', 'Optimum $\Delta_{low}$', 'Interpreter', 'latex', 'LineWidth', 1);
yline(0.439, '--k', 'Optimum $\Delta_{high}$', 'Interpreter', 'latex', 'LineWidth', 1);
xlim([1955 2025]);
% ylim([0 0.8]);
box on; xlabel('Time', 'Interpreter', 'latex');
ylabel('Entitlement Share $\frac{e}{e+g}$', 'Interpreter','latex');

fprintf('MSE Discretion %f \n', sqrt( mean( (dataRatio - sPlotRatioLikeData(2:end)).^2 ) ) );
