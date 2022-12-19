clear;

params = readmatrix('parameters.tab', 'FileType', 'text', 'Delimiter', '\t');
ix = 1;
nn = params(ix); ix = ix + 1;
RULE = params(ix); ix = ix + 1;
q = params(ix); ix = ix + 1;
beta = params(ix); ix = ix + 1;
theta = params(ix); ix = ix + 1;
Y = params(ix); ix = ix + 1;
DeltaY = params(ix); ix = ix + 1;
yR = params(ix); ix = ix + 1;
yP = params(ix); ix = ix + 1;
xBar = params(ix); ix = ix + 1;
errTol = params(ix); ix = ix + 1;
dchoice = params(ix); ix = ix + 1;
cDiscretion = params(ix); ix = ix + 1;
gDiscretion = params(ix); ix = ix + 1;
tauDiscretionR = params(ix); ix = ix + 1;
eDiscretionR = params(ix); ix = ix + 1;
tauDiscretionP = params(ix); ix = ix + 1;
eDiscretionP = params(ix); ix = ix + 1;
evalR = params(ix); ix = ix + 1;
evalP = params(ix); % ix = ix + 1;
clear params ix;

RULE_TAU = 1;
RULE_E = 2;
RULE_BOTH = 3;
RULE_G = 4;
RULE_ALL = 5;
RULE_NONE = 6;
ruleLabel = '';
switch RULE
  case RULE_TAU
    ruleLabel = '\tau Rule';
  case RULE_E
    ruleLabel = 'e Rule';
  case RULE_BOTH
    ruleLabel = '\tau and e Rule';
  case RULE_G
    ruleLabel = 'g Rule';
  case RULE_ALL
    ruleLabel = '\tau, e, and g Rule';
  case RULE_NONE
    ruleLabel = 'Discretion';
end

%
% Too big! Don't save & load these unless needed.
%
% polR = loadBinary('polR.bin', 'float64', [ nn, nn ]);
% polP = loadBinary('polP.bin', 'float64', [ nn, nn ]);

tau = loadBinary('tau.bin', 'float64', [ nn, 1 ]);
e = loadBinary('e.bin', 'float64', [ nn, 1 ]);
g = loadBinary('g.bin', 'float64', [ nn, 1 ]);
cR = loadBinary('cR.bin', 'float64', [ nn, 1 ]);
cP = loadBinary('cP.bin', 'float64', [ nn, 1 ]);

Wr = loadBinary('Wr.bin', 'float64', [ nn, 1 ]);
Vr = loadBinary('Vr.bin', 'float64', [ nn, 1 ]);
Wp = loadBinary('Wp.bin', 'float64', [ nn, 1 ]);
Vp = loadBinary('Vp.bin', 'float64', [ nn, 1 ]);

polR_cR = loadBinary('polRcR.bin', 'float64', [ nn, 1 ]);
polR_cP = loadBinary('polRcP.bin', 'float64', [ nn, 1 ]);
polR_tau = loadBinary('polRtau.bin', 'float64', [ nn, 1 ]);
polR_e = loadBinary('polRe.bin', 'float64', [ nn, 1 ]);
polR_g = loadBinary('polRg.bin', 'float64', [ nn, 1 ]);

polP_cR = loadBinary('polPcR.bin', 'float64', [ nn, 1 ]);
polP_cP = loadBinary('polPcP.bin', 'float64', [ nn, 1 ]);
polP_tau = loadBinary('polPtau.bin', 'float64', [ nn, 1 ]);
polP_e = loadBinary('polPe.bin', 'float64', [ nn, 1 ]);
polP_g = loadBinary('polPg.bin', 'float64', [ nn, 1 ]);

stationary = loadBinary('stationary.bin', 'float64', [ nn, 1 ]);

iWr = scatteredInterpolant(tau, e, Wr, 'linear', 'none');
iVr = scatteredInterpolant(tau, e, Vr, 'linear', 'none');
iWp = scatteredInterpolant(tau, e, Wp, 'linear', 'none');
iVp = scatteredInterpolant(tau, e, Vp, 'linear', 'none');

iR_cR = scatteredInterpolant(tau, e, polR_cR, 'linear', 'none');
iR_cP = scatteredInterpolant(tau, e, polR_cP, 'linear', 'none');
iR_tau = scatteredInterpolant(tau, e, polR_tau, 'linear', 'none');
iR_e = scatteredInterpolant(tau, e, polR_e, 'linear', 'none');
iR_g = scatteredInterpolant(tau, e, polR_g, 'linear', 'none');

iP_cR = scatteredInterpolant(tau, e, polP_cR, 'linear', 'none');
iP_cP = scatteredInterpolant(tau, e, polP_cP, 'linear', 'none');
iP_tau = scatteredInterpolant(tau, e, polP_tau, 'linear', 'none');
iP_e = scatteredInterpolant(tau, e, polP_e, 'linear', 'none');
iP_g = scatteredInterpolant(tau, e, polP_g, 'linear', 'none');

save results;
