% Copyright 2022 Joaquin Ambia Garrido (UT Austin)
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy 
% of this software and associated documentation files (the "Software"), to deal 
% in the Software without restriction, including without limitation the rights 
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
% of the Software, and to permit persons to whom the Software is furnished to do 
% so, subject to the following conditions:
% The above copyright notice and this permission notice shall be included in all 
% copies or substantial portions of the Software.
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
% INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
% PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
% HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION 
% OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
% SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

%%% This script shows multiple aspect of building a synthetic model:
% - Adding a well
% - Adding boundaries
% - Modifying earth model properties (including compositions)
% - Calculate physical properties
% - Running simulations

%%% Add well
csf = uu.rCSF;
well = csf.addNewWell;
well.name = 'Well 1';

%%% Set up modeling interval and bed boundaries
ml = [1350 1380];
well.ModelingLimits = ml;
a_md = ml(1) + 0.2 + rand*3;
while a_md(end) <= ml(2) - 3
    a_md = [a_md; a_md(end) + 0.2 + rand*3]; %#ok<AGROW>
end
well.rEM.addBB('md', a_md);


%%% Set up random composition

numLayers = numel(a_md) + 1;
% Porosity
a_phi = abs(randn(numLayers, 1)*0.04 + 0.22);
% CSh
a_CSh = abs(randn(numLayers, 1)*0.05 + 0.15);
% Shale Porosity
a_PhiSh = abs(randn(numLayers, 1)*0.01 + 0.05);
% Sw
a_Sw = abs(randn(numLayers, 1)*0.03 + 0.2);
% a_Sw = nan(numLayers, 1);
% for ii = 1:numLayers
%     a_Sw(ii) = randn*0.03 + 0.7 .* (ii/numLayers).^6 + 0.2;
% end
% One value properties
well.rEM.setProperty( ...
    'Porosity, Total', a_phi, ...
    'Shale Concentration', a_CSh, ...
    'Porosity, Shale', a_PhiSh, ...
    'Water Saturation, Total', a_Sw);

% Salinity
well.rEM.setProperty( ...
    'Salinity', 30000, ...
    'Shale Salinity', 30000);

% Matrix Composition
a_qz = nan(numLayers, 1);
a_qz(1) = 0.8 + rand*0.2;
for ii = 2:numLayers
    a_qz(ii) = a_qz(ii-1) + randn*0.04;
    while a_qz(ii) > 1
        a_qz(ii) = a_qz(ii) - 0.04;
    end
end
well.rEM.Compositions.setComposition('Matrix', cat(3, a_qz, 1-a_qz), ...
    'Matrix Components', {'Quartz', 'Albite'});

% Shale Composition
a_kao = ones(numLayers, 1);
well.rEM.Compositions.setComposition('Shale Solid Composition', a_kao, ...
    'Shale Solid Components', {'Kaolinite'});

% Fluid Composition
a_meth = ones(numLayers, 1);
well.rEM.Compositions.setComposition('Fluid Composition', a_meth, ...
    'Fluid Components', {'CH4'});

%%% Calculate physical properties

% Temperature
well.EM_Calculators.Temperature_GeothermalGradient.calculate;

% Pressure
well.EM_Calculators.PorePressure.RefPress = 3600;
well.EM_Calculators.PorePressure.calculatePorePress;

% Water Resistivity
well.EM_Calculators.WaterResistivity.calculate;

% Resistivity
well.EM_Calculators.Resistivities.calculate;
defRC = csf.rRockClassDatabase.s_rockClasses.Default;
defRC.ar_model = Calculators.E_ResistivityModel.Archie;
defRC.ar_a = 1.2;
defRC.ar_m = 2.1;
defRC.ar_n = 2.5;
well.EM_Calculators.Resistivities.calculate;

% Nuclear properties
well.EM_Calculators.Nuclear.InSituSalinity = true;
well.EM_Calculators.Nuclear.InSituPressure = true;
well.EM_Calculators.Nuclear.calculate;

% Sonic
csf.rComponentDatabase.s_SolidDB.Albite.slowness_p = 40;
csf.rComponentDatabase.s_SolidDB.Albite.slowness_s = 90;
csf.rComponentDatabase.s_SolidDB.Kaolinite.slowness_p = 230;
csf.rComponentDatabase.s_SolidDB.Kaolinite.slowness_s = 350;
defRC.rEffectiveMediumTheory.UseFromComponentDatabase = 'Slownesses';
well.EM_Calculators.EffectiveMediumTheory.calculate;


%%% Simulate

% Resistivity       ARC, 1-D sA, 825, si = 0.1524 m
resSim = well.Simulators.Resistivity;
resSim.run;

% Nuclear           UT_Longhorn_WL
nucSim = well.Simulators.Nuclear;
nucSim.DoGR = true;
nucSim.DoDensity = true;
nucSim.DoPEF = true;
nucSim.DoNeutron = true;
nucSim.SamplingRate = 0.1524;
nucSim.runSim;

% Sonic         Sonic Scanner
sonicSim = well.Simulators.SonicLog;
sonicSim.run;


%%% Add noise to simulations

% Resistivity       resSim.run;
rLS = well.getLogSet('ARC (1-D UT SA) Simulated');
for rL = rLS.a_logsNoDepth
    a_data = Units.convert(rL.rawData, ...
        Units.Resistivity.ohmm, Units.Resistivity.mS_m);
    a_data = a_data + a_data .* randn(size(a_data)) .* 0.03 + ...
        randn(size(a_data)) .* 0.5;
    rL.rawData = Units.convert(a_data, ...
        Units.Resistivity.mS_m, Units.Resistivity.ohmm);
end

% GR
rLS = well.getLogSet('Gamma Ray (Longhorn Wireline) Simulated Logs');
rL = rLS.getLog('ECGR');
a_data = rL.rawData;
a_data = a_data + a_data .* randn(size(a_data)) .* 0.02;
rL.rawData = a_data;

% Density
rLS = well.getLogSet('Density (Longhorn Wireline) Simulated Logs');
rL = rLS.getLog('\rho_\alpha');
a_data = rL.rawData;
a_data = a_data + randn(size(a_data)) .* 0.015;
rL.rawData = a_data;

% PEF
rLS = well.getLogSet('PEF (Longhorn Wireline) Simulated Logs');
rL = rLS.getLog('PEF');
a_data = rL.rawData;
a_data = a_data + a_data .* randn(size(a_data)) .* 0.005;
rL.rawData = a_data;

% NPHI
rLS = well.getLogSet('Neutron (Longhorn Wireline) Simulated Logs');
rL = rLS.getLog('NPHI');
a_data = rL.rawData;
a_data = a_data + a_data .* randn(size(a_data)) .* 0.02;
rL.rawData = a_data;

% Sonic P       sonicSim.run;
rLS = well.getLogSet('Sonic (Longhorn)');
rL = rLS.getLog('DTP');
a_data = rL.rawData;
a_data = a_data + a_data .* randn(size(a_data)) .* 0.01;
rL.rawData = a_data;
% S
rL = rLS.getLog('DTS');
a_data = rL.rawData;
a_data = a_data + a_data .* randn(size(a_data)) .* 0.01;
rL.rawData = a_data;


rLS = well.getLogSet('CAL');
rL = rLS.getLog('CAL');
a_data = rL.rawData;
a_data = a_data + randn(size(a_data)) .* 0.06;
rL.rawData = a_data;

%%% Compound of logs
% Resistivity, use as log set to start getting all logs together
rLSLogs = well.getLogSet('ARC (1-D UT SA) Simulated');

% GR
rLS = well.getLogSet('Gamma Ray (Longhorn Wireline) Simulated Logs');
for ii = numel(rLS.a_logsNoDepth) : -1 : 1
    rL = rLS.a_logsNoDepth(ii);
    if ~strcmp(rL.name, 'ECGR')
        rLS.removeLogs(rL)
    end
end
rLSLogs.updateFromLogSet(rLS);

% Density
rLS = well.getLogSet('Density (Longhorn Wireline) Simulated Logs');
for ii = numel(rLS.a_logsNoDepth) : -1 : 1
    rL = rLS.a_logsNoDepth(ii);
    if ~strcmp(rL.name, '\rho_\alpha')
        rLS.removeLogs(rL)
    else
        rL.name = 'RHO';
    end
end
rLSLogs.updateFromLogSet(rLS);

% PEF
rLS = well.getLogSet('PEF (Longhorn Wireline) Simulated Logs');
for ii = numel(rLS.a_logsNoDepth) : -1 : 1
    rL = rLS.a_logsNoDepth(ii);
    if ~strcmp(rL.name, 'PEF')
        rLS.removeLogs(rL)
    end
end
rLSLogs.updateFromLogSet(rLS);

% NPHI
rLS = well.getLogSet('Neutron (Longhorn Wireline) Simulated Logs');
for ii = numel(rLS.a_logsNoDepth) : -1 : 1
    rL = rLS.a_logsNoDepth(ii);
    if ~strcmp(rL.name, 'NPHI')
        rLS.removeLogs(rL)
    end
end
rLSLogs.updateFromLogSet(rLS);

% Sonic P       sonicSim.run;
rLS = well.getLogSet('Sonic (Longhorn)');
rLSLogs.updateFromLogSet(rLS);


%%% Mineralogy Logs
rLS = well.getLogSet('Key-Well Mineralogy');
for rL = rLS.a_logsNoDepth
    a_data = rL.rawData;
    a_data = a_data + a_data .* randn(size(a_data)) .* 0.001 + ...
        randn(size(a_data)) .* 0.0007;
    rL.rawData = a_data;
end

