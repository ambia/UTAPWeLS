% Copyright 2026 Joaquin Ambia Garrido (UT Austin)
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
csf = uu.CSF;
well = csf.Wells.addElement;
well.Name = 'Well 1';

%%% Set up modeling interval and bed boundaries
ml = [1350 1380];
well.ModelingLimits = ml;
a_md = ml(1) + 0.2 + rand*3;
while a_md(end) <= ml(2) - 3
    a_md = [a_md; a_md(end) + 0.2 + rand*3]; %#ok<AGROW>
end
well.EM.Geometry.addBB('md', a_md);


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
% One value properties
props = well.EM.Properties;
props.setProperty( ...
    'Porosity, Total', a_phi, ...
    'Shale Concentration', a_CSh, ...
    'Porosity, Shale', a_PhiSh, ...
    'Water Saturation, Total', a_Sw);

% Salinity
props.setProperty( ...
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
props.Compositions.setComposition('Matrix', cat(3, a_qz, 1-a_qz), ...
    'Matrix Components', {'Quartz', 'Albite'});

% Shale Composition
a_kao = ones(numLayers, 1);
props.Compositions.setComposition('Shale Solid Composition', a_kao, ...
    'Shale Solid Components', {'Kaolinite'});

% Fluid Composition
a_meth = ones(numLayers, 1);
props.Compositions.setComposition('Fluid Composition', a_meth, ...
    'Fluid Components', {'Methane'});

%%% Calculate physical properties

% Temperature
well.EMCalculators.Temperature_Gradient.run;

% Pressure
well.EMCalculators.Pore_Pressure.RefPressure = 3600;
well.EMCalculators.Pore_Pressure.run;

% Water Resistivity
well.EMCalculators.Rw_Cw_T.run;

% Resistivity rock class model
defRC = csf.RockClassDB.getElement('Default');
resDefRC = defRC.Resistivity;
resDefRC.Model = 'Archie';
resDefRC.a = 1.2;
resDefRC.m = 2.1;
resDefRC.n = 2.5;

% Resistivity calculator
well.EMCalculators.Resistivity.run;

% Nuclear properties
nucCalc = well.EMCalculators.Nuclear;
nucCalc.InSituSalinity = true;
nucCalc.InSituPressure = true;
nucCalc.calculate;


% ==============================================================================
%%% Simulate

% Resistivity       ARC, 1-D sA, 825, si = 0.1524 m
resSim = well.Simulators.Resistivity;
resSim.Tool = 'ARC';
resSim.Variant = '1-D UT semi-Analytical';
resSim.Settings.ARC_Model = '825';
resSim.run;

% Nuclear           UT_Longhorn_WL
nucSim = well.Simulators.Nuclear;
nucSim.DoGR = true;
nucSim.DoDensity = true;
nucSim.DoPEF = true;
nucSim.DoNeutron = true;
nucSim.SamplingRate = 0.1524;
nucSim.run;



% ==============================================================================
%%% Add noise to simulations

% Resistivity       resSim.run;
rLS = well.getLogSet('ARC (1-D UT SA) Simulated');
for rL = rLS.a_logsNoDepth
    a_data = ut.units.convert(rL.RawData, ...
        ut.units.Resistivity.ohmm, ut.units.Resistivity.mS_m);
    a_data = a_data + a_data .* randn(size(a_data)) .* 0.03 + ...
        randn(size(a_data)) .* 0.5;
    rL.RawData = ut.units.convert(a_data, ...
        ut.units.Resistivity.mS_m, ut.units.Resistivity.ohmm);
end

% GR
rLS = well.getLogSet('Gamma Ray (Longhorn Wireline) Simulated Logs');
rL = rLS.getLog('ECGR');
a_data = rL.RawData;
a_data = a_data + a_data .* randn(size(a_data)) .* 0.02;
rL.RawData = a_data;

% Density
rLS = well.getLogSet('Density (Longhorn Wireline) Simulated Logs');
rL = rLS.getLog('\rho_\alpha');
a_data = rL.RawData;
a_data = a_data + randn(size(a_data)) .* 0.015;
rL.RawData = a_data;

% PEF
rLS = well.getLogSet('PEF (Longhorn Wireline) Simulated Logs');
rL = rLS.getLog('PEF');
a_data = rL.RawData;
a_data = a_data + a_data .* randn(size(a_data)) .* 0.005;
rL.RawData = a_data;

% NPHI
rLS = well.getLogSet('Neutron (Longhorn Wireline) Simulated Logs');
rL = rLS.getLog('NPHI');
a_data = rL.RawData;
a_data = a_data + a_data .* randn(size(a_data)) .* 0.02;
rL.RawData = a_data;



%===============================================================================
%%% Combine all logs into a single set

% Resistivity, use as log set to start getting all logs together
rLSLogs = well.getLogSet('ARC (1-D UT SA) Simulated');

% GR
rLS = well.getLogSet('Gamma Ray (Longhorn Wireline) Simulated Logs');
for ii = numel(rLS.a_logsNoDepth) : -1 : 1
    rL = rLS.a_logsNoDepth(ii);
    if ~strcmp(rL.Name, 'ECGR')
        rLS.removeLogs(rL)
    end
end
rLSLogs.updateLogs(rLS);

% Density
rLS = well.getLogSet('Density (Longhorn Wireline) Simulated Logs');
for ii = numel(rLS.a_logsNoDepth) : -1 : 1
    rL = rLS.a_logsNoDepth(ii);
    if ~strcmp(rL.Name, '\rho_\alpha')
        rLS.removeLogs(rL)
    else
        rL.Name = 'RHO';
    end
end
rLSLogs.updateLogs(rLS);

% PEF
rLS = well.getLogSet('PEF (Longhorn Wireline) Simulated Logs');
for ii = numel(rLS.a_logsNoDepth) : -1 : 1
    rL = rLS.a_logsNoDepth(ii);
    if ~strcmp(rL.Name, 'PEF')
        rLS.removeLogs(rL)
    end
end
rLSLogs.updateLogs(rLS);

% NPHI
rLS = well.getLogSet('Neutron (Longhorn Wireline) Simulated Logs');
for ii = numel(rLS.a_logsNoDepth) : -1 : 1
    rL = rLS.a_logsNoDepth(ii);
    if ~strcmp(rL.Name, 'NPHI')
        rLS.removeLogs(rL)
    end
end
rLSLogs.updateLogs(rLS);


