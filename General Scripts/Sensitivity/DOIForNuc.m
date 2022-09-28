function DOIForNuc(utap, rWell, simType, logName)
% Copyright 2022 Bruce Klappauf (UT Austin)
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

% Compare DOI of Nuclear simulation vs invasion.
% Method:
%   Set up a well with multiple layers, each with a different invasion
%   radius. As the radius moves from the borehole to the deep formation, we
%   expect to see the response change from that of borehole/invasion fluids to
%   to that of formation fluids.
%
% utap: mandatory input of UTAPWeLS session object
% rWell: optional input of well to use. If empty this will make new well.
%   Otherwise it will delete all the layers in the well.
% simType: 'WL' or 'LWD'

% call with
%  DOIForNuc(uu, [], 'WL')
% =========================================================================
% Define values
% =========================================================================

% Now we first set hard coded values that we may want to change or to have
% as inputs later.  We put them here so they are easy to find.

zThick = 2; % Layer Thickness (m)
nRadThick = 21; % nubmer of zone size steps
maxRadThick = .5; % (m)
startMD = 1000;
% use default porosity
phiRef = .2;
SwF = .01;
SwI = 1;
Cw = 1000; %ppm

radThick = linspace(0, maxRadThick, nRadThick);
% TIP
% Define names like this and never type a constant more than once.
% 4 reasons: 
%   1) Easier to change a value since just do it in one place
%   2) Reduce errors from missed instances of the value somewhere
%   3) Easier to understand code with meaningful name than just a number.
%   4) Easier to find the place to change it if all likely values at top.

% Set resistivity tool and log to test
if nargin<4 || isempty(logName)    
    logName = 'NPHI';
end
   

if startsWith(upper(simType), 'WL')
    simTool = 'UT_Longhorn_WL';
    if logName(1)=='N'
        logSet = 'Neutron (Longhorn Wireline) Simulated Logs';
    else        
        logSet = 'Density (Longhorn Wireline) Simulated Logs';
    end
else
    simTool = 'UT_Longhorn_LWD';
    if logName(1)=='N'
        logSet = 'Neutron (Longhorn LWD) Simulated Logs';
    else        
        logSet = 'Density (Longhorn LWD) Simulated Logs';
    end
end
% =========================================================================
% Handle inputs and defaults
% =========================================================================

% assume we have a UTAPWeLS session started called 'utap'
rCSF = utap.rCSF;
% make new well if needed
if nargin<2 || isempty(rWell)   
    rW = rCSF.addNewWell('SensWell1');
else
    rW = rWell;  % Just to make variable easier to type
end

bhRad = rW.rEM.radB(1,1);

% =========================================================================
% Create the basic earth model
% =========================================================================
% We want nRadThick layers 
rW.ModelingLimits = [startMD, startMD + (nRadThick)*zThick];
rW.rEM.deleteBB('mdSegment', [-inf, inf]) % clear all current layers
dU = rW.rDisplayUnitSys.Distance; % get display distance unit for use below
% ref layer tops
bbTop = startMD:zThick:(startMD + (nRadThick)*zThick);
nBB = numel(bbTop);
% Now add all initial layers.
rW.rEM.addBB('md', bbTop);
% Note layer indices are indices for layer below BB index.
% Layer indices and value arrays should be vertical.
allLyrIdx = (2:nBB)';

% set porosity for all layers. Just need a single value if all the same.
rW.rEM.Compositions.setComposition('Fluid', {'CH4', 1}, ...
                                   'zoneIdx', 2:rW.rEM.numRadB)
rW.rEM.setProperty('propName', 'Porosity, Total', 'value', phiRef, ...
                   'layerIdx', allLyrIdx);
rW.rEM.setProperty('propName', 'Water Saturation, Total', 'value', SwF, ...
                   'layerIdx', allLyrIdx);               
rW.rEM.setProperty('propName', 'Salinity', 'value', Cw, ...
                   'layerIdx', allLyrIdx, 'zoneIdx', 'all');
         
% set different invasion radius for each layer
for lyrIdx = 3:(numel(radThick)+1)
    disp(lyrIdx);             % Show which iteration we are on
    rW.rEM.addRad('rad', bhRad+radThick(lyrIdx-1), ...
                'radUnits', 'm', 'layerIdx', lyrIdx)     
end

rW.rEM.setProperty('propName', 'Water Saturation, Total', 'value', SwI, ...
                   'layerIdx', allLyrIdx(2:end),...
                   'zoneIdx', 2);
% =========================================================================
% Get the simulator and set it up
% =========================================================================
               
% Get simulator and set up
rSim = rW.Simulators.Nuclear;
rSim.SimTool = simTool;
rSim.DoNeutron = true;
rSim.DoDensity = true;
rSim.RunNucCalcFirst = true;



% =========================================================================
% Run Sim and plot data 
% =========================================================================
rSim.runSim
a = getSimData(rW, logName, allLyrIdx, dU);
fmtnVal = a(1);
waterVal = a(end);

% make new plot figure
radThickCm = radThick*100;
unitline = ones(1, numel(radThick));

f1 = figure;
plot(radThickCm, a, 'b');
ax = f1.CurrentAxes;
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.YLabel.String = logName;
ax.XLabel.String = 'Invasion Radius (cm)';
ttl = 'Depth of Investigation';
ax.Title.String = ttl;
ax.NextPlot = 'add';


plot(radThickCm, unitline*fmtnVal, 'r--');
text(radThickCm(1)+2, fmtnVal+.005, 'Formation Value');

plot(radThickCm, unitline*waterVal, 'r--');
text(radThickCm(1)+2,waterVal+.005, 'Water (Invaded) Value');




% =========================================================================
% Helper Functions (used above)
% =========================================================================

function a = getSimData(rW, logName, allLyrIdx, dU)
    % This gets the simulated log we defined and extracts the values at the
    % center of the ref and val layers
    LS = rW.getLogSet(logSet);
    log = LS.getLog(logName);
    logDepths = Units.convert(LS.depthData, LS.depthDataUnits, dU);
    testDpth = rW.rEM.getMidLayerDepths(allLyrIdx);
    testVals = interp1(logDepths, log.rawData, testDpth);
    a = testVals(:);   
end

end

