function ThicknessSensitivityForRes(utap, rWell, resRef)
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

% Compare change in simulation value vs change in layer value for different
% layer thicknesses.
% utap: mandatory input of UTAPWeLS session object
% rWell: optional input of well to use. If empty this will make new well.
%   Otherwise it will delete all the layers in the well.
% resRef: Reference Water Resistivity. This is the water resistivity of the
%   background layers which provide a refernce for the change in
%   resistivity in the thin layers.  The amount of change is hard coded
%   below. 
%   Note this is an extra unnecessary step since we could set the
%   resistivities directly but we want to show the use of the calculator.

% Method:
% For the comparison we will alternate thick reference layers "Ref" in 
% which we will keep a constant base resistivity, and relatively thin 
% "Var" layers, in which we will vary both the thickness and the 
% resistivity.  The resistivity will have a logarithmic increase of 
% from resRef to nRng*resRef.

% =========================================================================
% Define values
% =========================================================================

% Now we first set hard coded values that we may want to change or to have
% as inputs later.  We put them here so they are easy to find.

nVar = 10; % Number of var layers. Different resistivity for each layer
nRng = 1000; % ie highest resistivity will be nRng*resRef
nThick = 6; % number of thickness variations to look at
dThick = 2; % Thicknesses will vary from dThick to nThick*dThick (m)
thickMin = 1; % Minimum thickness to test (m)
secThick = 15; % Section thickness (m) (ie size of 1 odd + 1 even layer)
phi = .25; % Porosity of all layers (this is default)
startMD = 1000;

% TIP
% Define names like this and never type a constant more than once.
% 4 reasons: 
%   1) Easier to change a value since just do it in one place
%   2) Reduce errors from missed instances of the value somewhere
%   3) Easier to understand code with meaningful name than just a number.
%   4) Easier to find the place to change it if all likely values at top.

% Set resistivity tool and log to test
toolName = 'ARC';
variantName = '1-D UT semi-Analytical';
logName = 'A40H';

% =========================================================================
% Handle inputs and defaults
% =========================================================================

% assume we have a UTAPWeLS session started called 'utap'
rCSF = utap.CSF;

% use default resRef if needed
if nargin<3 || isempty(resRef)   
    resRef = .025;
end

% make new well if needed
if nargin<2 || isempty(rWell) 
    wellName = sprintf('SensWellRes %.2f', resRef);
    rW = rCSF.Wells.addElement(wellName);
    wellName = rW.Name; % In case there was a duplicate name.
else
    rW = rWell;  % Just to make variable easier to type
    wellName = rW.Name;
end
% =========================================================================
% Create the basic earth model
% =========================================================================
% We want 2*nVar +1 layers where the odd layers will be reference layers 
% and the even layers will have variable resistivities for and the current
% test thickness. 

rW.ModelingLimits = [startMD, startMD + (nVar+1)*secThick];
rW.EM.Geometry.deleteBB('mdSegment', [-inf, inf]) % clear all current layers
dU = rW.rDispUS.Distance; % get display distance unit for use below
% ref layer tops
bbRef = startMD:secThick:(startMD + secThick*(nVar+1));
% make initial variable layers that will be moved each thickness loop
bbVar = bbRef(2:(nVar+1)) - thickMin;
allBB = sort([bbRef, bbVar]);
nBB = numel(allBB);
% Now add all initial layers.
rW.EM.Geometry.addBB('md', allBB);
% Note layer indices are indices for layer below BB index.
% Layer indices and value arrays should be vertical.
allLyrIdx = (2:nBB)'; % layers BETWEEN first and last BB
varLyrIdx = (3:2:nBB)';
% set porosity for all layers. Just need a single value if all the same.
rW.EM.Properties.setProperty('propName', 'Porosity, Total', 'value', phi, ...
                   'layerIdx', allLyrIdx);
% set initial water resistivity               
rW.EM.Properties.setProperty('propName', 'Water Resistivity', 'value', resRef, ...
                   'layerIdx', allLyrIdx);  
% set the variable resistivities in the even layers               
resVar = logspace(log10(resRef), log10(nRng*resRef), nVar);
rW.EM.Properties.setProperty('propName', 'Water Resistivity', 'value', resVar', ...
                   'layerIdx', varLyrIdx);
         
% =========================================================================
% Get the simulator and set it up
% =========================================================================
               
% Get simulator and set up
rSim = rW.Simulators.Resistivity;
rSim.Tool = toolName;
rSim.Variant = variantName;

% =========================================================================
% Calculate the data sets
% =========================================================================

% Want to plot log vs res for different thickenesses, 
% so get Earth Model Resistivities. 
% These will not change in thickness loop.
rW.EMCalculators.Archies.run;  
resVarii = rW.EM.Properties.getProperty(...
                'propName', 'Resistivity (Parallel or Homogeneous)', ...
             	'layerIdx', varLyrIdx);
% Data array containing Earth Model resistivity different in var layers
a_EMDiff = log10(resVarii) - log10(resVarii(1));

% Make additional storage for each thickness (and var layer resistivity)
a_thick = nan(1, nThick);
a_dRes = nan(nVar, nThick);
a_MaxRes = nan(nVar, nThick);

% In a loop, get the simulated resistivity curve for each thicknes, and
% record the center value at each layer to compare to the earth model
for d = 1:nThick
    disp(d);             % Show which iteration we are on
    dMD = dThick*(d-1);  % first step is 0 (use initial layer positions)
    newMD = bbVar - dMD;
    a_thick(d) = thickMin + dMD;  % Record each var layer thickness
    
    % Move the boundaries of the var layers to change thickness
    for ii = 1:numel(newMD)
        rW.EM.Geometry.moveBB('idx', ii*2, 'md', newMD(ii))
    end    
     
    % Run the sim
    rSim.run
    % get log data as array of ref layers and var layers in columns
    a = getSimData(rSim, rW, logName, nVar, varLyrIdx, dU);
    % get the diff of the log of resistivities
    a_dRes(:, d) = diff(log10(a), [], 2);
    % get the raw value of the log of resistivities
    a_MaxRes(:, d) = log10(a(:,2));
end

% make new plot figure
f1 = figure;
f1.NumberTitle = 'off';
plot([0,a_EMDiff(end)+.5], [0, a_EMDiff(end)+.5], 'r--');
labelii = sprintf('%.1f m Thick', a_thick(1));
text(a_EMDiff(end), a_dRes(end,1), labelii);
% text(a_EMDiff(end), a_MaxRes(end,1), labelii);
ax = f1.CurrentAxes;
ttl = sprintf('Well:%s:: Thickness Sensitivity (R_ref = %.3f)', ...
              wellName, resVarii(1));
ax.Title.String = ttl;
f1.Name = ttl;
f1.Tag = ttl;
ax.NextPlot = 'add';
ax.XLabel.String = 'log(R_{EMVar}) - log(R_{EMRef})';
ax.YLabel.String = 'log(R_{SimVar}) - log(R_{SimRef})';
% ax.YLabel.String = 'log(R_{SimVar})';

for ii = 1:nThick
    plot(a_EMDiff, a_dRes(:,ii));
%     plot(a_EMDiff, a_MaxRes(:,ii));
    labelii = sprintf('%.1f m Thick', a_thick(ii));
    text(a_EMDiff(end), a_dRes(end,ii), labelii);
%     text(a_EMDiff(end), a_MaxRes(end,ii), labelii);
end



% =========================================================================
% Helper Functions (used above)
% =========================================================================

function a = getSimData(rRes, rW, logname, nRes, lidxVar, dU)
    % This gets the simulated log we defined and extracts the values at the
    % center of the ref and val layers
    LS = rW.getLogSet(rRes.OutputLogSet);
    log = LS.getLog(logname);
    logDepths = ut.units.convert(LS.DepthData, LS.DepthDataUnits, dU);
    testDpth = rW.EM.Geometry.getMidLayerDepths(lidxVar(1:nRes));
    refDpth = rW.EM.Geometry.getMidLayerDepths(lidxVar(1:nRes)-1);
    testVals = interp1(logDepths, log.RawData, testDpth);
    refVals = interp1(logDepths, log.RawData, refDpth);
    a = [refVals(:), testVals(:)];   
end

end

