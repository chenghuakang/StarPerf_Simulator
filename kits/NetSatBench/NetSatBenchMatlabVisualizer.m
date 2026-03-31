function NetSatBenchMatlabVisualizer(constellationXmlFile, hdf5File, userXmlFile, gatewayXmlFile, varargin)
%NETSATBENCHMATLABVISUALIZER
%
% Visualize one selected shell from the generated HDF5 positions so that
% satellite ordering matches the Python StarPerf generation order.
%
% Inputs:
%   constellationXmlFile : XML describing the constellation
%   hdf5File             : generated HDF5 file with /position/shellN/timeslotM
%   userXmlFile          : XML describing users
%   gatewayXmlFile       : XML describing gateways
%
% Name-value parameters:
%   "SelectedShell"            : shell index to visualize (default: 1)
%   "AddUserAccess"            : create sat-user access objects from HDF5-selected pairs
%   "AddGatewayAccess"         : create sat-gateway access objects from HDF5-selected pairs
%   "AddISL"                   : create ISL access objects from HDF5-selected pairs
%   "InterPlaneOffset"         : slot offset for adjacent-plane links
%   "UserMinElevationAngle"    : user elevation mask [deg]
%   "GatewayMinElevationAngle" : gateway elevation mask [deg]
%   "StartTime"                : scenario start time
%   "StopTime"                 : scenario stop time; if empty, derived from HDF5
%   "SampleTime"               : fallback scenario sample time [s] if HDF5 /info dT is absent
%   "CacheFile"                : MAT file path for cached visualization inputs
%   "UseCache"                 : load and reuse cached visualization inputs
%   "ShowDetails"              : viewer detail flag

%% ------------------------------------------------------------------------
% Options
% -------------------------------------------------------------------------
p = inputParser;

addParameter(p, "SelectedShell", 1, @(x)isnumeric(x) && isscalar(x) && x >= 1);

addParameter(p, "AddUserAccess", false, @(x)islogical(x) || isnumeric(x));
addParameter(p, "AddGatewayAccess", false, @(x)islogical(x) || isnumeric(x));
addParameter(p, "AddISL", false, @(x)islogical(x) || isnumeric(x));
addParameter(p, "InterPlaneOffset", 0, @(x)isnumeric(x) && isscalar(x));
addParameter(p, "UserMinElevationAngle", 25, @(x)isnumeric(x) && isscalar(x));
addParameter(p, "GatewayMinElevationAngle", 25, @(x)isnumeric(x) && isscalar(x));
addParameter(p, "StartTime", datetime(2023,10,1,0,0,0), @(x)isdatetime(x));
addParameter(p, "StopTime", [], @(x)isempty(x) || isdatetime(x));
addParameter(p, "SampleTime", 30, @(x)isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, "CacheFile", "", @(x)ischar(x) || isstring(x));
addParameter(p, "UseCache", true, @(x)islogical(x) || isnumeric(x));
addParameter(p, "ShowDetails", false, @(x)islogical(x) || isnumeric(x));

parse(p, varargin{:});

selectedShell     = round(p.Results.SelectedShell);

addUserAccess     = logical(p.Results.AddUserAccess);
addGatewayAccess  = logical(p.Results.AddGatewayAccess);
addISL            = logical(p.Results.AddISL);

interPlaneOffset  = round(p.Results.InterPlaneOffset);

userMinEl         = p.Results.UserMinElevationAngle;
gatewayMinEl      = p.Results.GatewayMinElevationAngle;

startTime         = p.Results.StartTime;
stopTime          = p.Results.StopTime;
sampleTime        = p.Results.SampleTime;
cacheFile         = string(p.Results.CacheFile);
useCache          = logical(p.Results.UseCache);
showDetails       = logical(p.Results.ShowDetails);

usingDefaults = string(p.UsingDefaults);

if strlength(cacheFile) == 0
    cacheFile = defaultCacheFile(hdf5File, selectedShell);
end

cacheHit = useCache && isfile(cacheFile);
cacheMode = "rebuild";

%% ------------------------------------------------------------------------
% Read or load visualization inputs
% -------------------------------------------------------------------------
if cacheHit
    S = load(cacheFile);
    if hasCachedScenario(S, selectedShell, addUserAccess, addGatewayAccess, addISL, startTime, stopTime)
        cacheMode = "scenario";
        sc = S.sc;
        sat = S.sat;
        usrGS = S.usrGS;
        gwGS = S.gwGS;
        if isfield(S, 'islAccesses')
            islAccesses = S.islAccesses; %#ok<NASGU>
        end
        if isfield(S, 'userAccesses')
            userAccesses = S.userAccesses; %#ok<NASGU>
        end
        if isfield(S, 'gatewayAccesses')
            gatewayAccesses = S.gatewayAccesses; %#ok<NASGU>
        end
        constData = S.constData;
        shellData = S.shellData;
        altitude_km = S.altitude_km;
        inclination_deg = S.inclination_deg;
        phase_shift = S.phase_shift;
        numPlanes = S.numPlanes;
        satsPerPlane = S.satsPerPlane;
        totalSatellites = S.totalSatellites;
        sampleTime = S.sampleTime;
        timeOfArrival = S.timeOfArrival;
        latDeg = S.latDeg;
        lonDeg = S.lonDeg;
        altMeters = S.altMeters;
        usrData = S.usrData;
        gwData = S.gwData;
        userSatPairs = S.userSatPairs;
        gatewaySatPairs = S.gatewaySatPairs;
        islPairs = S.islPairs;
        userMinEl = S.userMinEl;
        gatewayMinEl = S.gatewayMinEl;
        if isfield(S, 'effectiveStopTime')
            stopTime = S.effectiveStopTime;
        elseif isempty(stopTime)
            stopTime = startTime + seconds(timeOfArrival(end));
        end
        fprintf('Loaded visualization scenario cache from %s.\n', cacheFile);
    else
        fprintf('Ignoring non-scenario visualization cache at %s and rebuilding it.\n', cacheFile);
    end
end

if cacheMode == "rebuild"
    constData = readConstellationXml(constellationXmlFile);

    if selectedShell > constData.NumberOfShells
        error('SelectedShell=%d, but the XML contains only %d shell(s).', ...
            selectedShell, constData.NumberOfShells);
    end

    shellData = constData.Shells(selectedShell);

    altitude_km      = shellData.Altitude_km;
    inclination_deg  = shellData.Inclination_deg;
    phase_shift      = shellData.PhaseShift;
    numPlanes        = shellData.NumPlanes;
    satsPerPlane     = shellData.SatsPerPlane;
    totalSatellites  = numPlanes * satsPerPlane;

    [sampleTimeFromHdf5, timeOfArrival, latDeg, lonDeg, altMeters] = readShellPositionsFromHdf5( ...
        hdf5File, selectedShell, sampleTime);
    sampleTime = sampleTimeFromHdf5;

    minimumElevationFromHdf5 = readMinimumElevationFromHdf5(hdf5File);
    if ~isnan(minimumElevationFromHdf5)
        if any(usingDefaults == "UserMinElevationAngle")
            userMinEl = minimumElevationFromHdf5;
        end
        if any(usingDefaults == "GatewayMinElevationAngle")
            gatewayMinEl = minimumElevationFromHdf5;
        end
    end

    usrData = readTerminalXml(userXmlFile);
    gwData  = readTerminalXml(gatewayXmlFile);
    userSatPairs = false(totalSatellites, numel(usrData));
    gatewaySatPairs = false(totalSatellites, numel(gwData));
    islPairs = [];

    if addISL
        islPairs = readISLPairsFromHdf5(hdf5File, selectedShell, totalSatellites);
    end
    if addUserAccess && ~isempty(usrData)
        userSatPairs = readTerminalSatPairsFromHdf5( ...
            hdf5File, selectedShell, totalSatellites, totalSatellites + numel(gwData), numel(usrData));
    end
    if addGatewayAccess && ~isempty(gwData)
        gatewaySatPairs = readTerminalSatPairsFromHdf5( ...
            hdf5File, selectedShell, totalSatellites, totalSatellites, numel(gwData));
    end

end

if isempty(stopTime)
    stopTime = startTime + seconds(timeOfArrival(end));
end

islAccesses = [];
userAccesses = [];
gatewayAccesses = [];

%% ------------------------------------------------------------------------
% Create scenario
% -------------------------------------------------------------------------
if cacheMode ~= "scenario"
    sc = satelliteScenario(startTime, stopTime, sampleTime);
end

%% ------------------------------------------------------------------------
% Create satellites from HDF5 trajectories
% -------------------------------------------------------------------------
if cacheMode ~= "scenario"
    sat = [];
    for satIdx = 1:totalSatellites
        traj = geoTrajectory( ...
            [latDeg(:, satIdx), lonDeg(:, satIdx), altMeters(:, satIdx)], ...
            timeOfArrival);
        sat = [sat; createTrajectoryObject(sc, traj, sprintf("Sat%d", satIdx))]; %#ok<AGROW>
        fprintf('Created satellite %d with trajectory from HDF5.\n', satIdx);
    end

    try
        for satIdx = 1:numel(sat)
            sat(satIdx).MarkerSize = 4;
        end
    catch
    end
end

%% ------------------------------------------------------------------------
% Add users as ground stations with constant elevation mask
% -------------------------------------------------------------------------
if cacheMode ~= "scenario"
    usrGS = [];
    if ~isempty(usrData)
        usrLat  = [usrData.Latitude]';
        usrLon  = [usrData.Longitude]';
        usrName = string({usrData.Id})';

        [azU, elU] = constantElevationMask(userMinEl);

        usrGS = groundStation(sc, usrLat, usrLon, ...
            Name=usrName, ...
            MaskAzimuthEdges=azU, ...
            MaskElevationAngle=elU);
        
        for k = 1:numel(usrGS)
            usrGS(k).MarkerColor = [0 0.447 0.741];
        end
        
        fprintf('Created %d users from XML.\n', numel(usrGS));
    end
end


%% ------------------------------------------------------------------------
% Add gateways as ground stations with constant elevation mask
% -------------------------------------------------------------------------
if cacheMode ~= "scenario"
    gwGS = [];
    if ~isempty(gwData)
        gwLat  = [gwData.Latitude]';
        gwLon  = [gwData.Longitude]';
        gwName = string({gwData.Id})';

        [azG, elG] = constantElevationMask(gatewayMinEl);

        gwGS = groundStation(sc, gwLat, gwLon, ...
            Name=gwName, ...
            MaskAzimuthEdges=azG, ...
            MaskElevationAngle=elG);

        for k = 1:numel(gwGS)
            gwGS(k).MarkerColor = [0.850 0.325 0.098];
        end
        
        fprintf('Created %d gateways from XML.\n', numel(gwGS));
    end
end


%% ------------------------------------------------------------------------
% Generate ISL topology from HDF5 delay matrices
% -------------------------------------------------------------------------
if addISL
    if cacheMode ~= "scenario"
        islAccesses = [];
        for k = 1:size(islPairs,1)
            islAccesses = [islAccesses; access(sat(islPairs(k,1)), sat(islPairs(k,2)))]; %#ok<AGROW>
        end
    end
    fprintf('Loaded %d ISL pairs from HDF5.\n', size(islPairs,1));
end

%% ------------------------------------------------------------------------
% Satellite-user access
% -------------------------------------------------------------------------
if addUserAccess && ~isempty(usrGS)
    if cacheMode ~= "scenario"
        userAccesses = [];
        for u = 1:numel(usrGS)
            satIndices = find(userSatPairs(:, u)).';
            for satIdx = satIndices
                userAccesses = [userAccesses; access(sat(satIdx), usrGS(u))]; %#ok<AGROW>
            end
            fprintf('User %s linked to %d satellites from HDF5.\n', string({usrData(u).Id}), numel(satIndices));
        end
    end
    fprintf('Loaded HDF5-driven user-satellite link sets for %d users.\n', numel(usrGS));
end

%% ------------------------------------------------------------------------
% Satellite-gateway access
% -------------------------------------------------------------------------
if addGatewayAccess && ~isempty(gwGS)
    if cacheMode ~= "scenario"
        gatewayAccesses = [];
        for g = 1:numel(gwGS)
            satIndices = find(gatewaySatPairs(:, g)).';
            for satIdx = satIndices
                gatewayAccesses = [gatewayAccesses; access(sat(satIdx), gwGS(g))]; %#ok<AGROW>
            end
            fprintf('Gateway %s linked to %d satellites from HDF5.\n', string({gwData(g).Id}), numel(satIndices));
        end
    end
    fprintf('Loaded HDF5-driven gateway-satellite link sets for %d gateways.\n', numel(gwGS));
end

if cacheMode ~= "scenario"
    effectiveStopTime = stopTime;
    cacheVersion = 2; %#ok<NASGU>
    cacheConfig = struct( ...
        'SelectedShell', selectedShell, ...
        'AddUserAccess', addUserAccess, ...
        'AddGatewayAccess', addGatewayAccess, ...
        'AddISL', addISL, ...
        'StartTime', startTime, ...
        'StopTime', stopTime);

    try
        save(cacheFile, ...
            'cacheVersion', 'cacheConfig', ...
            'constData', 'shellData', ...
            'altitude_km', 'inclination_deg', 'phase_shift', 'numPlanes', 'satsPerPlane', 'totalSatellites', ...
            'sampleTime', 'timeOfArrival', 'latDeg', 'lonDeg', 'altMeters', ...
            'usrData', 'gwData', 'userSatPairs', 'gatewaySatPairs', 'islPairs', ...
            'userMinEl', 'gatewayMinEl', 'effectiveStopTime', ...
            'sc', 'sat', 'usrGS', 'gwGS', ...
            'islAccesses', 'userAccesses', 'gatewayAccesses');
        fprintf('Saved full visualization scenario cache to %s.\n', cacheFile);
    catch cacheErr
        warning('Could not save full visualization scenario cache to %s: %s', cacheFile, cacheErr.message);
    end
end

%% ------------------------------------------------------------------------
% Viewer
% -------------------------------------------------------------------------
satelliteScenarioViewer(sc, ShowDetails=showDetails);
play(sc);

%% ------------------------------------------------------------------------
% Summary
% -------------------------------------------------------------------------
fprintf('\nSelected shell summary\n');
fprintf('----------------------\n');
fprintf('Selected shell             : %d / %d\n', selectedShell, constData.NumberOfShells);
fprintf('Altitude [km]              : %.2f\n', altitude_km);
fprintf('Inclination [deg]          : %.2f\n', inclination_deg);
fprintf('Planes                     : %d\n', numPlanes);
fprintf('Satellites per plane       : %d\n', satsPerPlane);
fprintf('Total satellites           : %d\n', totalSatellites);
fprintf('Users                      : %d\n', numel(usrData));
fprintf('Gateways                   : %d\n', numel(gwData));
fprintf('ISL pairs                  : %d\n', size(islPairs,1));
fprintf('User min elevation [deg]   : %.2f\n', userMinEl);
fprintf('Gateway min elevation [deg]: %.2f\n', gatewayMinEl);

end

%% =========================================================================
% Helper functions
% =========================================================================

function [maskAz, maskEl] = constantElevationMask(minElevationDeg)
maskAz = [0 360];
maskEl = minElevationDeg;
end

function cacheFile = defaultCacheFile(hdf5File, selectedShell)
[folder, baseName, ~] = fileparts(char(hdf5File));
cacheFile = string(fullfile(folder, sprintf('%s_shell%d_vizcache.mat', baseName, selectedShell)));
end

function tf = hasCachedScenario(S, selectedShell, addUserAccess, addGatewayAccess, addISL, startTime, stopTime)
requiredFields = { ...
    'cacheVersion', 'cacheConfig', ...
    'constData', 'shellData', ...
    'altitude_km', 'inclination_deg', 'phase_shift', 'numPlanes', 'satsPerPlane', 'totalSatellites', ...
    'sampleTime', 'timeOfArrival', 'latDeg', 'lonDeg', 'altMeters', ...
    'usrData', 'gwData', 'userSatPairs', 'gatewaySatPairs', 'islPairs', ...
    'userMinEl', 'gatewayMinEl', ...
    'sc', 'sat', 'usrGS', 'gwGS'};
tf = all(isfield(S, requiredFields));
if ~tf
    return;
end

cfg = S.cacheConfig;
requiredMatch = ...
    isfield(cfg, 'SelectedShell') && isequal(cfg.SelectedShell, selectedShell) && ...
    isfield(cfg, 'AddUserAccess') && isequal(logical(cfg.AddUserAccess), logical(addUserAccess)) && ...
    isfield(cfg, 'AddGatewayAccess') && isequal(logical(cfg.AddGatewayAccess), logical(addGatewayAccess)) && ...
    isfield(cfg, 'AddISL') && isequal(logical(cfg.AddISL), logical(addISL))

if ~requiredMatch
    tf = false;
    return;
end

if addISL && ~isfield(S, 'islAccesses')
    tf = false;
    return;
end
if addUserAccess && ~isfield(S, 'userAccesses')
    tf = false;
    return;
end
if addGatewayAccess && ~isfield(S, 'gatewayAccesses')
    tf = false;
end
end

function constData = readConstellationXml(xmlFile)
doc = xmlread(xmlFile);
root = doc.getDocumentElement;

constData = struct();
constData.NumberOfShells = round(getNodeValue(root, 'number_of_shells'));

shellTemplate = struct( ...
    'Altitude_km', NaN, ...
    'OrbitCycle_s', NaN, ...
    'Inclination_deg', NaN, ...
    'PhaseShift', NaN, ...
    'NumPlanes', NaN, ...
    'SatsPerPlane', NaN);

constData.Shells = repmat(shellTemplate, 1, constData.NumberOfShells);

for s = 1:constData.NumberOfShells
    shellTag = sprintf('shell%d', s);
    shellNodeList = doc.getElementsByTagName(shellTag);
    if shellNodeList.getLength == 0
        error('Missing <%s> in constellation XML.', shellTag);
    end

    shellNode = shellNodeList.item(0);
    constData.Shells(s).Altitude_km     = getNodeValue(shellNode, 'altitude');
    constData.Shells(s).OrbitCycle_s    = getOptionalNodeValue(shellNode, 'orbit_cycle', NaN);
    constData.Shells(s).Inclination_deg = getNodeValue(shellNode, 'inclination');
    constData.Shells(s).PhaseShift      = round(getNodeValue(shellNode, 'phase_shift'));
    constData.Shells(s).NumPlanes       = round(getNodeValue(shellNode, 'number_of_orbit'));
    constData.Shells(s).SatsPerPlane    = round(getNodeValue(shellNode, 'number_of_satellite_per_orbit'));
end
end

function data = readTerminalXml(xmlFile)
doc = xmlread(xmlFile);
root = doc.getDocumentElement;
children = root.getChildNodes;

data = struct( ...
    'Id', {}, ...
    'Latitude', {}, ...
    'Longitude', {}, ...
    'Name', {}, ...
    'AntennaCount', {}, ...
    'UplinkGHz', {}, ...
    'DownlinkGHz', {} );

idx = 0;
for i = 0:children.getLength-1
    node = children.item(i);
    if node.getNodeType ~= node.ELEMENT_NODE
        continue;
    end
    idx = idx + 1;
    data(idx).Id        = char(node.getNodeName);
    data(idx).Latitude  = getNodeValue(node, 'Latitude');
    data(idx).Longitude = getNodeValue(node, 'Longitude');
end
end

function [sampleTime, timeOfArrival, latDeg, lonDeg, altMeters] = readShellPositionsFromHdf5(hdf5File, selectedShell, fallbackSampleTime)
groupPath = sprintf('/position/shell%d', selectedShell);
groupInfo = h5info(hdf5File, groupPath);
datasetNames = {groupInfo.Datasets.Name};
sampleTime = readSampleTimeFromHdf5(hdf5File, fallbackSampleTime);
slotNumbers = zeros(size(datasetNames));

for k = 1:numel(datasetNames)
    token = regexp(datasetNames{k}, '\d+', 'match', 'once');
    if isempty(token)
        error('Invalid timeslot dataset name: %s', datasetNames{k});
    end
    slotNumbers(k) = str2double(token);
end

[~, order] = sort(slotNumbers);
datasetNames = datasetNames(order);

firstData = h5read(hdf5File, sprintf('%s/%s', groupPath, datasetNames{1}));
firstData = convertPositionDataset(firstData);
numTimeslots = numel(datasetNames);
numSatellites = size(firstData, 1);

latDeg = zeros(numTimeslots, numSatellites);
lonDeg = zeros(numTimeslots, numSatellites);
altMeters = zeros(numTimeslots, numSatellites);

lonCells = cell(numTimeslots, 1);
latCells = cell(numTimeslots, 1);
altCells = cell(numTimeslots, 1);
parfor t = 1:numTimeslots
    raw = h5read(hdf5File, sprintf('%s/%s', groupPath, datasetNames{t}));
    pos = convertPositionDataset(raw);
    if size(pos, 1) ~= numSatellites || size(pos, 2) ~= 3
        error('Unexpected position dataset size in %s/%s', groupPath, datasetNames{t});
    end
    lonCells{t} = pos(:, 1).';
    latCells{t} = pos(:, 2).';
    altCells{t} = pos(:, 3).' * 1000.0;
end
for t = 1:numTimeslots
    lonDeg(t, :) = lonCells{t};
    latDeg(t, :) = latCells{t};
    altMeters(t, :) = altCells{t};
end

timeOfArrival = (0:numTimeslots-1).' * sampleTime;
end

function islPairs = readISLPairsFromHdf5(hdf5File, selectedShell, totalSatellites)
groupPath = sprintf('/delay/shell%d', selectedShell);
groupInfo = h5info(hdf5File, groupPath);
datasetNames = {groupInfo.Datasets.Name};
slotNumbers = zeros(size(datasetNames));

for k = 1:numel(datasetNames)
    token = regexp(datasetNames{k}, '\d+', 'match', 'once');
    if isempty(token)
        error('Invalid timeslot dataset name: %s', datasetNames{k});
    end
    slotNumbers(k) = str2double(token);
end

[~, order] = sort(slotNumbers);
datasetNames = datasetNames(order);

pairMask = false(totalSatellites, totalSatellites);
pairMasks = cell(numel(datasetNames), 1);
parfor t = 1:numel(datasetNames)
    delayMatrix = double(h5read(hdf5File, sprintf('%s/%s', groupPath, datasetNames{t})));
    satBlock = delayMatrix(1:totalSatellites, 1:totalSatellites);
    pairMasks{t} = (satBlock > 0);
end
for t = 1:numel(datasetNames)
    pairMask = pairMask | pairMasks{t};
end

pairMask = triu(pairMask, 1);
[rows, cols] = find(pairMask);
islPairs = [rows, cols];
end

function terminalSatPairs = readTerminalSatPairsFromHdf5(hdf5File, selectedShell, totalSatellites, terminalStartIndex, numTerminals)
groupPath = sprintf('/delay/shell%d', selectedShell);
groupInfo = h5info(hdf5File, groupPath);
datasetNames = {groupInfo.Datasets.Name};
slotNumbers = zeros(size(datasetNames));

for k = 1:numel(datasetNames)
    token = regexp(datasetNames{k}, '\d+', 'match', 'once');
    if isempty(token)
        error('Invalid timeslot dataset name: %s', datasetNames{k});
    end
    slotNumbers(k) = str2double(token);
end

[~, order] = sort(slotNumbers);
datasetNames = datasetNames(order);

terminalSatPairs = false(totalSatellites, numTerminals);
if numTerminals == 0
    return;
end
terminalIndices = terminalStartIndex + (1:numTerminals);
terminalPairMasks = cell(numel(datasetNames), 1);
parfor t = 1:numel(datasetNames)
    delayMatrix = double(h5read(hdf5File, sprintf('%s/%s', groupPath, datasetNames{t})));
    terminalBlock = delayMatrix(1:totalSatellites, terminalIndices);
    terminalPairMasks{t} = (terminalBlock > 0);
end
for t = 1:numel(datasetNames)
    terminalSatPairs = terminalSatPairs | terminalPairMasks{t};
end
end

function sampleTime = readSampleTimeFromHdf5(hdf5File, fallbackSampleTime)
sampleTime = fallbackSampleTime;
try
    sampleTime = double(h5readatt(hdf5File, '/info', 'dT'));
catch
end
end

function minimumElevationDeg = readMinimumElevationFromHdf5(hdf5File)
minimumElevationDeg = NaN;
try
    minimumElevationDeg = double(h5readatt(hdf5File, '/info', 'minimum_elevation_deg'));
catch
end
end

function obj = createTrajectoryObject(sc, traj, name)
try
    obj = satellite(sc, traj, Name=name);
    return;
catch
end

try
    obj = satellite(sc, traj);
    try
        obj.Name = string(name);
    catch
    end
    return;
catch
end

obj = platform(sc, traj, Name=name);
end

function pos = convertPositionDataset(raw)
if isnumeric(raw)
    if ndims(raw) == 2 && size(raw, 1) == 3
        pos = double(raw).';
    else
        pos = double(raw);
    end
    return;
end

if isstring(raw)
    pos = str2double(raw);
    return;
end

if iscell(raw)
    pos = str2double(string(raw));
    return;
end

if ischar(raw)
    pos = str2double(string(cellstr(raw)));
    return;
end

error('Unsupported HDF5 dataset type for position data.');
end

function val = getNodeValue(parentNode, tagName)
txt = getNodeText(parentNode, tagName);
val = str2double(char(txt));
end

function val = getOptionalNodeValue(parentNode, tagName, defaultValue)
nodeList = parentNode.getElementsByTagName(tagName);
if nodeList.getLength == 0
    val = defaultValue;
else
    val = str2double(char(nodeList.item(0).getTextContent));
end
end

function txt = getNodeText(parentNode, tagName)
nodeList = parentNode.getElementsByTagName(tagName);
if nodeList.getLength == 0
    error('Missing tag <%s> in XML.', tagName);
end
txt = nodeList.item(0).getTextContent;
end
