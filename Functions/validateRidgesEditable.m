function outVal = validateRidgesEditable(outRidges, SandTable, varargin)
% {
%      Deanna Lanier 4.28.2026
% 
        % validateRidgesEditable (adapted from the validateRidgesInteractive. 
        %
        % Interactive validation + manual editing of ridge membership.
        %
        % Main workflow:
        %   - Shows all SAND peaks as gray dots when SandTable is provided.
        %   - Highlights one ridge at a time in yellow with black line.
        %   - Command window uses numeric choices to decide whether to accept, reject, or edit the ridge.
        %   - Edit mode lets you click directly on the plot to add, remove, or replace
        %     ridge points.
        %
        % Recommended:
        %   outVal = validateRidgesEditable( ...
        %       pHIT.outRidges, SandTable, ...
        %       'PeakTableField','Con_AF', ...
        %       'ppmCol','freq_ppm');
        %
        % Notes:
        %   - A ridge can have at most one peak from the same reordered sample.
        %   - In edit mode:
        %       Left-click gray/SAND peak area  = add or replace nearest SAND peak.
        %       Left-click existing ridge point = remove that point.
        %       Press Enter/Return             = finish editing current ridge.
        %   - After editing, choose Accept/Reject/Stop for the corrected ridge.

%% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('outRidges', @isstruct);
p.addOptional('SandTable', [], @(x) isempty(x) || istable(x) || isstruct(x));
p.addParameter('PeakTableField', 'Misaligned_Table', @(s)ischar(s)||isstring(s));
p.addParameter('ppmCol', 'Aligned_PPM', @(s)ischar(s)||isstring(s));
p.addParameter('AmpCol', 'amplitude', @(s)ischar(s)||isstring(s));
p.addParameter('ClickTolerancePPM', 0.005, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('ShowAllPeaks', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('MarkerSizeAll', 12, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('MarkerAlphaAll', 0.22, @(v)isnumeric(v)&&isscalar(v)&&v>=0&&v<=1);
p.addParameter('CurrentPointSize', 70, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('RemoveClickDistance', 0.018, @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.parse(outRidges, SandTable, varargin{:});
S = p.Results;

peakField = char(S.PeakTableField);
ppmCol    = char(S.ppmCol);
ampCol    = char(S.AmpCol);
SandTable = S.SandTable;

%% ---------------- Basic checks ----------------
if ~isfield(outRidges, 'ridgeSummaryTable') || isempty(outRidges.ridgeSummaryTable)
    error('outRidges.ridgeSummaryTable is missing or empty.');
end
if ~isfield(outRidges, 'ridgePointsTable') || isempty(outRidges.ridgePointsTable)
    error('outRidges.ridgePointsTable is missing or empty.');
end

sumTbl = outRidges.ridgeSummaryTable;
ptTbl  = outRidges.ridgePointsTable;

if ~ismember('RidgeID', sumTbl.Properties.VariableNames) || ...
   ~ismember('RidgeID', ptTbl.Properties.VariableNames)
    error('Both ridgeSummaryTable and ridgePointsTable must contain RidgeID.');
end

sampleCol = findSampleColumn(ptTbl);
if isempty(sampleCol)
    error('ridgePointsTable must contain SampleOrderId or SampleOrderId.');
end
if ~ismember('PeakPPM', ptTbl.Properties.VariableNames)
    error('ridgePointsTable must contain PeakPPM.');
end

% Standardize internally to SampleOrderId while preserving SampleOrderId if present.
if ~ismember('SampleOrderId', ptTbl.Properties.VariableNames)
    ptTbl.SampleOrderId = ptTbl.(sampleCol);
end
if ~ismember('SampleOrderId', ptTbl.Properties.VariableNames)
    ptTbl.SampleOrderId = ptTbl.SampleOrderId;
end

% Make sure useful columns exist.
if ~ismember('Accepted', sumTbl.Properties.VariableNames)
    sumTbl.Accepted = false(height(sumTbl),1);
else
    sumTbl.Accepted(:) = false;
end

if ~ismember('Accepted', ptTbl.Properties.VariableNames)
    ptTbl.Accepted = false(height(ptTbl),1);
else
    ptTbl.Accepted(:) = false;
end

if ~ismember('PeakNumber', ptTbl.Properties.VariableNames)
    ptTbl.PeakNumber = ones(height(ptTbl),1);
end

if ~ismember('OriginalSampleId', ptTbl.Properties.VariableNames)
    ptTbl.OriginalSampleId = nan(height(ptTbl),1);
end

if ~ismember('PeakRow', ptTbl.Properties.VariableNames)
    ptTbl.PeakRow = nan(height(ptTbl),1);
end

% Order mapping.
if isfield(outRidges, 'order') && ~isempty(outRidges.order)
    order = outRidges.order(:);
elseif isfield(outRidges, 'reorderOut') && isfield(outRidges.reorderOut,'order')
    order = outRidges.reorderOut.order(:);
else
    order = [];
end

% Load all SAND peak dots if available.
haveSand = ~isempty(SandTable);
peakCells = [];
allPeakDotTable = table();

if haveSand
    if isempty(order)
        error('SandTable editing requires outRidges.order or outRidges.reorderOut.order.');
    end
    peakCells = getPeakCells(SandTable, peakField);
    allPeakDotTable = buildAllPeakDotTable(peakCells, order, ppmCol, ampCol);
end

%% ---------------- Ridge order ----------------
nRidges = height(sumTbl);
if ismember('MinPPM', sumTbl.Properties.VariableNames)
    [~, idSort] = sort(sumTbl.MinPPM);
elseif ismember('RegionMinPPM', sumTbl.Properties.VariableNames)
    [~, idSort] = sort(sumTbl.RegionMinPPM);
else
    idSort = (1:nRidges)';
end
ridgeOrder = sumTbl.RidgeID(idSort);

%% ---------------- Figure setup ----------------
[fig, axBot] = setupValidationFigure(outRidges);
set(fig, 'Visible','on');
figure(fig);

% Draw all peak dots.
cla(axBot);
hold(axBot,'on');

if haveSand && S.ShowAllPeaks && ~isempty(allPeakDotTable)
    scatter(axBot, allPeakDotTable.PeakPPM, allPeakDotTable.SampleOrderId, ...
        S.MarkerSizeAll, [0.45 0.45 0.45], 'filled', ...
        'MarkerFaceAlpha', S.MarkerAlphaAll, ...
        'MarkerEdgeAlpha', S.MarkerAlphaAll);
else
    scatter(axBot, ptTbl.PeakPPM, ptTbl.SampleOrderId, ...
        S.MarkerSizeAll, [0.45 0.45 0.45], 'filled', ...
        'MarkerFaceAlpha', S.MarkerAlphaAll, ...
        'MarkerEdgeAlpha', S.MarkerAlphaAll);
end

set(axBot, 'XDir','reverse', 'YDir','normal');
xlabel(axBot, sprintf('%s (ppm)', ppmCol), 'Interpreter','none');
ylabel(axBot, 'Sample index (reordered)');
title(axBot, 'Ridge validation/editing');
grid(axBot,'on');
hold(axBot,'off');

%% ---------------- Print compact command guide once ----------------
fprintf('\nRIDGE VALIDATION OPTIONS\n');
fprintf('  1 = Accept current ridge\n');
fprintf('  2 = Reject current ridge\n');
fprintf('  3 = Edit current ridge manually\n');
fprintf('  4 = Stop and save progress\n\n');
fprintf('EDIT MODE\n');
fprintf('  Left-click an existing highlighted ridge point to remove it.\n');
fprintf('  Left-click a gray/SAND peak to add it to the current ridge.\n');
fprintf('  If the ridge already has a point for that sample, the clicked peak replaces it.\n');
fprintf('  Press Enter/Return while in the plot to finish editing that ridge.\n\n');

%% ---------------- Validation loop ----------------
for k = 1:numel(ridgeOrder)
    rid = ridgeOrder(k);

    iSum = find(sumTbl.RidgeID == rid, 1);
    if isempty(iSum)
        continue;
    end

    while true
        maskPts = (ptTbl.RidgeID == rid);

        % Highlight current ridge.
        [hLine, hPts, oldTitleStr] = drawCurrentRidge(axBot, ptTbl, sumTbl, iSum, rid, S);

        figure(fig);
        drawnow;

        fprintf('Ridge %d: ', rid);
        choiceNum = askNumericChoice(1:4, 1);

        deleteIfValid(hLine);
        deleteIfValid(hPts);
        if isgraphics(axBot)
            title(axBot, oldTitleStr);
        end
        drawnow;

        if choiceNum == 4
            acceptedIDs = sumTbl.RidgeID(sumTbl.Accepted);
            outVal = packOutput(outRidges, sumTbl, ptTbl, acceptedIDs);
            return;

        elseif choiceNum == 1
            sumTbl.Accepted(iSum) = true;
            ptTbl.Accepted(ptTbl.RidgeID == rid) = true;
            break;

        elseif choiceNum == 2
            sumTbl.Accepted(iSum) = false;
            ptTbl.Accepted(ptTbl.RidgeID == rid) = false;
            break;

        elseif choiceNum == 3
            if ~haveSand
                warning('Manual editing requires SandTable input. Continuing without edit.');
                continue;
            end

            [ptTbl, sumTbl] = editRidgeByClicking(ptTbl, sumTbl, iSum, rid, axBot, ...
                allPeakDotTable, peakCells, order, ppmCol, ampCol, S);
            % Go back to validation choice after editing.
        end
    end
end

acceptedIDs = sumTbl.RidgeID(sumTbl.Accepted);

%% ---------------- Pack outputs ----------------
outVal = packOutput(outRidges, sumTbl, ptTbl, acceptedIDs);

end

%% =====================================================================
function [ptTbl, sumTbl] = editRidgeByClicking(ptTbl, sumTbl, iSum, rid, axBot, ...
    allPeakDotTable, peakCells, order, ppmCol, ampCol, S)

fig = ancestor(axBot, 'figure');

fprintf('\nEditing Ridge %d.\n', rid);
fprintf('Click highlighted point = remove. Click gray/SAND peak = add/replace. Press Enter in plot when done.\n');

while true
    [hLine, hPts, oldTitleStr] = drawCurrentRidge(axBot, ptTbl, sumTbl, iSum, rid, S);
    drawnow;

    figure(fig);
    axes(axBot); %#ok<LAXES>

    try
        [xClick, yClick, button] = ginput(1);
    catch
        button = [];
        xClick = [];
        yClick = [];
    end

    deleteIfValid(hLine);
    deleteIfValid(hPts);
    if isgraphics(axBot), title(axBot, oldTitleStr); end

    % Enter/return ends edit mode in MATLAB ginput.
    if isempty(button) || isempty(xClick)
        fprintf('Done editing Ridge %d.\n\n', rid);
        [ptTbl, sumTbl] = updateSummaryForRidge(ptTbl, sumTbl, iSum, rid);
        break;
    end

    if button ~= 1
        fprintf('Done editing Ridge %d.\n\n', rid);
        [ptTbl, sumTbl] = updateSummaryForRidge(ptTbl, sumTbl, iSum, rid);
        break;
    end

    % First check if click is close to an existing ridge point -> remove.
    [isExisting, rowToRemove] = nearestExistingRidgePoint(ptTbl, rid, xClick, yClick, axBot, S.RemoveClickDistance);
    if isExisting
        fprintf('Removed point from Ridge %d: sample %d, ppm %.5f.\n', ...
            rid, ptTbl.SampleOrderId(rowToRemove), ptTbl.PeakPPM(rowToRemove));
        ptTbl(rowToRemove,:) = [];
        [ptTbl, sumTbl] = updateSummaryForRidge(ptTbl, sumTbl, iSum, rid);
        continue;
    end

    % Otherwise add/replace nearest SAND peak.
    sOrd = round(yClick);
    sOrd = max(1, min(numel(order), sOrd));

    [newRow, ok, msg] = makeRidgePointRowFromClick(ptTbl, rid, sOrd, xClick, ...
        peakCells, order, ppmCol, ampCol, S.ClickTolerancePPM);

    if ~ok
        fprintf('Could not add point: %s\n', msg);
        continue;
    end

    existing = find(ptTbl.RidgeID == rid & ptTbl.SampleOrderId == sOrd);
    if ~isempty(existing)
        fprintf('Replaced existing point for Ridge %d, sample %d.\n', rid, sOrd);
        ptTbl(existing,:) = [];
    else
        fprintf('Added point to Ridge %d: sample %d, ppm %.5f.\n', ...
            rid, sOrd, newRow.PeakPPM);
    end

    ptTbl = [ptTbl; newRow];
    [ptTbl, sumTbl] = updateSummaryForRidge(ptTbl, sumTbl, iSum, rid);
end
end

%% =====================================================================
function [tf, rowToRemove] = nearestExistingRidgePoint(ptTbl, rid, xClick, yClick, axBot, tolNorm)
rows = find(ptTbl.RidgeID == rid);
tf = false;
rowToRemove = NaN;

if isempty(rows)
    return;
end

samp = ptTbl.SampleOrderId(rows);
ppmv = ptTbl.PeakPPM(rows);

xl = xlim(axBot);
yl = ylim(axBot);
dx = abs(ppmv - xClick) ./ max(abs(diff(xl)), eps);
dy = abs(samp - yClick) ./ max(abs(diff(yl)), eps);
dist = sqrt(dx.^2 + dy.^2);

[dmin, j] = min(dist);
if dmin <= tolNorm
    tf = true;
    rowToRemove = rows(j);
end
end

%% =====================================================================
function [newRow, ok, msg] = makeRidgePointRowFromClick(ptTbl, rid, sOrd, xClick, ...
    peakCells, order, ppmCol, ampCol, tolPPM)

ok = false;
msg = '';
newRow = ptTbl(1,:);
newRow(1,:) = [];

origId = order(sOrd);
if origId < 1 || origId > numel(peakCells)
    msg = 'Clicked sample is outside the order mapping.';
    return;
end

Ti = peakCells{origId};
if isempty(Ti) || ~istable(Ti) || ~ismember(ppmCol, Ti.Properties.VariableNames)
    msg = sprintf('No usable peak table or ppm column for original sample %d.', origId);
    return;
end

ppmv = Ti.(ppmCol)(:);
if isempty(ppmv)
    msg = sprintf('No peaks in original sample %d.', origId);
    return;
end

[d, id] = min(abs(ppmv - xClick));
if isempty(d) || d > tolPPM
    msg = sprintf('Nearest peak is %.5f ppm away, greater than ClickTolerancePPM %.5f.', d, tolPPM);
    return;
end

newRow = makeDefaultRowLike(ptTbl);

for v = 1:numel(ptTbl.Properties.VariableNames)
    nm = ptTbl.Properties.VariableNames{v};
    switch nm
        case 'RidgeID'
            newRow.(nm) = rid;
        case {'SampleOrderId','SampleOrderId'}
            newRow.(nm) = sOrd;
        case 'OriginalSampleId'
            newRow.(nm) = origId;
        case 'PeakRow'
            newRow.(nm) = id;
        case 'PeakPPM'
            newRow.(nm) = ppmv(id);
        case 'ReconPPM'
            newRow.(nm) = ppmv(id);
        case 'MapDistancePPM'
            newRow.(nm) = abs(ppmv(id)-xClick);
        case 'Amplitude'
            if ismember(ampCol, Ti.Properties.VariableNames)
                newRow.(nm) = Ti.(ampCol)(id);
            else
                newRow.(nm) = NaN;
            end
        case 'Accepted'
            newRow.(nm) = false;
        case 'Inlier'
            newRow.(nm) = true;
        case 'PeakNumber'
            newRow.(nm) = 1;
        case 'Source'
            newRow.(nm) = string('manual_add_validation');
        case 'Mapped'
            newRow.(nm) = true;
    end
end

ok = true;
end

%% =====================================================================
function row = makeDefaultRowLike(T)
row = T(1,:);
for v = 1:numel(T.Properties.VariableNames)
    nm = T.Properties.VariableNames{v};
    col = T.(nm);
    if islogical(col)
        row.(nm) = false;
    elseif isnumeric(col)
        row.(nm) = NaN;
    elseif isstring(col)
        row.(nm) = "";
    elseif iscellstr(col) || iscell(col)
        row.(nm) = {''};
    elseif iscategorical(col)
        row.(nm) = categorical(missing);
    else
        try
            row.(nm) = missing;
        catch
        end
    end
end
end

%% =====================================================================
function [ptTbl, sumTbl] = updateSummaryForRidge(ptTbl, sumTbl, iSum, rid)
mask = ptTbl.RidgeID == rid;
if ~any(mask)
    vars = sumTbl.Properties.VariableNames;
    if ismember('NumPoints', vars), sumTbl.NumPoints(iSum) = 0; end
    if ismember('NumSamples', vars), sumTbl.NumSamples(iSum) = 0; end
    if ismember('CoverageFrac', vars), sumTbl.CoverageFrac(iSum) = 0; end
    if ismember('MappedCoverage', vars), sumTbl.MappedCoverage(iSum) = 0; end
    return;
end

samp = ptTbl.SampleOrderId(mask);
ppmv = ptTbl.PeakPPM(mask);

if ismember('NumPoints', sumTbl.Properties.VariableNames)
    sumTbl.NumPoints(iSum) = sum(mask);
end
if ismember('NumSamples', sumTbl.Properties.VariableNames)
    sumTbl.NumSamples(iSum) = numel(unique(samp));
end
if ismember('MinPPM', sumTbl.Properties.VariableNames)
    sumTbl.MinPPM(iSum) = min(ppmv);
end
if ismember('MaxPPM', sumTbl.Properties.VariableNames)
    sumTbl.MaxPPM(iSum) = max(ppmv);
end
if ismember('RegionMinPPM', sumTbl.Properties.VariableNames)
    sumTbl.RegionMinPPM(iSum) = min(ppmv);
end
if ismember('RegionMaxPPM', sumTbl.Properties.VariableNames)
    sumTbl.RegionMaxPPM(iSum) = max(ppmv);
end
if ismember('NMappedSandPoints', sumTbl.Properties.VariableNames)
    sumTbl.NMappedSandPoints(iSum) = sum(mask);
end
if ismember('CoverageFrac', sumTbl.Properties.VariableNames)
    span = max(samp) - min(samp) + 1;
    sumTbl.CoverageFrac(iSum) = numel(unique(samp)) / max(span,1);
end
if ismember('MappedCoverage', sumTbl.Properties.VariableNames)
    span = max(samp) - min(samp) + 1;
    sumTbl.MappedCoverage(iSum) = numel(unique(samp)) / max(span,1);
end
end

%% =====================================================================
function [hLine, hPts, oldTitleStr] = drawCurrentRidge(axBot, ptTbl, sumTbl, iSum, rid, S)
oldTitleObj = get(axBot,'Title');
if isgraphics(oldTitleObj) && isprop(oldTitleObj,'String')
    oldTitleStr = oldTitleObj.String;
else
    oldTitleStr = '';
end

maskPts = ptTbl.RidgeID == rid;
sampR = ptTbl.SampleOrderId(maskPts);
ppmR  = ptTbl.PeakPPM(maskPts);

if isempty(sampR)
    hLine = gobjects(1);
    hPts = gobjects(1);
    title(axBot, sprintf('Ridge %d | no points', rid), 'FontWeight','bold');
    return;
end

[sampR, ord] = sort(sampR);
ppmR = ppmR(ord);

figure(ancestor(axBot,'figure'));
axes(axBot); %#ok<LAXES>
hold(axBot,'on');
hLine = plot(axBot, ppmR, sampR, '-', 'Color','k', 'LineWidth', 2.5);
hPts = scatter(axBot, ppmR, sampR, S.CurrentPointSize, ...
    'MarkerFaceColor','y', 'MarkerEdgeColor','k', 'LineWidth',1.2);
hold(axBot,'off');

txtTitle = sprintf('Ridge %d | samples=%d, points=%d, span=[%.4f, %.4f] ppm', ...
    rid, ...
    getFieldOrDefault(sumTbl,iSum,'NumSamples',numel(unique(sampR))), ...
    getFieldOrDefault(sumTbl,iSum,'NumPoints',numel(sampR)), ...
    getFieldOrDefault(sumTbl,iSum,'MinPPM',min(ppmR)), ...
    getFieldOrDefault(sumTbl,iSum,'MaxPPM',max(ppmR)));

title(axBot, txtTitle, 'FontWeight','bold');
drawnow;
end

%% =====================================================================
function [fig, axBot] = setupValidationFigure(outRidges)
axBot = [];
fig = [];

if isfield(outRidges, 'axBot') && isgraphics(outRidges.axBot)
    axBot = outRidges.axBot;
    fig = ancestor(axBot, 'figure');
end

if isempty(fig) || ~isgraphics(fig)
    if isfield(outRidges, 'fig') && isgraphics(outRidges.fig)
        fig = outRidges.fig;
        ax = findobj(fig, 'Type','axes');
        if isempty(ax)
            axBot = axes('Parent',fig);
        else
            axBot = ax(1);
        end
    else
        fig = figure('Color','w','Name','Ridge validation and editing');
        axBot = axes('Parent',fig);
    end
end
end

%% =====================================================================
function allPeakDotTable = buildAllPeakDotTable(peakCells, order, ppmCol, ampCol)
rows = {};
for sOrd = 1:numel(order)
    origId = order(sOrd);
    if origId < 1 || origId > numel(peakCells)
        continue;
    end
    Ti = peakCells{origId};
    if isempty(Ti) || ~istable(Ti) || ~ismember(ppmCol, Ti.Properties.VariableNames)
        continue;
    end

    ppmv = Ti.(ppmCol)(:);
    if ismember(ampCol, Ti.Properties.VariableNames)
        ampv = Ti.(ampCol)(:);
    else
        ampv = nan(height(Ti),1);
    end

    for r = 1:numel(ppmv)
        rows(end+1,:) = {sOrd, origId, r, ppmv(r), ampv(r)}; %#ok<AGROW>
    end
end

if isempty(rows)
    allPeakDotTable = table();
else
    allPeakDotTable = cell2table(rows, ...
        'VariableNames', {'SampleOrderId','OriginalSampleId','PeakRow','PeakPPM','Amplitude'});
end
end

%% =====================================================================
function peakCells = getPeakCells(SandTable, fieldName)
if istable(SandTable)
    if ~ismember(fieldName, SandTable.Properties.VariableNames)
        error('SandTable does not contain variable "%s".', fieldName);
    end
    peakCells = SandTable.(fieldName);
elseif isstruct(SandTable)
    if ~isfield(SandTable, fieldName)
        error('SandTable struct does not contain field "%s".', fieldName);
    end
    peakCells = SandTable.(fieldName);
else
    error('SandTable must be a table or struct.');
end
if ~iscell(peakCells)
    error('PeakTableField must be a cell array of per-sample tables.');
end
end

%% =====================================================================
function sampleCol = findSampleColumn(ptTbl)
if ismember('SampleOrderId', ptTbl.Properties.VariableNames)
    sampleCol = 'SampleOrderId';
elseif ismember('SampleOrderId', ptTbl.Properties.VariableNames)
    sampleCol = 'SampleOrderId';
else
    sampleCol = '';
end
end

%% =====================================================================
function val = getFieldOrDefault(tbl, id, varName, defaultVal)
if ismember(varName, tbl.Properties.VariableNames)
    v = tbl.(varName);
    if id <= numel(v)
        val = v(id);
        return;
    end
end
val = defaultVal;
end

%% =====================================================================
function deleteIfValid(h)
if isempty(h), return; end
try
    h = h(isvalid(h));
    if ~isempty(h), delete(h); end
catch
end
end

%% =====================================================================
function choiceNum = askNumericChoice(validChoices, defaultChoice)
while true
    resp = input('Select option: ', 's');

    if isempty(strtrim(resp))
        choiceNum = defaultChoice;
        return;
    end

    val = str2double(resp);
    if ~isnan(val) && ismember(round(val), validChoices)
        choiceNum = round(val);
        return;
    end

    fprintf('Invalid choice. Select one of: ');
    fprintf('%d ', validChoices);
    fprintf('\n');
end
end

%% =====================================================================
function outVal = packOutput(outRidges, sumTbl, ptTbl, acceptedIDs)
outVal = outRidges;
outVal.ridgeSummaryTable = sumTbl;
outVal.ridgePointsTable = ptTbl;
outVal.acceptedRidgeIDs = acceptedIDs;
end
