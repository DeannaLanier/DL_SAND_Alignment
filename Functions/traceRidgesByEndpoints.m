function outRidgesManual = traceRidgesByEndpoints(ppm, reorderOut, SandTable, varargin)
% traceRidgesByEndpoints
%
% Interactive Sicong-style endpoint/window ridge detector adapted for SAND.
%
% This version is designed to match the older Sicong responding-region logic
% more closely:
%
%   User clicks two endpoints -> this defines a ppm range and sample span.
%   User enters NRidges -> expected number of ridge tracks/signals in the range.
%
%   For the selected window:
%       1) Detect the highest local maxima in each spectrum within the bounds.
%       2) Run DBSCAN on [ppm position, sample number].
%       3) For NRidges == 1:
%             use the largest DBSCAN cluster as the on-curve peak set,
%             fit a polynomial curve,
%             use a prediction interval to re-pick one peak per sample.
%       4) For NRidges > 1:
%             first try to separate tracks by DBSCAN clusters;
%             if DBSCAN merges close curves, fall back to ordered peak columns,
%             then fit/refine each ridge separately.
%       5) Map final reconstructed ridge ppm positions back to SAND peak rows.
%
% Output is compatible with:
%       outVal = validateRidgesInteractive(...);
%       outAligned = alignAndReconstructValidatedRidges(...);
%
% Recommended call:
%   outRidgesManual = traceRidgesByEndpoints( ...
%       ppm, pHIT.out, SandTable, ...
%       'PeakTableField','Peaks_Condensed_RR', ...
%       'ppmCol','freq_ppm', ...
%       'InterFactor',10000, ...
%       'MapTolerancePPM',0.0025);

%% ---------------- Parse inputs ----------------
p = inputParser;
p.addRequired('ppm', @(v)isnumeric(v)&&isvector(v));
p.addRequired('reorderOut', @isstruct);
p.addRequired('SandTable', @(x)istable(x)||isstruct(x));

p.addParameter('PeakTableField', 'Peaks_Condensed_RR', @(s)ischar(s)||isstring(s));
p.addParameter('ppmCol', 'freq_ppm', @(s)ischar(s)||isstring(s));
p.addParameter('AmpCol', 'amplitude', @(s)ischar(s)||isstring(s));
p.addParameter('DecayCol', 'decay_hz', @(s)ischar(s)||isstring(s));

p.addParameter('InterFactor', 500, @(v)isnumeric(v)&&isscalar(v));
p.addParameter('MapTolerancePPM', 0.0025, @(v)isnumeric(v)&&isscalar(v));
p.addParameter('RangePaddingPPM', 0.000, @(v)isnumeric(v)&&isscalar(v));
p.addParameter('PromptForNumPeaks', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('DefaultNumPeaks', 1, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
p.addParameter('PredictionLevel', 0.95, @(v)isnumeric(v)&&isscalar(v)&&v>0&&v<1);
p.addParameter('PolyOrder', 3, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
p.addParameter('DBSCANMinPts', 3, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
p.addParameter('UseLargestClusterForSingle', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('ShowDBSCANFigure', false, @(b)islogical(b)&&isscalar(b));

p.addParameter('MakeFigure', true, @(b)islogical(b)&&isscalar(b));
p.addParameter('FigureVisible', 'on', @(s)ischar(s)||isstring(s));
p.addParameter('EnableZoomBeforeSelect', true, @(b)islogical(b)&&isscalar(b));

p.parse(ppm, reorderOut, SandTable, varargin{:});
S = p.Results;

ppm = ppm(:)';
X_reo = reorderOut.X_reordered;
order = reorderOut.order(:);
[nSamples, nPoints] = size(X_reo);

if numel(ppm) ~= nPoints
    error('ppm length must match columns of reorderOut.X_reordered.');
end
if ~S.MakeFigure
    error('Interactive tracing requires MakeFigure=true.');
end

peakCells = getPeakCells(SandTable, char(S.PeakTableField));
sandCache = cacheSandPeaks(peakCells, char(S.ppmCol), char(S.AmpCol), char(S.DecayCol));

%% ---------------- Figure ----------------
fig = figure('Color','w', 'Visible', S.FigureVisible, ...
    'Name','Sicong-style SAND endpoint/window ridge tracer', ...
    'NumberTitle','off');

tl = tiledlayout(fig, 2, 1);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';
axTop = nexttile(tl, 1);
axBot = nexttile(tl, 2);

plotBaseViews(axTop, axBot, ppm, X_reo, order, sandCache, S.InterFactor);
title(axTop, {'Stacked reordered spectra', ...
    'Zoom/pan first. Press SPACE/ENTER or Ready, then click START and END of the ppm/sample window.'});

setappdata(fig, 'EndpointTracerAction', '');
readyBtn = uicontrol(fig, 'Style','pushbutton', ...
    'Units','normalized', 'Position',[0.73 0.955 0.25 0.035], ...
    'String','Ready to click endpoints', ...
    'FontWeight','bold', ...
    'Callback', @(~,~) setReadyAndResume(fig));

doneBtn = uicontrol(fig, 'Style','pushbutton', ...
    'Units','normalized', 'Position',[0.60 0.955 0.12 0.035], ...
    'String','Done', ...
    'Callback', @(~,~) setDoneAndResume(fig));

set(fig, 'WindowKeyPressFcn', @(src,event) endpointKeyPress(src,event));

ridgeSummaryRows = {};
ridgePointRows = {};
reconPointRows = {};
nextRidgeID = 1;

%% ---------------- Interactive loop ----------------
while isvalid(fig)
    figure(fig);
    axes(axTop); %#ok<LAXES>

    fprintf('\nRidge set starting at %d: zoom/pan first, then press SPACE/ENTER or Ready.\n', nextRidgeID);

    setappdata(fig, 'EndpointTracerAction', '');
    if isvalid(readyBtn), set(readyBtn, 'Enable','on', 'Visible','on'); end
    if isvalid(doneBtn),  set(doneBtn,  'Enable','on', 'Visible','on'); end

    if S.EnableZoomBeforeSelect
        try
            zoom(fig, 'on');
            pan(fig, 'off');
        catch
        end
        uiwait(fig);
        if ~isvalid(fig), break; end
        action = getappdata(fig, 'EndpointTracerAction');
        if strcmp(action, 'done')
            break;
        end
        try
            zoom(fig, 'off');
            pan(fig, 'off');
        catch
        end
    end

    if isvalid(readyBtn), set(readyBtn, 'Enable','off'); end

    fprintf('Click START and END of the ppm/sample window.\n');
    try
        [xClick, yClick, button] = ginput(2);
    catch
        break;
    end

    if ~isvalid(fig), break; end
    if numel(xClick) < 2 || any(button ~= 1)
        break;
    end

    % Keep existing y-click behavior: y coordinate determines sample range.
    sClick = round(yClick ./ S.InterFactor + 1);
    sClick = max(1, min(nSamples, sClick));
    ppmClick = xClick(:);
    sClick = sClick(:);

    NRidges = S.DefaultNumPeaks;
    if S.PromptForNumPeaks
        answer = inputdlg({'Expected number of ridge tracks/signals in this selected window:'}, ...
            'Expected ridges', [1 60], {num2str(S.DefaultNumPeaks)});
        if isempty(answer)
            fprintf('Selection canceled.\n');
            continue;
        end
        NRidges = str2double(answer{1});
        if isnan(NRidges) || NRidges < 1
            NRidges = S.DefaultNumPeaks;
        end
        NRidges = round(NRidges);
    end

    % Show selected window, not a forced line.
    hold(axTop, 'on');
    hAnchor = plot(axTop, xClick, yClick, 'ks', ...
        'LineWidth', 1.2, 'MarkerFaceColor','y', 'MarkerSize', 8);
    yl = ylim(axTop);
    xlo = min(ppmClick); xhi = max(ppmClick);
    hPatch = patch(axTop, [xlo xhi xhi xlo], [yl(1) yl(1) yl(2) yl(2)], ...
        [1 0.9 0.55], 'FaceAlpha', 0.10, 'EdgeColor', [0.95 0.55 0], 'LineWidth', 1.0);
    hold(axTop, 'off');
    drawnow;

    try
        [trace, diagInfo] = traceOneWindowSicongExact(ppm, X_reo, ppmClick, sClick, NRidges, S);
    catch ME
        warning('Tracing failed for this selected window: %s', ME.message);
        trace = table();
        diagInfo = struct('RMSE', NaN, 'Notes', string(ME.message));
    end

    if isempty(trace) || height(trace) == 0
        deleteIfValid(hAnchor); deleteIfValid(hPatch);
        msgbox('No ridge traced for this selected range. Accepted ridges so far are still saved.', ...
            'No ridge traced', 'warn');
        choice = questdlg('No ridge was detected. What would you like to do?', ...
            'No ridge detected', ...
            'Try another range','Finish and save','Try another range');
        if isempty(choice) || strcmp(choice,'Finish and save')
            break;
        else
            continue;
        end
    end

    try
        mapped = mapTraceToSand(trace, order, sandCache, S.MapTolerancePPM);
    catch ME
        warning('SAND mapping failed for this traced ridge: %s', ME.message);
        mapped = makeUnmappedTable(trace, order);
    end

    peakNums = unique(trace.PeakNumber);
    peakNums = peakNums(:)';
    nTracksDetected = numel(peakNums);
    nMapped = sum(mapped.Mapped);

    fprintf('Window traced: %d reconstructed points; %d mapped SAND points; %d track(s).\n', ...
        height(trace), nMapped, nTracksDetected);

    % Preview each track with future RidgeID color.
    hold(axTop, 'on');
    hTraceTop = gobjects(0);
    for kk = 1:nTracksDetected
        ppn = peakNums(kk);
        ridPreview = nextRidgeID + kk - 1;
        c = ridgeColor(ridPreview);

        rows = trace.PeakNumber == ppn;
        tr = sortrows(trace(rows,:), 'SampleOrderId');
        yPlot = X_reo(sub2ind(size(X_reo), tr.SampleOrderId, tr.PPMIndex)) + ...
            (tr.SampleOrderId-1)*S.InterFactor;
        hTraceTop(end+1) = plot(axTop, tr.PeakPPM_Recon, yPlot, 'o-', 'Color', c, ...
            'MarkerFaceColor', c, 'MarkerEdgeColor','k', 'LineWidth', 1.3, 'MarkerSize', 5); %#ok<AGROW>
    end
    hold(axTop, 'off');

    hold(axBot, 'on');
    hTraceBot = gobjects(0);
    for kk = 1:nTracksDetected
        ppn = peakNums(kk);
        ridPreview = nextRidgeID + kk - 1;
        c = ridgeColor(ridPreview);

        if any(mapped.Mapped & mapped.PeakNumber == ppn)
            rows = mapped.Mapped & mapped.PeakNumber == ppn;
            mt = sortrows(mapped(rows,:), 'SampleOrderId');
            hTraceBot(end+1) = plot(axBot, mt.PeakPPM, mt.SampleOrderId, 'o-', ...
                'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor','k', ...
                'LineWidth', 1.3, 'MarkerSize', 5); %#ok<AGROW>
        end

        if any(~mapped.Mapped & mapped.PeakNumber == ppn)
            rows = ~mapped.Mapped & mapped.PeakNumber == ppn;
            hTraceBot(end+1) = plot(axBot, trace.PeakPPM_Recon(rows), trace.SampleOrderId(rows), ...
                'x--', 'Color', c, 'LineWidth', 1.0); %#ok<AGROW>
        end
    end
    hold(axBot, 'off');
    drawnow;

    statusText = sprintf(['Selected window traced using Sicong-style region logic.\n\n', ...
        'Requested ridge tracks: %d\nDetected ridge tracks: %d\nReconstructed points: %d\nMapped SAND points: %d\n\n', ...
        'Accept these as %d separate ridge(s)?'], ...
        NRidges, nTracksDetected, height(trace), nMapped, nTracksDetected);

    choice = questdlg(statusText, sprintf('Ridge set starting at %d', nextRidgeID), ...
        'Accept + next','Reject + next','Accept + done','Accept + next');

    if isempty(choice)
        choice = 'Accept + next';
    end

    acceptThis = strcmp(choice, 'Accept + next') || strcmp(choice, 'Accept + done');

    if acceptThis
        for kk = 1:nTracksDetected
            ppn = peakNums(kk);
            thisRidgeID = nextRidgeID;

            tr = trace(trace.PeakNumber == ppn,:);
            mp = mapped(mapped.PeakNumber == ppn,:);

            nSpanTrack = max(tr.SampleOrderId) - min(tr.SampleOrderId) + 1;
            nMappedTrack = sum(mp.Mapped);
            covTrack = nMappedTrack / max(nSpanTrack, 1);

            ridgeSummaryRows(end+1,:) = {thisRidgeID, min(tr.PeakPPM_Recon), max(tr.PeakPPM_Recon), ...
                min(tr.SampleOrderId), max(tr.SampleOrderId), height(tr), nMappedTrack, covTrack, ...
                diagInfo.RMSE, ppmClick(1), sClick(1), ppmClick(2), sClick(2), 1, ...
                'endpoint_window_exact_region'}; %#ok<AGROW>

            for r = 1:height(tr)
                reconPointRows(end+1,:) = {thisRidgeID, tr.SampleOrderId(r), tr.SampleOrderId(r), 1, ...
                    tr.PeakPPM_Recon(r), tr.PPMIndex(r), tr.ReconIntensity(r), ...
                    tr.FitPPM(r), tr.ResidualPPM(r), tr.ClusterLabel(r)}; %#ok<AGROW>
            end

            for r = 1:height(mp)
                if mp.Mapped(r)
                    ridgePointRows(end+1,:) = {thisRidgeID, mp.SampleOrderId(r), mp.SampleOrderId(r), ...
                        mp.OriginalSampleId(r), 1, mp.PeakRow(r), mp.PeakPPM(r), mp.ReconPPM(r), ...
                        mp.MapDistancePPM(r), mp.Amplitude(r), mp.DecayHz(r)}; %#ok<AGROW>
                end
            end

            nextRidgeID = nextRidgeID + 1;
        end
    else
        fprintf('Rejected ridge set. Removing preview overlay.\n');
        deleteIfValid(hAnchor); deleteIfValid(hPatch); deleteIfValid(hTraceTop); deleteIfValid(hTraceBot);
    end

    if strcmp(choice, 'Accept + done')
        break;
    end
end

ridgeSummaryTable = cell2table_safe(ridgeSummaryRows, {'RidgeID','RegionMinPPM','RegionMaxPPM', ...
    'StartSampleOrderId','EndSampleOrderId','NReconPoints','NMappedSandPoints','MappedCoverage', ...
    'FitRMSE','Click1PPM','Click1SampleOrderId','Click2PPM','Click2SampleOrderId','NumPeaksExpected','Source'});

ridgePointsTable = cell2table_safe(ridgePointRows, {'RidgeID','SampleOrderId','SampleOrderIdx','OriginalSampleId', ...
    'PeakNumber','PeakRow','PeakPPM','ReconPPM','MapDistancePPM','Amplitude','DecayHz'});

reconRidgePointsTable = cell2table_safe(reconPointRows, {'RidgeID','SampleOrderId','SampleOrderIdx','PeakNumber', ...
    'PeakPPM_Recon','PPMIndex','ReconIntensity','FitPPM','ResidualPPM','ClusterLabel'});

outRidgesManual = struct();
outRidgesManual.ridgeSummaryTable = ridgeSummaryTable;
outRidgesManual.ridgePointsTable = ridgePointsTable;
outRidgesManual.reconRidgePointsTable = reconRidgePointsTable;
outRidgesManual.order = order;
outRidgesManual.reorderOut = reorderOut;
outRidgesManual.ppm = ppm;
outRidgesManual.settings = S;

if height(ridgeSummaryTable) > 0
    outRidgesManual.acceptedRidgeIDs = unique(ridgeSummaryTable.RidgeID);
else
    outRidgesManual.acceptedRidgeIDs = [];
end

fprintf('\nSicong-style endpoint/window tracing complete: %d accepted ridge(s).\n', height(ridgeSummaryTable));
end

%% =====================================================================
function [trace, diagInfo] = traceOneWindowSicongExact(ppm, X_reo, ppmClick, sClick, NRidges, S)
[nSamples, ~] = size(X_reo);

sStart = min(sClick);
sEnd   = max(sClick);
if sStart == sEnd
    sStart = max(1, sStart-3);
    sEnd   = min(nSamples, sEnd+3);
end

ppmMin = min(ppmClick) - S.RangePaddingPPM;
ppmMax = max(ppmClick) + S.RangePaddingPPM;

idx = find(ppm >= ppmMin & ppm <= ppmMax);
if isempty(idx)
    trace = table();
    diagInfo = struct('RMSE', NaN);
    return;
end

sampleIdx = (sStart:sEnd)';
sampleLocalN = numel(sampleIdx);
y1 = (1:sampleLocalN)';

XRs = X_reo(sampleIdx, idx);
ppms = ppm(idx).';

% ---------- Sicong Step a: Find peaks in user range ----------
[locV, ~] = localFindPeaksSicong(XRs, 'rspd', NRidges, []);
Kfound = size(locV,2);
if Kfound == 0
    trace = table();
    diagInfo = struct('RMSE', NaN);
    return;
end
if Kfound < NRidges
    NRidges = Kfound;
end

x1 = ppms(locV);  % n x Kfound ppm peak positions

% ---------- Sicong Step b: DBSCAN clustering ----------
try
    ClsMat = [x1(:), repmat(y1, size(x1,2), 1)];
    cls = localSicongCluster(ClsMat, x1, y1, sampleLocalN, ternary(S.ShowDBSCANFigure,'on','off'), S.DBSCANMinPts);
catch ME
    warning('Cluster error: %s. Using ordered peak columns.', ME.message);
    cls = repmat(1:size(x1,2), sampleLocalN, 1);
end
if isvector(cls)
    cls = reshape(cls, size(x1));
end

if NRidges == 1
    [trace, rmse] = processSingleRidge(ppm, X_reo, sampleIdx, ppms, XRs, x1, cls, y1, idx, S);
else
    [trace, rmse] = processMultiRidge(ppm, X_reo, sampleIdx, ppms, XRs, x1, cls, y1, idx, NRidges, S);
end

diagInfo = struct('RMSE', rmse);
end

%% =====================================================================
function [trace, rmse] = processSingleRidge(ppm, X_reo, sampleIdx, ppms, XRs, x1, cls, y1, idx, S)
sampleLocalN = numel(sampleIdx);

% Use largest non-noise cluster as on-curve points.
clsVec = cls(:);
xVec = x1(:);
yMat = repmat(y1, 1, size(x1,2));
yVec = yMat(:);

valid = isfinite(xVec) & isfinite(yVec);
xVec = xVec(valid);
yVec = yVec(valid);
clsVec = clsVec(valid);

labs = unique(clsVec);
labs(labs <= 0) = [];

if isempty(labs)
    xFit = xVec;
    yFit = yVec;
    labUse = 1;
else
    counts = arrayfun(@(z)sum(clsVec==z), labs);
    [~,m] = max(counts);
    labUse = labs(m);
    on = clsVec == labUse;
    xFit = xVec(on);
    yFit = yVec(on);
end

if numel(xFit) < 2
    trace = table();
    rmse = NaN;
    return;
end

[fitFcn, rmse] = localPolyFitFunction(yFit, xFit, S.PolyOrder);
[predCenter, predLow, predHigh] = localPredInterval(fitFcn, y1, yFit, xFit, S.PredictionLevel);

rowsOut = repickOneTrackWithinIntervals(ppm, X_reo, sampleIdx, predCenter, predLow, predHigh, ...
    min(ppms), max(ppms), 1, labUse);

trace = cell2table_safe(rowsOut, {'SampleOrderId','SampleOrderIdx','PeakNumber','PeakPPM_Recon','PPMIndex', ...
    'ReconIntensity','FitPPM','ResidualPPM','ClusterLabel'});
end

%% =====================================================================
function [trace, rmse] = processMultiRidge(ppm, X_reo, sampleIdx, ppms, XRs, x1, cls, y1, idx, NRidges, S)
sampleLocalN = numel(sampleIdx);
rowsOut = {};
rmseVals = [];

% Determine if DBSCAN separated curves well enough.
labs = unique(cls(:));
labs(labs <= 0) = [];

clusterTracks = {};
if numel(labs) >= NRidges
    % Use the largest NRidges clusters as candidate tracks.
    counts = arrayfun(@(z)sum(cls(:)==z), labs);
    [~,ord] = sort(counts, 'descend');
    labsUse = labs(ord(1:NRidges));

    for k = 1:numel(labsUse)
        lab = labsUse(k);
        mask = cls == lab;
        xFit = x1(mask);
        yMat = repmat(y1,1,size(x1,2));
        yFit = yMat(mask);
        good = isfinite(xFit) & isfinite(yFit);
        xFit = xFit(good);
        yFit = yFit(good);
        if numel(xFit) >= 2
            clusterTracks{end+1} = struct('xFit',xFit,'yFit',yFit,'label',lab); %#ok<AGROW>
        end
    end
end

% If DBSCAN merged close curves or did not form enough tracks, fall back to
% ordered peak columns. This mirrors the practical multiplet behavior:
% expected peaks/ridges define the number of tracks, and each track can only
% take one peak per sample.
if numel(clusterTracks) < NRidges
    clusterTracks = {};
    for k = 1:min(NRidges,size(x1,2))
        xFit = x1(:,k);
        yFit = y1;
        good = isfinite(xFit) & isfinite(yFit);
        clusterTracks{end+1} = struct('xFit',xFit(good),'yFit',yFit(good),'label',k); %#ok<AGROW>
    end
end

for k = 1:min(NRidges,numel(clusterTracks))
    tr = clusterTracks{k};
    if numel(tr.xFit) < 2
        continue;
    end

    [fitFcn, fitRMSE] = localPolyFitFunction(tr.yFit, tr.xFit, S.PolyOrder);
    rmseVals(end+1) = fitRMSE; %#ok<AGROW>

    [predCenter, predLow, predHigh] = localPredInterval(fitFcn, y1, tr.yFit, tr.xFit, S.PredictionLevel);

    rowsK = repickOneTrackWithinIntervals(ppm, X_reo, sampleIdx, predCenter, predLow, predHigh, ...
        min(ppms), max(ppms), k, tr.label);

    rowsOut = [rowsOut; rowsK]; %#ok<AGROW>
end

trace = cell2table_safe(rowsOut, {'SampleOrderId','SampleOrderIdx','PeakNumber','PeakPPM_Recon','PPMIndex', ...
    'ReconIntensity','FitPPM','ResidualPPM','ClusterLabel'});

if ~isempty(trace)
    trace = enforceOnePeakPerSamplePerTrack(trace);
    trace = sortrows(trace, {'PeakNumber','SampleOrderId'});
end

if isempty(rmseVals)
    rmse = NaN;
else
    rmse = median(rmseVals, 'omitnan');
end
end

%% =====================================================================
function rowsOut = repickOneTrackWithinIntervals(ppm, X_reo, sampleIdx, predCenter, predLow, predHigh, ppmMin, ppmMax, peakNum, clab)
rowsOut = {};
for s = 1:numel(sampleIdx)
    lo = min(predLow(s), predHigh(s));
    hi = max(predLow(s), predHigh(s));

    % Keep repicking constrained to original selected user ppm window.
    lo = max(lo, ppmMin);
    hi = min(hi, ppmMax);

    ii1 = localMatchPPM(lo, ppm);
    ii2 = localMatchPPM(hi, ppm);
    iiLo = max(1, min(ii1, ii2));
    iiHi = min(numel(ppm), max(ii1, ii2));

    y = X_reo(sampleIdx(s), iiLo:iiHi);
    if isempty(y)
        continue;
    end

    [pks, loc] = findpeaks(y, 'SortStr','descend');
    if isempty(loc)
        [~, loc] = max(y);
    else
        loc = loc(1);
    end

    globalIdx = iiLo + loc - 1;
    peakPPM = ppm(globalIdx);
    fitPPM = predCenter(s);

    rowsOut(end+1,:) = {sampleIdx(s), sampleIdx(s), peakNum, peakPPM, globalIdx, ...
        X_reo(sampleIdx(s), globalIdx), fitPPM, peakPPM-fitPPM, clab}; %#ok<AGROW>
end
end

%% =====================================================================
function trace = enforceOnePeakPerSamplePerTrack(trace)
if isempty(trace) || height(trace)==0
    return;
end

keep = false(height(trace),1);
tracks = unique(trace.PeakNumber);
for t = tracks(:)'
    rowsT = find(trace.PeakNumber == t);
    samples = unique(trace.SampleOrderId(rowsT));
    for s = samples(:)'
        rows = rowsT(trace.SampleOrderId(rowsT) == s);
        if numel(rows) == 1
            keep(rows) = true;
        else
            [~,bestLocal] = min(abs(trace.ResidualPPM(rows)));
            keep(rows(bestLocal)) = true;
        end
    end
end
trace = trace(keep,:);
end

%% =====================================================================
function [locV,pks1] = localFindPeaksSicong(Xin, modeName, NumPeak, x4)
% Local implementation of the behavior needed from ReorderAlign_FindPeaks.
% Returns exactly NumPeak peak indices per sample when possible.
if isnumeric(Xin)
    sample = size(Xin,1);
    getRow = @(i) Xin(i,:);
elseif iscell(Xin)
    sample = numel(Xin);
    getRow = @(i) Xin{i};
else
    error('Xin must be numeric matrix or cell array.');
end

locCell = cell(sample,1);
pksCell = cell(sample,1);

for ind = 1:sample
    y = getRow(ind);

    if isempty(y)
        loc = ones(1,NumPeak);
        pks = nan(1,NumPeak);
    else
        [pks, loc] = findpeaks(y, 'SortStr','descend');

        if isempty(loc)
            [~, loc] = max(y);
            pks = y(loc);
        end

        % Keep strongest NumPeak, then sort by location so multiplet/parallel
        % peaks become ordered columns.
        if numel(loc) < NumPeak
            for k = (numel(loc)+1):NumPeak
                if strcmpi(modeName, 'multi') && ~isempty(x4) && numel(x4) >= sample
                    if (x4(end)-x4(1))*(sample/2-ind) > 0
                        loc(k) = 1; %#ok<AGROW>
                    else
                        loc(k) = numel(y); %#ok<AGROW>
                    end
                else
                    loc(k) = round(numel(y)/2); %#ok<AGROW>
                end
            end
        end

        loc = loc(1:NumPeak);
        loc = sort(loc(:))';
        pks = y(loc);
    end

    locCell{ind} = loc(:)';
    pksCell{ind} = pks(:)';
end

locV = cell2mat(locCell);
pks1 = cell2mat(pksCell);
end

%% =====================================================================
function cls = localSicongCluster(ClsMat, x1, y1, sample, FigDisp, minPts)
% Same spirit/formula as ReorderAlign_Cluster.
try
    ClsScaled = scale(ClsMat,'auto');
catch
    mu = mean(ClsMat,1,'omitnan');
    sig = std(ClsMat,0,1,'omitnan');
    sig(sig==0 | ~isfinite(sig)) = 1;
    ClsScaled = (ClsMat - mu)./sig;
end

d = size(ClsScaled,2);
rangeVals = max(ClsScaled)-min(ClsScaled);
rangeVals(rangeVals <= 0 | ~isfinite(rangeVals)) = eps;
epdist = exp(sum(log(rangeVals))/d) * ...
    (3*sqrt(pi*d)/size(ClsMat,1))^(1/d) * ...
    sqrt(d/(2*exp(1)*pi));

try
    cls = dbscan(ClsScaled, epdist, minPts);
catch ME
    warning('DBSCAN failed: %s. Treating all candidates as one cluster.', ME.message);
    cls = ones(size(ClsScaled,1),1);
end

if size(ClsMat,1) > length(y1)
    cls = reshape(cls, size(x1));
end

if strcmp(FigDisp,'on')
    ucls = unique(cls);
    cmap = jet(length(ucls));
    figure; hold on
    if size(ClsMat,1) > length(y1)
        for ind = 1:length(ucls)
            for i = 1:size(x1,2)
                thisMask = cls(1+sample*(i-1):sample*i) == ucls(ind);
                scatter(x1(thisMask,i), y1(thisMask), 100, ...
                    'MarkerFaceColor', cmap(ind,:), 'MarkerEdgeColor','none');
            end
        end
    else
        for ind = 1:length(ucls)
            scatter(x1(cls==ucls(ind),1), y1(cls==ucls(ind)), 100, ...
                'MarkerFaceColor', cmap(ind,:), 'MarkerEdgeColor','none');
        end
    end
    hold off
    set(gca,'XDir','rev')
    title('DBSCAN-clustering result');
end
end

%% =====================================================================
function [fitFcn, rmse] = localPolyFitFunction(x, y, polyOrder)
x = x(:); y = y(:);
good = isfinite(x) & isfinite(y);
x = x(good); y = y(good);
po = min(polyOrder, max(numel(unique(x))-1, 1));
if numel(unique(x)) <= po
    po = 1;
end
try
    coeff = polyfit(x, y, po);
    fitFcn = @(t) polyval(coeff, t);
    pred = fitFcn(x);
    rmse = sqrt(mean((y-pred).^2, 'omitnan'));
catch
    coeff = polyfit(x, y, 1);
    fitFcn = @(t) polyval(coeff, t);
    pred = fitFcn(x);
    rmse = sqrt(mean((y-pred).^2, 'omitnan'));
end
end

%% =====================================================================
function [center, low, high] = localPredInterval(fitFcn, yAll, yFit, xFit, level)
center = fitFcn(yAll);

% Try Curve Fitting Toolbox predint for a closer Sicong match.
try
    ft = fit(yFit(:), xFit(:), 'poly3');
    intv = predint(ft, yAll(:), level, 'observation', 'on');
    low = intv(:,1);
    high = intv(:,2);
    center = ft(yAll(:));
catch
    predFit = fitFcn(yFit);
    res = xFit(:) - predFit(:);
    sigma = 1.4826 * median(abs(res - median(res,'omitnan')), 'omitnan');
    if ~isfinite(sigma) || sigma == 0
        sigma = sqrt(mean(res.^2, 'omitnan'));
    end
    if ~isfinite(sigma) || sigma == 0
        sigma = 0.003;
    end
    k = 2;
    if level >= 0.99
        k = 2.8;
    elseif level <= 0.90
        k = 1.65;
    end
    low = center(:) - k*sigma;
    high = center(:) + k*sigma;
end
end

%% =====================================================================
function mapped = mapTraceToSand(trace, order, sandCache, tolPPM)
rows = {};
for i = 1:height(trace)
    sOrd = trace.SampleOrderId(i);
    origIdx = order(sOrd);
    reconPPM = trace.PeakPPM_Recon(i);
    peakNum = trace.PeakNumber(i);

    cache = sandCache{origIdx};
    if isempty(cache.ppm)
        rows(end+1,:) = {sOrd, sOrd, origIdx, peakNum, NaN, NaN, reconPPM, NaN, NaN, NaN, false}; %#ok<AGROW>
        continue;
    end

    [d,j] = min(abs(cache.ppm - reconPPM));

    if isempty(d) || d > tolPPM
        rows(end+1,:) = {sOrd, sOrd, origIdx, peakNum, NaN, NaN, reconPPM, d, NaN, NaN, false}; %#ok<AGROW>
    else
        rows(end+1,:) = {sOrd, sOrd, origIdx, peakNum, cache.row(j), cache.ppm(j), reconPPM, d, cache.amp(j), cache.decay(j), true}; %#ok<AGROW>
    end
end

mapped = cell2table_safe(rows, {'SampleOrderId','SampleOrderIdx','OriginalSampleId','PeakNumber','PeakRow','PeakPPM', ...
    'ReconPPM','MapDistancePPM','Amplitude','DecayHz','Mapped'});
end

%% =====================================================================
function mapped = makeUnmappedTable(trace, order)
rows = {};
for i = 1:height(trace)
    sOrd = trace.SampleOrderId(i);
    if sOrd >= 1 && sOrd <= numel(order)
        origIdx = order(sOrd);
    else
        origIdx = NaN;
    end
    rows(end+1,:) = {sOrd, sOrd, origIdx, trace.PeakNumber(i), NaN, NaN, trace.PeakPPM_Recon(i), NaN, NaN, NaN, false}; %#ok<AGROW>
end

mapped = cell2table_safe(rows, {'SampleOrderId','SampleOrderIdx','OriginalSampleId','PeakNumber','PeakRow','PeakPPM', ...
    'ReconPPM','MapDistancePPM','Amplitude','DecayHz','Mapped'});
end

%% =====================================================================
function plotBaseViews(axTop, axBot, ppm, X_reo, order, sandCache, InterFactor)
[nSamples, ~] = size(X_reo);

hold(axTop,'on');
for i = nSamples:-1:1
    plot(axTop, ppm, X_reo(i,:) + (i-1)*InterFactor, 'k', 'LineWidth', 0.55);
end
hold(axTop,'off');
set(axTop, 'XDir','reverse');
ylabel(axTop, 'Intensity (offset)');
grid(axTop,'on');

allPPM = [];
allSamp = [];
for sOrd = 1:nSamples
    origIdx = order(sOrd);
    if origIdx <= numel(sandCache) && ~isempty(sandCache{origIdx}.ppm)
        allPPM = [allPPM; sandCache{origIdx}.ppm(:)]; %#ok<AGROW>
        allSamp = [allSamp; sOrd*ones(numel(sandCache{origIdx}.ppm),1)]; %#ok<AGROW>
    end
end

scatter(axBot, allPPM, allSamp, 12, [0.65 0.65 0.66], 'filled', ...
    'MarkerFaceAlpha', 0.75, 'MarkerEdgeAlpha', 0.75);
set(axBot, 'XDir','reverse', 'YDir','normal', 'YLim',[0.5 nSamples+0.5]);
xlabel(axBot, 'ppm');
ylabel(axBot, 'Sample index (reordered)');
title(axBot, 'Peak-table dots in reordered sample space');
grid(axBot,'on');
linkaxes([axTop axBot], 'x');
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
    error('SandTable must be table or struct.');
end

if ~iscell(peakCells)
    error('PeakTableField must contain a cell array of per-sample tables.');
end
end

%% =====================================================================
function sandCache = cacheSandPeaks(peakCells, ppmCol, ampCol, decayCol)
n = numel(peakCells);
sandCache = cell(n,1);
for i = 1:n
    Ti = peakCells{i};
    cache = struct('ppm', [], 'amp', [], 'decay', [], 'row', []);
    if isempty(Ti) || ~istable(Ti) || ~ismember(ppmCol, Ti.Properties.VariableNames)
        sandCache{i} = cache;
        continue;
    end

    cache.ppm = Ti.(ppmCol)(:);
    cache.row = (1:height(Ti))';

    if ismember(ampCol, Ti.Properties.VariableNames)
        cache.amp = Ti.(ampCol)(:);
    else
        cache.amp = nan(height(Ti),1);
    end

    if ismember(decayCol, Ti.Properties.VariableNames)
        cache.decay = Ti.(decayCol)(:);
    else
        cache.decay = nan(height(Ti),1);
    end

    sandCache{i} = cache;
end
end

%% =====================================================================
function T = cell2table_safe(C, names)
if isempty(C)
    T = cell2table(cell(0,numel(names)), 'VariableNames', names);
    return;
end
T = cell2table(C, 'VariableNames', names);
end

%% =====================================================================
function idx = localMatchPPM(val, ppm)
[~, idx] = min(abs(ppm(:) - val));
end

%% =====================================================================
function setReadyAndResume(fig)
if isvalid(fig)
    setappdata(fig, 'EndpointTracerAction', 'ready');
    uiresume(fig);
end
end

function setDoneAndResume(fig)
if isvalid(fig)
    setappdata(fig, 'EndpointTracerAction', 'done');
    uiresume(fig);
end
end

function endpointKeyPress(fig, event)
if ~isvalid(fig), return; end
switch lower(event.Key)
    case {'space','return','enter'}
        setappdata(fig, 'EndpointTracerAction', 'ready');
        uiresume(fig);
    case {'escape','q'}
        setappdata(fig, 'EndpointTracerAction', 'done');
        uiresume(fig);
end
end

function c = ridgeColor(k)
palette = lines(12);
c = palette(mod(k-1,size(palette,1))+1,:);
end

function out = ternary(cond,a,b)
if cond, out = a; else, out = b; end
end

function deleteIfValid(h)
if isempty(h), return; end
try
    h = h(isvalid(h));
    if ~isempty(h), delete(h); end
catch
end
end
