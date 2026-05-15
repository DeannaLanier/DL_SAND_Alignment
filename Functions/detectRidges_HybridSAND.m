function out = detectRidges_HybridSAND(ppm, reorderOut, SandTable, varargin)
% {Deanna Lanier 4.27.2026
% 
% detectRidgesGlobalPlusEndpoints_HybridSAND
%
% Combined semi-automatic ridge workflow:
%   1) Run detectRidgesGlobal first.
%   2) Let the user review/zoom the global detected ridges.
%   3) Optionally add more ridges interactively using the
%      endpoint/window selector:
%           traceRidgesByEndpoints
%   4) Merge global + manual endpoint/window ridges into one output struct
%      that can go directly into:
%           outVal = validateRidgesInteractive(out);
%           outAligned = alignAndReconstructValidatedRidges(outVal,...)
%
% IMPORTANT:
%   This wrapper uses the endpoint function, not the older
%   literal endpoint-path tracer. In this, the two clicks define the
%   ppm range and sample span; the user is prompted for the expected number
%   of peaks/ridges in that selected window.
%
% BASIC USAGE
%   outRidges = detectRidgesGlobalPlusEndpoints_HybridSAND( ...
%       ppm, pHIT.out, SandTable, ...
%       'PeakTableField', 'Con_AF', ...
%       'ppmCol', 'freq_ppm', ...
%       'ppmTol', 0.02, ...
%       'MaxSampleGap', 1, ...
%       'MinSamplesPerRidge', 22, ...
%       'MinPointsPerRidge', 22, ...
%       'CurveType', 'poly3', ...
%       'InterFactor', 100);
%
% OPTIONAL:
%   'EndpointTracerFunction' : name of endpoint/window function.
%                              default 'traceRidgesByEndpoints'
%   'EndpointExtraArgs'      : cell array of extra args passed to the endpoint
%                              function, e.g.
%                              {'PromptForNumPeaks',true,'PredictionLevel',0.95}

    p = inputParser;
    p.addRequired('ppm', @(v)isnumeric(v)&&isvector(v));
    p.addRequired('reorderOut', @isstruct);
    p.addRequired('SandTable', @(x)istable(x)||isstruct(x));

    % Shared/global options
    p.addParameter('PeakTableField', 'Misaligned_Table', @(s)ischar(s)||isstring(s));
    p.addParameter('ppmCol', 'Aligned_PPM', @(s)ischar(s)||isstring(s));
    p.addParameter('ppmTol', 0.01, @(v)isnumeric(v)&&isscalar(v)&&v>0);
    p.addParameter('MaxSampleGap', 1, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
    p.addParameter('MinSamplesPerRidge', 4, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
    p.addParameter('MinPointsPerRidge', 4, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
    p.addParameter('CurveType', 'poly3', @(s)ischar(s)||isstring(s));
    p.addParameter('InterFactor', [], @(v)isempty(v)||(isnumeric(v)&&isscalar(v)));
    p.addParameter('MakePlots', true, @(b)islogical(b)&&isscalar(b));

    % Endpoint options
    p.addParameter('AddEndpointRidges', true, @(b)islogical(b)&&isscalar(b));
    p.addParameter('AskBeforeEndpoint', true, @(b)islogical(b)&&isscalar(b));
    p.addParameter('ReviewGlobalBeforeEndpoint', true, @(b)islogical(b)&&isscalar(b));
    p.addParameter('MapTolerancePPM', 0.0025, @(v)isnumeric(v)&&isscalar(v));
    p.addParameter('EndpointTracerFunction', 'traceRidgesByEndpoints', @(s)ischar(s)||isstring(s));

    % Pass-through extras
    p.addParameter('GlobalExtraArgs', {}, @(c)iscell(c));
    p.addParameter('EndpointExtraArgs', {}, @(c)iscell(c));

    p.parse(ppm, reorderOut, SandTable, varargin{:});
    S = p.Results;

    peakField = char(S.PeakTableField);
    ppmCol    = char(S.ppmCol);
    endpointFcnName = char(S.EndpointTracerFunction);

    % ---------------- 1) Run global detector ----------------
    globalArgs = { ...
        'PeakTableField', peakField, ...
        'ppmCol', ppmCol, ...
        'ppmTol', S.ppmTol, ...
        'MaxSampleGap', S.MaxSampleGap, ...
        'MinSamplesPerRidge', S.MinSamplesPerRidge, ...
        'MinPointsPerRidge', S.MinPointsPerRidge, ...
        'CurveType', char(S.CurveType), ...
        'MakePlots', S.MakePlots};

    if ~isempty(S.InterFactor)
        globalArgs = [globalArgs, {'InterFactor', S.InterFactor}]; %#ok<AGROW>
    end
    globalArgs = [globalArgs, S.GlobalExtraArgs];

    fprintf('\nRunning global ridge detection...\n');
    globalOut = detectRidgesGlobal(ppm, reorderOut, SandTable, globalArgs{:});

    if isfield(globalOut,'ridgeSummaryTable') && ~isempty(globalOut.ridgeSummaryTable)
        nGlobal = height(globalOut.ridgeSummaryTable);
    else
        nGlobal = 0;
    end
    fprintf('Global detection returned %d ridge(s).\n', nGlobal);

    % Let the user inspect/zoom the detected global ridges before deciding
    % whether to add endpoint/window ridges.
    if S.MakePlots && S.ReviewGlobalBeforeEndpoint && isfield(globalOut,'fig') && isgraphics(globalOut.fig)
        figure(globalOut.fig);
        fprintf('\nGlobal ridges are displayed. Zoom/pan the figure now.\n');
        fprintf('When you are done reviewing, click the Command Window and press Enter.\n');
        input('Press Enter after reviewing global ridges: ', 's');
    end

    endpointOut = [];
    runEndpoint = S.AddEndpointRidges;

    % ---------------- 2) Ask whether to add endpoint/window ridges ----------------
    if runEndpoint && S.AskBeforeEndpoint
        if S.MakePlots && isfield(globalOut,'fig') && isgraphics(globalOut.fig)
            figure(globalOut.fig);
        end
        choice = questdlg( ...
            sprintf('Global detection found %d ridge(s). Add ridges?', nGlobal), ...
            'Add endpoint/window ridges?', ...
            'Yes','No','Yes');
        runEndpoint = strcmp(choice, 'Yes');
    end

    if runEndpoint
        fprintf('\nStarting ridge tracing.\n');
        fprintf('Use it to add ridges missed by global detection.\n');

        if exist(endpointFcnName,'file') ~= 2
            error(['Endpoint tracer function "%s" was not found on the MATLAB path. ', ...
                   'Make sure traceRidgesByEndpoints.m is on the path.'], endpointFcnName);
        end

        endpointArgs = { ...
            'PeakTableField', peakField, ...
            'ppmCol', ppmCol, ...
            'MapTolerancePPM', S.MapTolerancePPM, ...
            'MakeFigure', true};

        if ~isempty(S.InterFactor)
            endpointArgs = [endpointArgs, {'InterFactor', S.InterFactor}]; %#ok<AGROW>
        end

        endpointArgs = [endpointArgs, S.EndpointExtraArgs];

        endpointFcn = str2func(endpointFcnName);
        endpointOut = endpointFcn(ppm, reorderOut, SandTable, endpointArgs{:});

        if isstruct(endpointOut) && isfield(endpointOut,'ridgeSummaryTable') && ~isempty(endpointOut.ridgeSummaryTable)
            fprintf('Endpoint/window tracing returned %d accepted ridge(s).\n', height(endpointOut.ridgeSummaryTable));
        else
            fprintf('Endpoint/window tracing returned 0 accepted ridge(s).\n');
        end
    end

    % ---------------- 3) Merge outputs ----------------
    if isempty(endpointOut) || ~isstruct(endpointOut) || ...
            ~isfield(endpointOut,'ridgeSummaryTable') || isempty(endpointOut.ridgeSummaryTable)
        out = globalOut;
        out.globalOut = globalOut;
        out.endpointOut = endpointOut;
        out.SourceWorkflow = 'global_only';
        return;
    end

    out = mergeGlobalAndEndpointOutputs(globalOut, endpointOut, reorderOut);
    out.globalOut = globalOut;
    out.endpointOut = endpointOut;
    out.SourceWorkflow = 'global_plus_endpoint_window';

    % ---------------- 4) Overlay endpoint ridges on global plot ----------------
    if isfield(out,'axBot') && isgraphics(out.axBot)
        if nGlobal > 0
            offset = max(globalOut.ridgeSummaryTable.RidgeID);
        else
            offset = 0;
        end
        overlayEndpointRidges(out.axBot, endpointOut, offset);
    end

    fprintf('\nCombined output contains %d total ridge(s): %d global + %d endpoint/window.\n', ...
        height(out.ridgeSummaryTable), nGlobal, height(endpointOut.ridgeSummaryTable));
end

% =====================================================================
function out = mergeGlobalAndEndpointOutputs(globalOut, endpointOut, reorderOut)
    gSum = globalOut.ridgeSummaryTable;
    gPts = globalOut.ridgePointsTable;
    eSum = endpointOut.ridgeSummaryTable;
    ePts = endpointOut.ridgePointsTable;

    if isempty(gSum)
        offset = 0;
    else
        offset = max(gSum.RidgeID);
    end

    % Renumber endpoint ridges.
    oldIDs = eSum.RidgeID;
    newIDs = oldIDs + offset;
    for i = 1:numel(oldIDs)
        ePts.RidgeID(ePts.RidgeID == oldIDs(i)) = newIDs(i);
    end
    eSum.RidgeID = newIDs;

    % Harmonize point tables.
    gPts = ensurePointColumns(gPts, reorderOut, 'global');
    ePts = ensurePointColumns(ePts, reorderOut, 'endpoint_window');
    allPtVars = unique([gPts.Properties.VariableNames, ePts.Properties.VariableNames], 'stable');
    gPts = addMissingVars(gPts, allPtVars);
    ePts = addMissingVars(ePts, allPtVars);
    gPts = gPts(:, allPtVars);
    ePts = ePts(:, allPtVars);
    mergedPts = [gPts; ePts];

    % Harmonize summary tables.
    gSum = ensureSummaryColumns(gSum, gPts, 'global');
    eSum = ensureSummaryColumns(eSum, ePts, 'endpoint_window');
    allSumVars = unique([gSum.Properties.VariableNames, eSum.Properties.VariableNames], 'stable');
    gSum = addMissingVars(gSum, allSumVars);
    eSum = addMissingVars(eSum, allSumVars);
    gSum = gSum(:, allSumVars);
    eSum = eSum(:, allSumVars);
    mergedSum = [gSum; eSum];

    out = globalOut;
    out.ridgeSummaryTable = mergedSum;
    out.ridgePointsTable = mergedPts;
    out.order = reorderOut.order(:);
    if isfield(reorderOut,'X_reordered')
        out.X_reordered = reorderOut.X_reordered;
    end

    % Merge ridgeCurves if possible. Endpoint curves are represented simply
    % using reconRidgePointsTable when available.
    if isfield(globalOut,'ridgeCurves')
        out.ridgeCurves = globalOut.ridgeCurves;
    else
        out.ridgeCurves = struct([]);
    end
    if isfield(endpointOut,'reconRidgePointsTable') && ~isempty(endpointOut.reconRidgePointsTable)
        rcNew = endpointReconToCurveStruct(endpointOut.reconRidgePointsTable, offset);
        out.ridgeCurves = [out.ridgeCurves, rcNew];
    end
end

% =====================================================================
function T = ensurePointColumns(T, reorderOut, sourceName)
    if isempty(T)
        return;
    end

    if ~ismember('OriginalSampleId', T.Properties.VariableNames)
        order = reorderOut.order(:);
        orig = nan(height(T),1);
        if ismember('SampleOrderId', T.Properties.VariableNames)
            s = T.SampleOrderId;
            good = s >= 1 & s <= numel(order);
            orig(good) = order(s(good));
        end
        T.OriginalSampleId = orig;
    end

    if ~ismember('PeakNumber', T.Properties.VariableNames)
        T.PeakNumber = ones(height(T),1);
    end

    if ~ismember('Residual', T.Properties.VariableNames)
        if ismember('MapDistancePPM', T.Properties.VariableNames)
            T.Residual = T.MapDistancePPM;
        else
            T.Residual = nan(height(T),1);
        end
    end

    if ~ismember('Inlier', T.Properties.VariableNames)
        T.Inlier = true(height(T),1);
    end

    if ~ismember('Source', T.Properties.VariableNames)
        T.Source = repmat(string(sourceName), height(T), 1);
    end
end

% =====================================================================
function T = ensureSummaryColumns(T, ptTbl, sourceName)
    if isempty(T)
        return;
    end

    if ~ismember('Source', T.Properties.VariableNames)
        T.Source = repmat(string(sourceName), height(T), 1);
    end

    if ~ismember('NumPoints', T.Properties.VariableNames)
        T.NumPoints = nan(height(T),1);
        for i = 1:height(T)
            T.NumPoints(i) = sum(ptTbl.RidgeID == T.RidgeID(i));
        end
    end

    if ~ismember('NumSamples', T.Properties.VariableNames)
        T.NumSamples = nan(height(T),1);
        for i = 1:height(T)
            mask = ptTbl.RidgeID == T.RidgeID(i);
            if any(mask) && ismember('SampleOrderId', ptTbl.Properties.VariableNames)
                T.NumSamples(i) = numel(unique(ptTbl.SampleOrderId(mask)));
            end
        end
    end

    if ~ismember('MinPPM', T.Properties.VariableNames)
        if ismember('RegionMinPPM', T.Properties.VariableNames)
            T.MinPPM = T.RegionMinPPM;
        else
            T.MinPPM = nan(height(T),1);
            for i = 1:height(T)
                mask = ptTbl.RidgeID == T.RidgeID(i);
                if any(mask), T.MinPPM(i) = min(ptTbl.PeakPPM(mask)); end
            end
        end
    end

    if ~ismember('MaxPPM', T.Properties.VariableNames)
        if ismember('RegionMaxPPM', T.Properties.VariableNames)
            T.MaxPPM = T.RegionMaxPPM;
        else
            T.MaxPPM = nan(height(T),1);
            for i = 1:height(T)
                mask = ptTbl.RidgeID == T.RidgeID(i);
                if any(mask), T.MaxPPM(i) = max(ptTbl.PeakPPM(mask)); end
            end
        end
    end

    if ~ismember('CoverageFrac', T.Properties.VariableNames)
        if ismember('MappedCoverage', T.Properties.VariableNames)
            T.CoverageFrac = T.MappedCoverage;
        else
            T.CoverageFrac = nan(height(T),1);
        end
    end

    if ~ismember('CurveType', T.Properties.VariableNames)
        T.CurveType = repmat("endpoint", height(T), 1);
    end

    if ~ismember('RMSE_ppm', T.Properties.VariableNames)
        if ismember('FitRMSE', T.Properties.VariableNames)
            T.RMSE_ppm = T.FitRMSE;
        else
            T.RMSE_ppm = nan(height(T),1);
        end
    end

    if ~ismember('Direction', T.Properties.VariableNames)
        T.Direction = repmat("n/a", height(T), 1);
    end
end

% =====================================================================
function T = addMissingVars(T, varNames)
    for i = 1:numel(varNames)
        v = varNames{i};
        if ~ismember(v, T.Properties.VariableNames)
            T.(v) = makeDefaultColumn(height(T), v);
        end
    end
end

% =====================================================================
function col = makeDefaultColumn(n, varName)
    stringVars = {'Direction','CurveType','Source'};
    logicalVars = {'Inlier','Accepted','Mapped'};
    if any(strcmp(varName, stringVars))
        col = strings(n,1);
    elseif any(strcmp(varName, logicalVars))
        col = false(n,1);
    else
        col = nan(n,1);
    end
end

% =====================================================================
function rc = endpointReconToCurveStruct(reconTbl, offset)
    if isempty(reconTbl)
        rc = struct([]);
        return;
    end
    ids = unique(reconTbl.RidgeID);
    rc = struct([]);
    for i = 1:numel(ids)
        oldID = ids(i);
        newID = oldID + offset;
        rows = reconTbl.RidgeID == oldID;
        x = reconTbl.SampleOrderId(rows);
        if ismember('PeakPPM_Recon', reconTbl.Properties.VariableNames)
            y = reconTbl.PeakPPM_Recon(rows);
        else
            y = nan(sum(rows),1);
        end
        [x, id] = sort(x);
        y = y(id);
        rc(end+1).RidgeID = newID; %#ok<AGROW>
        rc(end).sampleId = x;
        rc(end).ppm = y;
        rc(end).fitType = 'endpoint_window';
        rc(end).fitFcn = [];
        rc(end).xGrid = x;
        rc(end).yHat = y;
        rc(end).lowB = y;
        rc(end).hiB = y;
        rc(end).residuals = nan(size(y));
        rc(end).sigmaResidual = NaN;
        rc(end).bandWidth = NaN;
        rc(end).inlierMask = true(size(y));
        rc(end).outlierMask = false(size(y));
    end
end

% =====================================================================
function overlayEndpointRidges(axBot, endpointOut, offset)
    if isempty(endpointOut) || ~isfield(endpointOut,'ridgePointsTable') || isempty(endpointOut.ridgePointsTable)
        return;
    end
    ePts = endpointOut.ridgePointsTable;
    ids = unique(ePts.RidgeID);
    cmap = lines(max(1,numel(ids)));
    hold(axBot,'on');
    for i = 1:numel(ids)
        oldID = ids(i);
        newID = oldID + offset;
        rows = ePts.RidgeID == oldID;
        c = cmap(i,:);
        [~, ord] = sort(ePts.SampleOrderId(rows));
        ppmVals = ePts.PeakPPM(rows);
        sampVals = ePts.SampleOrderId(rows);
        ppmVals = ppmVals(ord);
        sampVals = sampVals(ord);
        plot(axBot, ppmVals, sampVals, 'o-', ...
            'Color', c, 'MarkerFaceColor', c, 'MarkerEdgeColor','k', ...
            'LineWidth', 1.6, 'MarkerSize', 5, ...
            'DisplayName', sprintf('endpoint ridge %d', newID));
    end
    hold(axBot,'off');
    drawnow;
end
