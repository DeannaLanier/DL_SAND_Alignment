function out = detectRidgesGlobal(ppm, reorderOut, alignedData, varargin)
% detectRidgesGlobal
%
% Globally detect ridges/curves across the reordered sample-by-ppm space
% by tracking peaks across samples, then:
%   - cleaning each ridge with a DBSCAN-based filter,
%   - fitting a user-chosen model:
%         'linear' : straight line
%         'poly3'  : cubic polynomial
%         'exp'    : exponential (exp1),
%   - computing a residual-based band and flagging off-curve outliers.
%
% This mirrors the logic used in:
%   - detectGuideCurvesDBSCAN
%   - detectRespondingCurvesFromGuide
%
% INPUTS
%   ppm         : 1 x N chemical shift vector
%   reorderOut  : struct with at least:
%                   - X_reordered  (nSamples x N)
%                   - order        (nSamples x 1, reordered -> original index)
%   alignedData : table or struct containing peak tables; by default
%                 expects field/variable 'Misaligned_Table' with a cell
%                 array of tables (one per sample) with a ppm column.
%
% NAME–VALUE PAIRS
%   'PeakTableField'    : field/var name with per-sample peak tables
%                         (default: 'Misaligned_Table')
%   'ppmCol'            : ppm column in peak tables (default: 'Aligned_PPM')
%
%   ---- Ridge construction ----
%   'ppmTol'            : maximum |ppm_current - ppm_last| to attach peak
%                         to an existing ridge (default: 0.01 ppm)
%   'MaxSampleGap'      : maximum sample gap to continue a ridge
%                         (default: 1 → only connect to previous sample)
%   'MinSamplesPerRidge': minimum distinct samples in a ridge (default: 4)
%   'MinPointsPerRidge' : minimum total points in a ridge (default: 4)
%
%   ---- Shape & banding ----
%   'CurveType'         : 'linear' | 'poly3' | 'exp' (default: 'linear')
%   'MonoTol'           : tolerance for monotonicity (ppm vs sample) (default: 1e-5)
%   'alphaBand'         : central band level for in-curve band (default: 0.99)
%   'minCoverageFrac'   : minimum fraction of samples spanned by inliers
%                         (default: 0.3)
%   'rmseMax'           : maximum allowed RMSE on inliers (ppm, default: 0.004)
%
%   ---- Plotting / GUI integration ----
%   'FontName'          : font for plots (default: 'Times New Roman')
%   'TickSize'          : tick label font size (default: 12)
%   'LabelSize'         : axis label font size (default: 14)
%   'TitleSize'         : title font size (default: 16)
%   'InterFactor'       : vertical spacing between stacked spectra in top plot
%                         (default: 0.02*max|X|)
%   'MakePlots'         : logical, whether to draw plots (default: true)
%   'AxesTop'           : existing axes handle for top FR plot (optional)
%   'AxesBot'           : existing axes handle for bottom dots/ridges (optional)
%
% OUTPUT STRUCT
%   out.X_reordered       : from reorderOut
%   out.order             : from reorderOut
%   out.ppm               : ppm
%   out.ridgePointsTable  : table with membership:
%                           RidgeID, SampleOrderId, PeakPPM, PeakRow,
%                           Residual, Inlier (logical)
%   out.ridgeSummaryTable : summary per ridge:
%                           RidgeID, NumPoints, NumSamples, CoverageFrac,
%                           MinPPM, MaxPPM, Direction, CurveType,
%                           SigmaResidual, BandWidth, RMSE_ppm
%   out.ridgeCurves       : struct array with fields:
%                           RidgeID, sampleId, ppm, fitType, fitFcn,
%                           xGrid, yHat, lowB, hiB, residuals,
%                           sigmaResidual, bandWidth, inlierMask, outlierMask
%   out.fig               : figure handle ([] if MakePlots=false and no axes)
%   out.axTop, out.axBot  : axes (may be [] if not plotted)

% ---------- Parse inputs ----------
p = inputParser;
p.addParameter('PeakTableField',    'Misaligned_Table', @(s)ischar(s)||isstring(s));
p.addParameter('ppmCol',            'Aligned_PPM',      @(s)ischar(s)||isstring(s));

% Ridge construction
p.addParameter('ppmTol',            0.01,               @(v)isnumeric(v)&&isscalar(v)&&v>0);
p.addParameter('MaxSampleGap',      1,                  @(v)isnumeric(v)&&isscalar(v)&&v>=1);
p.addParameter('MinSamplesPerRidge',4,                  @(v)isnumeric(v)&&isscalar(v)&&v>=1);
p.addParameter('MinPointsPerRidge', 4,                  @(v)isnumeric(v)&&isscalar(v)&&v>=1);

% Shape & banding
p.addParameter('CurveType',         'linear',           @(s)ischar(s)||isstring(s));
p.addParameter('MonoTol',           1e-5,               @(v)isnumeric(v)&&isscalar(v)&&v>=0);
p.addParameter('alphaBand',         0.99,               @(v)isnumeric(v)&&isscalar(v)&&v>0&&v<1);
p.addParameter('minCoverageFrac',   0.3,                @(v)isnumeric(v)&&isscalar(v)&&v>0&&v<=1);
p.addParameter('rmseMax',           0.004,              @(v)isnumeric(v)&&isscalar(v)&&v>0);

% Plot / GUI options
p.addParameter('FontName',          'Times New Roman',  @(s)ischar(s)||isstring(s));
p.addParameter('TickSize',          12,                 @(v)isnumeric(v)&&isscalar(v));
p.addParameter('LabelSize',         14,                 @(v)isnumeric(v)&&isscalar(v));
p.addParameter('TitleSize',         16,                 @(v)isnumeric(v)&&isscalar(v));
p.addParameter('InterFactor',       [],                 @(v)isnumeric(v)&&isscalar(v));
p.addParameter('MakePlots',         true,               @(v)islogical(v)&&isscalar(v));
p.addParameter('AxesTop',           [],                 @(h) isempty(h) || isgraphics(h,'axes'));
p.addParameter('AxesBot',           [],                 @(h) isempty(h) || isgraphics(h,'axes'));

p.parse(varargin{:});
S = p.Results;

ppm          = ppm(:)';  % ensure row
peakField    = char(S.PeakTableField);
ppmCol       = char(S.ppmCol);
ppmTol       = S.ppmTol;
maxGap       = S.MaxSampleGap;
minSampRidge = S.MinSamplesPerRidge;
minPtsRidge  = S.MinPointsPerRidge;
monoTol      = S.MonoTol;
curveType    = lower(string(S.CurveType));
alphaBand    = S.alphaBand;
minCovFrac   = S.minCoverageFrac;
rmseMax      = S.rmseMax;

fontName     = S.FontName;
tickSize     = S.TickSize;
labelSize    = S.LabelSize;
titleSize    = S.TitleSize;

% z-band for given alpha
if exist('norminv','file')
    zBand = norminv(0.5 + alphaBand/2);
else
    % approximate: alphaBand=0.99 → z ≈ 2.58
    zBand = 2.58;
end

% ---------- Unpack reorderOut ----------
if ~isfield(reorderOut, 'X_reordered') || ~isfield(reorderOut, 'order')
    error('reorderOut must contain fields X_reordered and order.');
end
X_reo = reorderOut.X_reordered;
order = reorderOut.order(:);

[nSamples, nPoints] = size(X_reo);
if numel(ppm) ~= nPoints
    error('ppm length (%d) must match columns of X_reordered (%d).', ...
        numel(ppm), nPoints);
end

% InterFactor (for top stacked plot, if drawn)
if isempty(S.InterFactor)
    InterFactor = 0.02 * max(abs(X_reo(:)));
else
    InterFactor = S.InterFactor;
end

% Precompute stacked spectra (even if not plotted)
X_stack = X_reo + repmat(InterFactor*(0:nSamples-1)',1,nPoints);

% ---------- Unpack alignedData & peak tables ----------
if istable(alignedData)
    if ~ismember(peakField, alignedData.Properties.VariableNames)
        error('alignedData table lacks variable "%s".', peakField);
    end
    Tcell = alignedData.(peakField);
elseif isstruct(alignedData)
    if ~isfield(alignedData, peakField)
        error('alignedData struct lacks field "%s".', peakField);
    end
    Tcell = alignedData.(peakField);
else
    error('alignedData must be a table or struct containing "%s".', peakField);
end

if numel(Tcell) < nSamples
    error('Peak table cell array (%d) has fewer entries than samples (%d).', ...
        numel(Tcell), nSamples);
end

% ---------- Collect per-sample peaks in reordered order ----------
perSamplePPM = cell(nSamples,1);
perSampleRow = cell(nSamples,1);

for sOrd = 1:nSamples
    origId = order(sOrd);
    Ti = Tcell{origId};
    if isempty(Ti) || ~ismember(ppmCol, Ti.Properties.VariableNames)
        perSamplePPM{sOrd} = [];
        perSampleRow{sOrd} = [];
        continue;
    end
    ppmv = Ti.(ppmCol);
    ppmv = ppmv(:);
    perSamplePPM{sOrd} = ppmv;
    perSampleRow{sOrd} = (1:numel(ppmv))';
end

% ---------- Build ridges by sequential linking ----------
ridges = struct('sampList',{}, 'ppmList',{}, 'rowList',{}, ...
                'lastSample',{}, 'lastPPM',{});

for sOrd = 1:nSamples
    ppmv = perSamplePPM{sOrd};
    rowv = perSampleRow{sOrd};
    if isempty(ppmv), continue; end

    for k = 1:numel(ppmv)
        pCurr   = ppmv(k);
        rowCurr = rowv(k);

        % Find candidate ridges whose last sample is within maxGap
        candId = [];
        for r = 1:numel(ridges)
            gap = sOrd - ridges(r).lastSample;
            if gap >= 1 && gap <= maxGap
                if abs(pCurr - ridges(r).lastPPM) <= ppmTol
                    % ensure this sample not already in the ridge
                    if ~ismember(sOrd, ridges(r).sampList)
                        candId(end+1) = r; %#ok<AGROW>
                    end
                end
            end
        end

        if isempty(candId)
            % start new ridge
            ridges(end+1).sampList   = sOrd;        %#ok<AGROW>
            ridges(end).ppmList      = pCurr;
            ridges(end).rowList      = rowCurr;
            ridges(end).lastSample   = sOrd;
            ridges(end).lastPPM      = pCurr;
        else
            % attach to the best matching ridge (closest ppm)
            bestR   = candId(1);
            bestDif = abs(pCurr - ridges(bestR).lastPPM);
            for r = candId(2:end)
                d = abs(pCurr - ridges(r).lastPPM);
                if d < bestDif
                    bestDif = d;
                    bestR   = r;
                end
            end

            ridges(bestR).sampList(end+1) = sOrd;
            ridges(bestR).ppmList(end+1)  = pCurr;
            ridges(bestR).rowList(end+1)  = rowCurr;
            ridges(bestR).lastSample      = sOrd;
            ridges(bestR).lastPPM         = pCurr;
        end
    end
end

% ---------- Convert ridges to tables & apply filters ----------
ridgeRows   = [];   % summary
ptRows      = [];   % membership
ridgeCurves = struct('RidgeID',{}, 'sampleId',{}, 'ppm',{}, ...
                     'fitType',{}, 'fitFcn',{}, ...
                     'xGrid',{}, 'yHat',{}, 'lowB',{}, 'hiB',{}, ...
                     'residuals',{}, 'sigmaResidual',{}, ...
                     'bandWidth',{}, 'inlierMask',{}, 'outlierMask',{});

for r = 1:numel(ridges)
    sampList = ridges(r).sampList(:);
    ppmList  = ridges(r).ppmList(:);
    rowList  = ridges(r).rowList(:);

    if isempty(sampList), continue; end

    % Sort by sample index
    [sampSorted, ord] = sort(sampList);
    ppmSorted = ppmList(ord);
    rowSorted = rowList(ord);

    % ---- DBSCAN-based noise removal along the ridge ----
    [keepMask, ~] = local_dbscan_clean(sampSorted, ppmSorted);
    if isempty(keepMask) || ~any(keepMask)
        continue;
    end

    sampClean = sampSorted(keepMask);
    ppmClean  = ppmSorted(keepMask);
    rowClean  = rowSorted(keepMask);

    nPtsRaw = numel(ppmClean);
    uSamp   = unique(sampClean);
    nSamp   = numel(uSamp);

    % Filter by minimum samples / points (after cleanup)
    if nSamp < minSampRidge || nPtsRaw < minPtsRidge
        continue;
    end

    coverageFrac = nSamp / nSamples;

    x = sampClean;
    y = ppmClean;

    % ---- Fit chosen curve type ----
    [yFitData, fitFcn] = local_fit_curve(x, y, curveType);
    if any(isnan(yFitData))
        continue;
    end

    res   = y - yFitData;
    % robust sigma (MAD-based)
    sigma = 1.4826 * mad(res, 1);
    if sigma == 0
        sigma = std(res);
    end
    bandW = zBand * sigma;

    inlierMask  = abs(res) <= bandW;
    outlierMask = ~inlierMask;

    if ~any(inlierMask)
        continue;
    end

    % Keep only in-curve points for defining the ridge
    xIn  = x(inlierMask);
    yIn  = y(inlierMask);
    rIn  = res(inlierMask);
    rowIn = rowClean(inlierMask);

    uSampIn   = unique(xIn);
    nSampIn   = numel(uSampIn);
    covInFrac = nSampIn / nSamples;
    rmse      = sqrt(mean(rIn.^2));

    if covInFrac < minCovFrac || rmse > rmseMax
        continue;
    end

    % ---- Monotonicity classification using inliers ----
    dy = diff(yIn);
    isMonoInc = all(dy >= -monoTol);
    isMonoDec = all(dy <=  monoTol);

    if isMonoInc && ~isMonoDec
        direction = "increasing";
    elseif isMonoDec && ~isMonoInc
        direction = "decreasing";
    elseif all(abs(dy) <= monoTol)
        direction = "flat";
    else
        direction = "nonmonotonic";
    end

    % ---- Evaluate curve on full sample range for plotting ----
    xGrid = (min(xIn):max(xIn))';
    yHat  = fitFcn(xGrid);
    lowB  = yHat - bandW;
    hiB   = yHat + bandW;

    ridgeID = r;

    % Summary row
    ridgeRows = [ridgeRows; table( ...
        ridgeID, nPtsRaw, nSampIn, covInFrac, ...
        min(yIn), max(yIn), direction, string(curveType), ...
        sigma, bandW, rmse, ...
        'VariableNames', { ...
            'RidgeID','NumPoints','NumSamples','CoverageFrac', ...
            'MinPPM','MaxPPM','Direction','CurveType', ...
            'SigmaResidual','BandWidth','RMSE_ppm'})]; %#ok<AGROW>

    % Membership rows (only inliers, but we keep residual & inlier flag)
    thisPts = table( ...
        repmat(ridgeID, numel(xIn),1), ...
        xIn, ...
        yIn, ...
        rowIn, ...
        rIn, ...
        true(numel(xIn),1), ...
        'VariableNames', {'RidgeID','SampleOrderId','PeakPPM', ...
                          'PeakRow','Residual','Inlier'});
    ptRows = [ptRows; thisPts]; %#ok<AGROW>

    % Store ridge curve
    ridgeCurves(end+1).RidgeID       = ridgeID;      %#ok<AGROW>
    ridgeCurves(end).sampleId       = xIn;
    ridgeCurves(end).ppm             = yIn;
    ridgeCurves(end).fitType         = curveType;
    ridgeCurves(end).fitFcn          = fitFcn;
    ridgeCurves(end).xGrid           = xGrid;
    ridgeCurves(end).yHat            = yHat;
    ridgeCurves(end).lowB            = lowB;
    ridgeCurves(end).hiB             = hiB;
    ridgeCurves(end).residuals       = rIn;
    ridgeCurves(end).sigmaResidual   = sigma;
    ridgeCurves(end).bandWidth       = bandW;
    ridgeCurves(end).inlierMask      = inlierMask;
    ridgeCurves(end).outlierMask     = outlierMask;
end

ridgeSummaryTable = ridgeRows;
ridgePointsTable  = ptRows;

% ---------- Build plotting axes (if requested) ----------
fig   = [];
axTop = [];
axBot = [];

if S.MakePlots
    % Collect all peaks for background dots
    allPPM = cell2mat(perSamplePPM(:));
    allSampId = cell2mat(arrayfun(@(sOrd) ...
        sOrd*ones(numel(perSamplePPM{sOrd}),1), ...
        (1:nSamples)', 'UniformOutput',false));

    % Case 1: user provided AxesTop/AxesBot → reuse them
    if ~isempty(S.AxesTop) || ~isempty(S.AxesBot)
        axTop = S.AxesTop;
        axBot = S.AxesBot;

        % Top axis
        if ~isempty(axTop)
            axes(axTop); %#ok<LAXES>
            cla(axTop);
            plot(axTop, ppm, X_stack', 'LineWidth', 1.0);
            set(axTop, 'XDir','reverse', 'FontName',fontName, 'FontSize',tickSize);
            ylabel(axTop, 'Intensity (offset)', 'FontName',fontName,'FontSize',labelSize);
            title(axTop, 'Stacked reordered spectra', 'FontName',fontName,'FontSize',titleSize);
            grid(axTop,'on');
        end

        % Bottom axis
        if ~isempty(axBot)
            axes(axBot); %#ok<LAXES>
            cla(axBot);

            scatter(axBot, allPPM, allSampId, 10, [0.45 0.45 0.45], 'filled', ...
                'MarkerFaceAlpha',0.6, 'MarkerEdgeAlpha',0.6);

            set(axBot, 'XDir','reverse', 'YDir','normal', ...
                'YLim',[0.5 nSamples+0.5], ...
                'FontName',fontName,'FontSize',tickSize);
            xlabel(axBot, sprintf('%s (ppm)', ppmCol), ...
                'FontName',fontName,'FontSize',labelSize);
            ylabel(axBot, 'Sample index (reordered)', ...
                'FontName',fontName,'FontSize',labelSize);
            title(axBot, sprintf('Global ridges (%s, banded)', curveType), ...
                'FontName',fontName,'FontSize',titleSize);
            grid(axBot,'on');
            hold(axBot,'on');

            cmap = lines(max(1, height(ridgeSummaryTable)));

            for i = 1:height(ridgeSummaryTable)
                rid = ridgeSummaryTable.RidgeID(i);
                col = cmap(i,:);
                rows = ridgePointsTable.RidgeID == rid;

                sampR = ridgePointsTable.SampleOrderId(rows);
                ppmR  = ridgePointsTable.PeakPPM(rows);

                % inliers (all points in ridgePointsTable are inliers by construction)
                scatter(axBot, ppmR, sampR, 36, col, 'filled', 'MarkerFaceAlpha',0.9);

                % plot curve and bands
                rcId = find([ridgeCurves.RidgeID] == rid, 1);
                if ~isempty(rcId)
                    xG   = ridgeCurves(rcId).xGrid;
                    yHat = ridgeCurves(rcId).yHat;
                    lowB = ridgeCurves(rcId).lowB;
                    hiB  = ridgeCurves(rcId).hiB;

                    plot(axBot, yHat, xG, '-', 'Color', col, 'LineWidth', 2.0);
                    plot(axBot, lowB, xG, '--', 'Color', col, 'LineWidth', 1.0);
                    plot(axBot, hiB,  xG, '--', 'Color', col, 'LineWidth', 1.0);
                end
            end

            hold(axBot,'off');

            if ~isempty(axTop)
                linkaxes([axTop, axBot],'x');
            end
        end

        if ~isempty(axTop)
            fig = ancestor(axTop, 'figure');
        elseif ~isempty(axBot)
            fig = ancestor(axBot, 'figure');
        end

    % Case 2: make standalone figure
    else
        fig = figure('Color','w', 'Name','Global ridge detection (banded)');
        tl  = tiledlayout(fig, 2, 1);
        tl.TileSpacing = 'compact';
        tl.Padding     = 'compact';

        % Top: stacked reordered spectra
        axTop = nexttile(tl, 1);
        plot(axTop, ppm, X_stack', 'LineWidth', 1.0);
        set(axTop, 'XDir','reverse', 'FontName',fontName, 'FontSize',tickSize);
        ylabel(axTop, 'Intensity (offset)', 'FontName',fontName,'FontSize',labelSize);
        title(axTop, 'Stacked reordered spectra', 'FontName',fontName,'FontSize',titleSize);
        grid(axTop,'on');

        % Bottom: all peaks + ridges
        axBot = nexttile(tl, 2);

        scatter(axBot, allPPM, allSampId, 10, [0.6 0.6 0.6], 'filled', ...
            'MarkerFaceAlpha',0.3, 'MarkerEdgeAlpha',0.3);

        set(axBot, 'XDir','reverse', 'YDir','normal', ...
            'YLim',[0.5 nSamples+0.5], ...
            'FontName',fontName,'FontSize',tickSize);
        xlabel(axBot, sprintf('%s (ppm)', ppmCol), ...
            'FontName',fontName,'FontSize',labelSize);
        ylabel(axBot, 'Sample index (reordered)', ...
            'FontName',fontName,'FontSize',labelSize);
        title(axBot, sprintf('Global ridges (%s, banded)', curveType), ...
            'FontName',fontName,'FontSize',titleSize);
        grid(axBot,'on');
        hold(axBot,'on');

        cmap = lines(max(1, height(ridgeSummaryTable)));

        for i = 1:height(ridgeSummaryTable)
            rid = ridgeSummaryTable.RidgeID(i);
            col = cmap(i,:);
            rows = ridgePointsTable.RidgeID == rid;

            sampR = ridgePointsTable.SampleOrderId(rows);
            ppmR  = ridgePointsTable.PeakPPM(rows);

            scatter(axBot, ppmR, sampR, 36, col, 'filled', 'MarkerFaceAlpha',0.9);

            rcId = find([ridgeCurves.RidgeID] == rid, 1);
            if ~isempty(rcId)
                xG   = ridgeCurves(rcId).xGrid;
                yHat = ridgeCurves(rcId).yHat;
                lowB = ridgeCurves(rcId).lowB;
                hiB  = ridgeCurves(rcId).hiB;

                plot(axBot, yHat, xG, '-', 'Color', col, 'LineWidth', 2.0);
                plot(axBot, lowB, xG, '--', 'Color', col, 'LineWidth', 1.0);
                plot(axBot, hiB,  xG, '--', 'Color', col, 'LineWidth', 1.0);
            end
        end

        hold(axBot,'off');
        linkaxes([axTop, axBot],'x');
    end
end

% ---------- Package outputs ----------
out = struct();
out.X_reordered       = X_reo;
out.order             = order;
out.ppm               = ppm;
out.ridgePointsTable  = ridgePointsTable;
out.ridgeSummaryTable = ridgeSummaryTable;
out.ridgeCurves       = ridgeCurves;
out.fig               = fig;
out.axTop             = axTop;
out.axBot             = axBot;

end

% =====================================================================
% Local helper: DBSCAN-based noise removal along a ridge
% (same spirit as ReorderAlign_Cluster)
% =====================================================================
function [keepMask, cls] = local_dbscan_clean(x, y)
    x = x(:);
    y = y(:);
    n = numel(x);

    if n < 3
        keepMask = true(n,1);
        cls      = ones(n,1);
        return;
    end

    ClsMat = [x, y];

    % z-score
    mu  = mean(ClsMat,1);
    sig = std(ClsMat,0,1);
    sig(sig==0) = 1;
    ClsScaled = (ClsMat - mu) ./ sig;

    d = size(ClsScaled,2);
    rangeProd = exp(sum(log(max(ClsScaled)-min(ClsScaled))) / d);
    epdist = rangeProd * (3*sqrt(pi*d)/n)^(1/d) * sqrt(d/(2*exp(1)*pi));

    cls = dbscan(ClsScaled, epdist, 3);

    u = unique(cls);
    u(u <= 0) = [];

    if isempty(u)
        keepMask = false(n,1);
        return;
    end

    counts = arrayfun(@(c) sum(cls==c), u);
    [~, idMax] = max(counts);
    mainCluster = u(idMax);

    keepMask = (cls == mainCluster);
end

% =====================================================================
% Local helper: fit chosen curve type to (x,y)
% =====================================================================
function [yFit, fitFcn] = local_fit_curve(x, y, curveType)
    x = x(:);
    y = y(:);
    n = numel(x);

    curveType = lower(string(curveType));
    yFit   = NaN(size(x));
    fitFcn = @(t) NaN(size(t));

    switch curveType
        case "linear"
            if n >= 2
                p = polyfit(x, y, 1);
                yFit = polyval(p, x);
                fitFcn = @(t) polyval(p, t);
            end

        case "poly3"
            if n >= 6
                p = polyfit(x, y, 3);
                yFit = polyval(p, x);
                fitFcn = @(t) polyval(p, t);
            elseif n >= 2
                % Fallback to linear if too few points
                p = polyfit(x, y, 1);
                yFit = polyval(p, x);
                fitFcn = @(t) polyval(p, t);
            end

        case "exp"
            if n >= 4
                try
                    ft = fittype('exp1');
                    expFit = fit(x, y, ft);
                    yFit   = expFit(x);
                    fitFcn = @(t) expFit(t);
                catch
                    if n >= 2
                        p = polyfit(x, y, 1);
                        yFit = polyval(p, x);
                        fitFcn = @(t) polyval(p, t);
                    end
                end
            elseif n >= 2
                p = polyfit(x, y, 1);
                yFit = polyval(p, x);
                fitFcn = @(t) polyval(p, t);
            end

        otherwise
            error('Unsupported CurveType "%s". Use "linear", "poly3", or "exp".', curveType);
    end
end
