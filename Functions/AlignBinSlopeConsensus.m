function consensusOut = AlignBinSlopeConsensus_v2(dataIn, varargin)
% AlignBinSlopeConsensus_v2
% -------------------------------------------------------------------------
% Wrapper around AlignBinSlope that runs multiple center samples, combines
% metadata, and separates:
%
%   1) UniqueBins
%      Well-aligned consensus bins only. By default this keeps only
%      RegionType = aligned or moderate.
%
%   2) IncompleteCandidates
%      Incomplete bins rescored using the FIRST bin width where the maximum
%      sample coverage was reached, instead of using MaxBinSize/STOP width.

%
% Example:
% consensusOut = AlignBinSlopeConsensus_v2(SandTable, ...
%     'PeakTableField', 'pHITRemain', ...
%     'PPMColName', 'freq_ppm', ...
%     'InitialBinSize', 0.001, ...
%     'MaxBinSize', 0.02, ...
%     'BinIncrement', 0.001, ...
%     'AmplitudeThreshold', 0.05, ...
%     'NormLowThreshold', 0.10, ...
%     'RegionTypesToKeep', {'aligned','moderate'}, ...
%     'MinUniqueCount', 5, ...
%     'MergeTolerancePPM', 0.002, ...
%     'Verbose', true);
%
% Outputs:
%   consensusOut.AllMeta              all rows from all center-sample runs
%   consensusOut.CandidateMeta        aligned/moderate rows used for UniqueBins
%   consensusOut.UniqueBins           merged well-aligned consensus bins
%   consensusOut.IncompleteMeta       all incomplete rows after min count filter
%   consensusOut.IncompleteCandidates merged/rescored incomplete candidate bins
%   consensusOut.Parameters           parameters used
%
% -------------------------------------------------------------------------

p = inputParser;
addRequired(p, 'dataIn');

% Core table options
addParameter(p, 'PeakTableField', 'Peaks');
addParameter(p, 'PPMColName', '', @(s) ischar(s) || isstring(s));
addParameter(p, 'CenterSampleIDs', []);

% Consensus filtering
addParameter(p, 'RegionTypesToKeep', {'aligned','moderate'}); % for UniqueBins only
addParameter(p, 'IncompleteRegionTypesToKeep', {'incomplete'});
addParameter(p, 'MinUniqueCount', 3);
addParameter(p, 'MergeTolerancePPM', 0.002);

% Incomplete rescoring behavior
addParameter(p, 'IncompleteMinAdjustedNormSlope', []); % optional filter; [] keeps all
addParameter(p, 'UseAdjustedIncompleteRange', true);   % use best-coverage width for incomplete table ranges

% Store run outputs?
addParameter(p, 'KeepCenterRuns', false);
addParameter(p, 'Verbose', true);

% Pass-through AlignBinSlope parameters
addParameter(p, 'InitialBinSize', 0.001);
addParameter(p, 'MaxBinSize', 0.02);
addParameter(p, 'BinIncrement', 0.001);
addParameter(p, 'AmplitudeThreshold', 0.05);
addParameter(p, 'NormHighThreshold', 0.70);
addParameter(p, 'NormLowThreshold', 0.30);

parse(p, dataIn, varargin{:});

peakField      = char(p.Results.PeakTableField);
ppmColName     = char(p.Results.PPMColName);
centerIDs      = p.Results.CenterSampleIDs;
keepTypes      = string(p.Results.RegionTypesToKeep);
incompleteKeep = string(p.Results.IncompleteRegionTypesToKeep);
mergeTol       = p.Results.MergeTolerancePPM;
minUniqueCount = p.Results.MinUniqueCount;
keepRuns       = p.Results.KeepCenterRuns;
verbose        = p.Results.Verbose;
initialBinSize = p.Results.InitialBinSize;
binIncrement   = p.Results.BinIncrement;
normHighThr    = p.Results.NormHighThreshold;
normLowThr     = p.Results.NormLowThreshold;
useAdjustedIncompleteRange = p.Results.UseAdjustedIncompleteRange;
incompleteMinAdjustedNormSlope = p.Results.IncompleteMinAdjustedNormSlope;

% Convert input if needed only to count samples.
if istable(dataIn)
    dataStruct = table2struct(dataIn);
else
    dataStruct = dataIn;
end

nSamples = numel(dataStruct);

if isempty(centerIDs)
    centerIDs = 1:nSamples;
end

allMeta = table();
centerRuns = cell(numel(centerIDs),1);

for i = 1:numel(centerIDs)

    cID = centerIDs(i);

    if verbose
        fprintf('\nRunning AlignBinSlope with center sample %d of %d...\n', cID, nSamples);
    end

    [dataTemp, metaTemp] = AlignBinSlope(dataIn, ...
        'PeakTableField', peakField, ...
        'CenterSampleID', cID, ...
        'PPMColName', ppmColName, ...
        'InitialBinSize', initialBinSize, ...
        'MaxBinSize', p.Results.MaxBinSize, ...
        'BinIncrement', binIncrement, ...
        'AmplitudeThreshold', p.Results.AmplitudeThreshold, ...
        'NormHighThreshold', normHighThr, ...
        'NormLowThreshold', normLowThr, ...
        'Verbose', false);

    metaTemp.CenterSampleID = repmat(cID, height(metaTemp), 1);

    % Add derived diagnostics needed for consensus and incomplete rescoring.
    metaTemp = addConsensusDiagnostics(metaTemp, nSamples, initialBinSize, binIncrement, normHighThr, normLowThr);

    allMeta = [allMeta; metaTemp]; %#ok<AGROW>

    if keepRuns
        centerRuns{i} = dataTemp;
    end
end

% Build center column for all rows, based on stop width.
allMeta.BinCenter = (allMeta.BinStart + allMeta.BinEnd) ./ 2;

% -------------------------------------------------------------------------
% 1) Well-aligned consensus bins: aligned/moderate only by default
% -------------------------------------------------------------------------
wellMask = ismember(string(allMeta.RegionType), keepTypes) & ...
           allMeta.StopUniqueCount >= minUniqueCount;

candidateMeta = allMeta(wellMask, :);

if verbose
    fprintf('\nWell-aligned candidate bins before merging: %d\n', height(candidateMeta));
end

if ~isempty(candidateMeta)
    candidateMeta = sortrows(candidateMeta, 'BinCenter', 'ascend');
end

uniqueBins = mergeConsensusBins(candidateMeta, mergeTol, false);

if verbose
    fprintf('Unique well-aligned consensus bins after merging: %d\n', height(uniqueBins));
end

% -------------------------------------------------------------------------
% 2) Incomplete candidates: score at first width where max coverage occurs
% -------------------------------------------------------------------------
incompleteMask = ismember(string(allMeta.RegionType), incompleteKeep) & ...
                 allMeta.BestCoverageUniqueCount >= minUniqueCount;

incompleteMeta = allMeta(incompleteMask, :);

if ~isempty(incompleteMinAdjustedNormSlope)
    incompleteMeta = incompleteMeta(incompleteMeta.AdjustedNormSlope_Incomplete >= incompleteMinAdjustedNormSlope, :);
end

% For incomplete rows, optionally use adjusted/best-coverage range instead
% of the full max-width incomplete range. This makes validation tighter.
if useAdjustedIncompleteRange && ~isempty(incompleteMeta)
    incompleteMeta.OriginalIncompleteBinStart = incompleteMeta.BinStart;
    incompleteMeta.OriginalIncompleteBinEnd   = incompleteMeta.BinEnd;
    incompleteMeta.OriginalIncompleteBinCenter = incompleteMeta.BinCenter;

    incompleteMeta.BinStart  = incompleteMeta.AdjustedBinStart_Incomplete;
    incompleteMeta.BinEnd    = incompleteMeta.AdjustedBinEnd_Incomplete;
    incompleteMeta.BinCenter = incompleteMeta.AdjustedBinCenter_Incomplete;
end

if verbose
    fprintf('Incomplete candidate bins before merging: %d\n', height(incompleteMeta));
end

if ~isempty(incompleteMeta)
    incompleteMeta = sortrows(incompleteMeta, 'BinCenter', 'ascend');
end

incompleteCandidates = mergeConsensusBins(incompleteMeta, mergeTol, true);

if verbose
    fprintf('Merged incomplete candidates after adjusted rescoring: %d\n', height(incompleteCandidates));
end

% -------------------------------------------------------------------------
% Output
% -------------------------------------------------------------------------
consensusOut = struct;
consensusOut.AllMeta = allMeta;
consensusOut.CandidateMeta = candidateMeta;
consensusOut.UniqueBins = uniqueBins;
consensusOut.IncompleteMeta = incompleteMeta;
consensusOut.IncompleteCandidates = incompleteCandidates;

consensusOut.Parameters = struct;
consensusOut.Parameters.PeakTableField = peakField;
consensusOut.Parameters.PPMColName = ppmColName;
consensusOut.Parameters.CenterSampleIDs = centerIDs;
consensusOut.Parameters.RegionTypesToKeep = keepTypes;
consensusOut.Parameters.IncompleteRegionTypesToKeep = incompleteKeep;
consensusOut.Parameters.MinUniqueCount = minUniqueCount;
consensusOut.Parameters.MergeTolerancePPM = mergeTol;
consensusOut.Parameters.InitialBinSize = initialBinSize;
consensusOut.Parameters.MaxBinSize = p.Results.MaxBinSize;
consensusOut.Parameters.BinIncrement = binIncrement;
consensusOut.Parameters.AmplitudeThreshold = p.Results.AmplitudeThreshold;
consensusOut.Parameters.NormHighThreshold = normHighThr;
consensusOut.Parameters.NormLowThreshold = normLowThr;
consensusOut.Parameters.UseAdjustedIncompleteRange = useAdjustedIncompleteRange;
consensusOut.Parameters.IncompleteMinAdjustedNormSlope = incompleteMinAdjustedNormSlope;

if keepRuns
    consensusOut.CenterRuns = centerRuns;
end

end

% =========================================================================
% Helper: add diagnostics to every AlignBinSlope metadata row
% =========================================================================
function metaTemp = addConsensusDiagnostics(metaTemp, nSamples, initialBinSize, binIncrement, normHighThr, normLowThr)

nRows = height(metaTemp);

stopUnique = zeros(nRows,1);
bestID = zeros(nRows,1);
bestWidth = zeros(nRows,1);
bestUnique = zeros(nRows,1);
bestCounts = zeros(nRows,1);
adjustedSlope = zeros(nRows,1);
adjustedNormSlope = zeros(nRows,1);
adjustedBinStart = zeros(nRows,1);
adjustedBinEnd = zeros(nRows,1);
adjustedBinCenter = zeros(nRows,1);
adjustedRegion = strings(nRows,1);

maxPossibleSlope = nSamples / initialBinSize;

for r = 1:nRows

    widths = metaTemp.BinWidths{r};
    uniqueCounts = metaTemp.UniqueCount{r};
    counts = metaTemp.Counts{r};
    refPPM = metaTemp.RefPPM(r);
    stopID = metaTemp.StopID(r);

    % Protect against accidental mismatch/truncation.
    stopID = min(stopID, numel(uniqueCounts));

    stopUnique(r) = uniqueCounts(stopID);

    % Best coverage should be computed only across attempted widths up to STOP.
    attemptedIDs = 1:stopID;
    attemptedUnique = uniqueCounts(attemptedIDs);

    maxU = max(attemptedUnique);

    % "As soon as it reached max samples" = first width where max coverage occurs.
    firstBestLocal = find(attemptedUnique == maxU, 1, 'first');
    bestID(r) = attemptedIDs(firstBestLocal);

    bestWidth(r) = widths(bestID(r));
    bestUnique(r) = uniqueCounts(bestID(r));
    bestCounts(r) = counts(bestID(r));

    % Adjusted slope:
    % If best coverage occurs at the initial bin, use width itself to avoid
    % division by ~0. Otherwise use added width, matching AlignBinSlope logic.
    if bestID(r) == 1
        denom = max(bestWidth(r), eps);
    else
        denom = max(bestWidth(r) - initialBinSize, binIncrement);
    end

    adjustedSlope(r) = bestUnique(r) / denom;
    adjustedNormSlope(r) = min(1, max(0, adjustedSlope(r) / maxPossibleSlope));

    adjustedBinStart(r) = refPPM - bestWidth(r)/2;
    adjustedBinEnd(r) = refPPM + bestWidth(r)/2;
    adjustedBinCenter(r) = refPPM;

    if adjustedNormSlope(r) >= normHighThr
        adjustedRegion(r) = "incomplete_aligned_candidate";
    elseif adjustedNormSlope(r) >= normLowThr
        adjustedRegion(r) = "incomplete_moderate_candidate";
    else
        adjustedRegion(r) = "incomplete_weak_candidate";
    end
end

metaTemp.StopUniqueCount = stopUnique;
metaTemp.BestCoverageID_Incomplete = bestID;
metaTemp.BestCoverageWidth_Incomplete = bestWidth;
metaTemp.BestCoverageUniqueCount = bestUnique;
metaTemp.BestCoverageTotalCount = bestCounts;
metaTemp.AdjustedSlope_Incomplete = adjustedSlope;
metaTemp.AdjustedNormSlope_Incomplete = adjustedNormSlope;
metaTemp.AdjustedRegionType_Incomplete = adjustedRegion;
metaTemp.AdjustedBinStart_Incomplete = adjustedBinStart;
metaTemp.AdjustedBinEnd_Incomplete = adjustedBinEnd;
metaTemp.AdjustedBinCenter_Incomplete = adjustedBinCenter;

end

% =========================================================================
% Helper: merge redundant bins and keep best representative row
% =========================================================================
function uniqueBins = mergeConsensusBins(candidateMeta, mergeTol, useIncompleteAdjustedRanking)

if isempty(candidateMeta)
    uniqueBins = candidateMeta;
    return;
end

candidateMeta = sortrows(candidateMeta, 'BinCenter', 'ascend');

n = height(candidateMeta);
assigned = false(n,1);
groups = {};

for i = 1:n
    if assigned(i)
        continue;
    end

    currentCenter = candidateMeta.BinCenter(i);
    groupIDs = find(abs(candidateMeta.BinCenter - currentCenter) <= mergeTol);

    % Expand group based on overlapping bin boundaries.
    binStart = min(candidateMeta.BinStart(groupIDs));
    binEnd   = max(candidateMeta.BinEnd(groupIDs));

    expanded = true;
    while expanded
        expanded = false;

        overlapIDs = find(candidateMeta.BinStart <= binEnd & candidateMeta.BinEnd >= binStart);
        newIDs = unique([groupIDs; overlapIDs]);

        if numel(newIDs) > numel(groupIDs)
            groupIDs = newIDs;
            binStart = min(candidateMeta.BinStart(groupIDs));
            binEnd   = max(candidateMeta.BinEnd(groupIDs));
            expanded = true;
        end
    end

    assigned(groupIDs) = true;
    groups{end+1} = groupIDs; %#ok<AGROW>
end

bestRows = table();

for g = 1:numel(groups)

    ids = groups{g};
    groupTbl = candidateMeta(ids,:);

    rankTbl = groupTbl;

    if useIncompleteAdjustedRanking && ismember('AdjustedNormSlope_Incomplete', rankTbl.Properties.VariableNames)
        % For incomplete candidates, prefer:
        %   1) highest best-coverage unique count
        %   2) highest adjusted norm slope
        %   3) narrowest adjusted/best-coverage width
        rankTbl.NegRankWidth = -rankTbl.BestCoverageWidth_Incomplete;

        rankTbl = sortrows(rankTbl, ...
            {'BestCoverageUniqueCount','AdjustedNormSlope_Incomplete','NegRankWidth'}, ...
            {'descend','descend','descend'});
    else
        % For aligned/moderate bins, prefer:
        %   1) highest stop unique count
        %   2) highest original norm slope
        %   3) narrowest stop width
        rankTbl.NegRankWidth = -rankTbl.StopWidth_ppm;

        rankTbl = sortrows(rankTbl, ...
            {'StopUniqueCount','NormSlope','NegRankWidth'}, ...
            {'descend','descend','descend'});
    end

    best = rankTbl(1,:);

    best.ConsensusBinStart = min(groupTbl.BinStart);
    best.ConsensusBinEnd   = max(groupTbl.BinEnd);
    best.ConsensusBinCenter = median(groupTbl.BinCenter);
    best.NumMergedBins = height(groupTbl);
    best.CenterSamplesRepresented = {unique(groupTbl.CenterSampleID)};

    if useIncompleteAdjustedRanking
        best.ConsensusBestCoverageUniqueCount = max(groupTbl.BestCoverageUniqueCount);
        best.ConsensusAdjustedNormSlope = max(groupTbl.AdjustedNormSlope_Incomplete);
        best.ConsensusAdjustedRegionTypes = {unique(string(groupTbl.AdjustedRegionType_Incomplete))};
    end

    bestRows = [bestRows; best]; %#ok<AGROW>
end

uniqueBins = bestRows;

if ismember('NegRankWidth', uniqueBins.Properties.VariableNames)
    uniqueBins.NegRankWidth = [];
end

uniqueBins = sortrows(uniqueBins, 'ConsensusBinCenter', 'ascend');

end
