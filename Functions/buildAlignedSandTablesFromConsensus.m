function [dataOut, alignmentSummary, featureTable] = buildAlignedSandTablesFromConsensus(dataIn, uniqueBins, promotedIncomplete, varargin)
% buildAlignedSandTablesFromConsensus
% -------------------------------------------------------------------------
% Builds Aligned_Table, Misaligned_Table, and All_Table from consensus bins.
%
% Intended workflow:
%   1) Run AlignBinSlopeConsensus_v2.
%   2) Visually validate incomplete candidates with validateIncompleteConsensusBins_v2.
%   3) Use UniqueBins + promotedIncomplete to create final aligned SAND tables.
%
% This function preserves the original ppm column, e.g. freq_ppm, and stores
% the aligned ppm in a separate column named Aligned_PPM.
%
% Outputs per sample:
%   .Aligned_Table      peaks selected by consensus/promoted bins
%   .Misaligned_Table   peaks not selected
%   .All_Table          full peak table with alignment annotations
%
% Example:
% promotedIncomplete = validatedIncomplete(validatedIncomplete.PromoteToAligned == true, :);
%
% [syntheticUrine.alignedData, alignmentSummary, featureTable] = ...
%     buildAlignedSandTablesFromConsensus(SandTable, ...
%         consensusOut.UniqueBins, ...
%         promotedIncomplete, ...
%         'PeakTableField', 'pHITRemain', ...
%         'PPMColName', 'freq_ppm', ...
%         'AlignedTableName', 'Aligned_Table', ...
%         'MisalignedTableName', 'Misaligned_Table', ...
%         'AllTableName', 'All_Table');
% -------------------------------------------------------------------------

p = inputParser;
addRequired(p, 'dataIn');
addRequired(p, 'uniqueBins');
addRequired(p, 'promotedIncomplete');

addParameter(p, 'PeakTableField', 'pHITRemain');
addParameter(p, 'PPMColName', 'freq_ppm');
addParameter(p, 'AlignedPPMCol', 'Aligned_PPM');
addParameter(p, 'AlignedFlagCol', 'Aligned');
addParameter(p, 'FeatureIDCol', 'AlignmentFeatureID');
addParameter(p, 'SourceCol', 'AlignmentSource');
addParameter(p, 'RegionTypeCol', 'AlignmentRegionType');

addParameter(p, 'AlignedTableName', 'Aligned_Table');
addParameter(p, 'MisalignedTableName', 'Misaligned_Table');
addParameter(p, 'AllTableName', 'All_Table');

addParameter(p, 'AlignStatistic', 'median');    % 'median', 'mean', or 'center'
addParameter(p, 'MultiplePeakRule', 'closest'); % 'closest', 'highest_amplitude', or 'skip'
addParameter(p, 'AmplitudeColName', 'amplitude');
addParameter(p, 'RequirePromoteFlag', true);    % only use promotedIncomplete rows where PromoteToAligned == true
addParameter(p, 'UseOnlyUnassignedPeaks', true);
addParameter(p, 'Verbose', true);

parse(p, dataIn, uniqueBins, promotedIncomplete, varargin{:});

peakField      = char(p.Results.PeakTableField);
ppmCol         = char(p.Results.PPMColName);
alignedPPMCol  = char(p.Results.AlignedPPMCol);
alignedFlagCol = char(p.Results.AlignedFlagCol);
featureIDCol   = char(p.Results.FeatureIDCol);
sourceCol      = char(p.Results.SourceCol);
regionTypeCol  = char(p.Results.RegionTypeCol);

alignedTableName    = char(p.Results.AlignedTableName);
misalignedTableName = char(p.Results.MisalignedTableName);
allTableName        = char(p.Results.AllTableName);

alignStatistic = lower(string(p.Results.AlignStatistic));
multipleRule   = lower(string(p.Results.MultiplePeakRule));
ampCol         = char(p.Results.AmplitudeColName);
verbose        = p.Results.Verbose;

% -------------------------------------------------------------------------
% Convert master input to struct array for row-wise editing.
% -------------------------------------------------------------------------
inputWasTable = istable(dataIn);
if inputWasTable
    data = table2struct(dataIn);
else
    data = dataIn;
end

nSamples = numel(data);

if nSamples == 0
    error('dataIn appears to contain no samples.');
end

% -------------------------------------------------------------------------
% Build one standardized feature table from UniqueBins + promotedIncomplete.
% -------------------------------------------------------------------------
featureTable = standardizeFeatureTables(uniqueBins, promotedIncomplete, p.Results.RequirePromoteFlag);

if isempty(featureTable)
    warning('No features were provided after filtering. Returning input with initialized alignment tables.');
end

% Sort by ppm center so the resulting feature IDs are stable.
if ~isempty(featureTable)
    featureTable = sortrows(featureTable, 'FeatureCenterPPM', 'ascend');
    featureTable.FeatureID = (1:height(featureTable))';
end

% -------------------------------------------------------------------------
% Initialize each sample table.
% -------------------------------------------------------------------------
for s = 1:nSamples
    tbl = getNestedTable(data(s), peakField);

    if ~istable(tbl)
        error('Sample %d field %s is not a table.', s, peakField);
    end

    if ~ismember(ppmCol, tbl.Properties.VariableNames)
        error('PPM column "%s" not found in sample %d table "%s".', ppmCol, s, peakField);
    end

    nRows = height(tbl);

    tbl.(alignedFlagCol) = zeros(nRows, 1);
    tbl.(alignedPPMCol)  = tbl.(ppmCol);
    tbl.(featureIDCol)   = strings(nRows, 1);
    tbl.(sourceCol)      = strings(nRows, 1);
    tbl.(regionTypeCol)  = strings(nRows, 1);

    data(s) = setNestedTable(data(s), peakField, tbl);
end

% -------------------------------------------------------------------------
% Assign peaks feature-by-feature.
% -------------------------------------------------------------------------
alignmentSummary = table();

for f = 1:height(featureTable)

    fLow    = featureTable.FeatureStartPPM(f);
    fHigh   = featureTable.FeatureEndPPM(f);
    fCenter = featureTable.FeatureCenterPPM(f);

    rangeLow  = min(fLow, fHigh);
    rangeHigh = max(fLow, fHigh);

    chosenSampleIDs = [];
    chosenRowIDs    = [];
    chosenPPMs      = [];
    multiPeakSamples = [];
    skippedSamples   = [];

    % First pass: identify one candidate peak per sample.
    for s = 1:nSamples
        tbl = getNestedTable(data(s), peakField);
        ppmVals = tbl.(ppmCol);

        inRange = ppmVals >= rangeLow & ppmVals <= rangeHigh;

        if p.Results.UseOnlyUnassignedPeaks
            inRange = inRange & tbl.(alignedFlagCol) == 0;
        end

        id = find(inRange);

        if isempty(id)
            continue;
        elseif numel(id) == 1
            chosenId = id;
        else
            multiPeakSamples(end+1,1) = s; %#ok<AGROW>

            if multipleRule == "closest"
                [~, j] = min(abs(ppmVals(id) - fCenter));
                chosenId = id(j);
            elseif multipleRule == "highest_amplitude"
                if ismember(ampCol, tbl.Properties.VariableNames)
                    [~, j] = max(tbl.(ampCol)(id));
                    chosenId = id(j);
                else
                    [~, j] = min(abs(ppmVals(id) - fCenter));
                    chosenId = id(j);
                end
            elseif multipleRule == "skip"
                skippedSamples(end+1,1) = s; %#ok<AGROW>
                continue;
            else
                error('MultiplePeakRule must be closest, highest_amplitude, or skip.');
            end
        end

        chosenSampleIDs(end+1,1) = s; %#ok<AGROW>
        chosenRowIDs(end+1,1)    = chosenId; %#ok<AGROW>
        chosenPPMs(end+1,1)      = ppmVals(chosenId); %#ok<AGROW>
    end

    if isempty(chosenPPMs)
        if verbose
            fprintf('Feature %d skipped: no peaks found in %.5f-%.5f ppm.\n', ...
                featureTable.FeatureID(f), rangeLow, rangeHigh);
        end
        alignmentSummary = [alignmentSummary; makeSummaryRow(featureTable(f,:), NaN, [], [], [], multiPeakSamples, skippedSamples)]; %#ok<AGROW>
        continue;
    end

    % Determine the aligned ppm for this feature.
    switch alignStatistic
        case "median"
            alignedPPM = median(chosenPPMs, 'omitnan');
        case "mean"
            alignedPPM = mean(chosenPPMs, 'omitnan');
        case "center"
            alignedPPM = fCenter;
        otherwise
            error('AlignStatistic must be median, mean, or center.');
    end

    % Second pass: commit assignments.
    for k = 1:numel(chosenSampleIDs)
        s     = chosenSampleIDs(k);
        rowID = chosenRowIDs(k);

        tbl = getNestedTable(data(s), peakField);

        tbl.(alignedFlagCol)(rowID) = 1;
        tbl.(alignedPPMCol)(rowID)  = alignedPPM;
        tbl.(featureIDCol)(rowID)   = string(featureTable.FeatureID(f));
        tbl.(sourceCol)(rowID)      = string(featureTable.FeatureSource(f));
        tbl.(regionTypeCol)(rowID)  = string(featureTable.FeatureRegionType(f));

        data(s) = setNestedTable(data(s), peakField, tbl);
    end

    alignmentSummary = [alignmentSummary; makeSummaryRow(featureTable(f,:), alignedPPM, chosenSampleIDs, chosenRowIDs, chosenPPMs, multiPeakSamples, skippedSamples)]; %#ok<AGROW>

    if verbose
        fprintf('Feature %d | %s | %.5f-%.5f ppm | %d samples | aligned ppm = %.5f\n', ...
            featureTable.FeatureID(f), string(featureTable.FeatureSource(f)), rangeLow, rangeHigh, numel(chosenSampleIDs), alignedPPM);
    end
end

% -------------------------------------------------------------------------
% Write final Aligned_Table, Misaligned_Table, and All_Table per sample.
% -------------------------------------------------------------------------
for s = 1:nSamples
    tbl = getNestedTable(data(s), peakField);

    data(s).(alignedTableName)    = tbl(tbl.(alignedFlagCol) == 1, :);
    data(s).(misalignedTableName) = tbl(tbl.(alignedFlagCol) == 0, :);
    data(s).(allTableName)        = tbl;
end

if inputWasTable
    dataOut = struct2table(data);
else
    dataOut = data;
end

end

% =========================================================================
% Helper: combine UniqueBins and promoted incomplete rows into one table.
% =========================================================================
function featureTable = standardizeFeatureTables(uniqueBins, promotedIncomplete, requirePromoteFlag)

featureTable = table();

% Unique consensus bins: high-confidence aligned/moderate features.
if ~isempty(uniqueBins) && istable(uniqueBins) && height(uniqueBins) > 0
    u = uniqueBins;
    tmp = table();
    [tmp.FeatureStartPPM, tmp.FeatureEndPPM, tmp.FeatureCenterPPM] = getRangeColumns(u, "unique");
    tmp.FeatureSource = repmat("UniqueBins", height(u), 1);
    tmp.FeatureRegionType = getRegionLabel(u, "RegionType", "unique_consensus");
    tmp.SourceRow = (1:height(u))';
    featureTable = [featureTable; tmp]; %#ok<AGROW>
end

% Promoted incomplete candidates: visually reviewed features.
if ~isempty(promotedIncomplete) && istable(promotedIncomplete) && height(promotedIncomplete) > 0
    inc = promotedIncomplete;

    if requirePromoteFlag && ismember('PromoteToAligned', inc.Properties.VariableNames)
        inc = inc(inc.PromoteToAligned == true, :);
    end

    if height(inc) > 0
        tmp = table();
        [tmp.FeatureStartPPM, tmp.FeatureEndPPM, tmp.FeatureCenterPPM] = getRangeColumns(inc, "incomplete");
        tmp.FeatureSource = repmat("PromotedIncomplete", height(inc), 1);

        if ismember('ValidationRegionType', inc.Properties.VariableNames)
            tmp.FeatureRegionType = string(inc.ValidationRegionType);
        elseif ismember('AdjustedRegionType_Incomplete', inc.Properties.VariableNames)
            tmp.FeatureRegionType = string(inc.AdjustedRegionType_Incomplete);
        elseif ismember('RegionType', inc.Properties.VariableNames)
            tmp.FeatureRegionType = string(inc.RegionType);
        else
            tmp.FeatureRegionType = repmat("promoted_incomplete", height(inc), 1);
        end

        tmp.SourceRow = (1:height(inc))';
        featureTable = [featureTable; tmp]; %#ok<AGROW>
    end
end

if ~isempty(featureTable)
    good = isfinite(featureTable.FeatureStartPPM) & isfinite(featureTable.FeatureEndPPM);
    featureTable = featureTable(good, :);
end

end

% =========================================================================
% Helper: choose range columns by priority.
% =========================================================================
function [startPPM, endPPM, centerPPM] = getRangeColumns(tbl, tableType)

vars = tbl.Properties.VariableNames;

if strcmp(tableType, "incomplete")
    priority = { ...
        {'ValidatedBinStart','ValidatedBinEnd','ValidatedBinCenter'}, ...
        {'AdjustedBinStart_Incomplete','AdjustedBinEnd_Incomplete','AdjustedBinCenter_Incomplete'}, ...
        {'ConsensusBinStart','ConsensusBinEnd','ConsensusBinCenter'}, ...
        {'BinStart','BinEnd','BinCenter'} ...
        };
else
    priority = { ...
        {'ConsensusBinStart','ConsensusBinEnd','ConsensusBinCenter'}, ...
        {'BinStart','BinEnd','BinCenter'}, ...
        {'ValidatedBinStart','ValidatedBinEnd','ValidatedBinCenter'} ...
        };
end

startPPM = [];
endPPM = [];
centerPPM = [];

for k = 1:numel(priority)
    cols = priority{k};
    if all(ismember(cols(1:2), vars))
        startPPM = tbl.(cols{1});
        endPPM   = tbl.(cols{2});

        if numel(cols) >= 3 && ismember(cols{3}, vars)
            centerPPM = tbl.(cols{3});
        else
            centerPPM = mean([startPPM, endPPM], 2, 'omitnan');
        end
        return;
    end
end

if ismember('RefPPM', vars)
    startPPM = tbl.RefPPM;
    endPPM = tbl.RefPPM;
    centerPPM = tbl.RefPPM;
    return;
end

error('Could not determine feature ppm range. Expected Validated, Adjusted, Consensus, or BinStart/BinEnd columns.');

end

% =========================================================================
% Helper: region labels.
% =========================================================================
function labels = getRegionLabel(tbl, preferredCol, fallback)
if ismember(preferredCol, tbl.Properties.VariableNames)
    labels = string(tbl.(preferredCol));
else
    labels = repmat(string(fallback), height(tbl), 1);
end
end

% =========================================================================
% Helper: table field access, including cell-wrapped tables.
% =========================================================================
function tbl = getNestedTable(sampleStruct, fieldName)
if ~isfield(sampleStruct, fieldName)
    error('Field "%s" not found in sample structure.', fieldName);
end
value = sampleStruct.(fieldName);
if iscell(value)
    if isempty(value)
        tbl = table();
    else
        tbl = value{1};
    end
else
    tbl = value;
end
end

function sampleStruct = setNestedTable(sampleStruct, fieldName, tbl)
oldValue = sampleStruct.(fieldName);
if iscell(oldValue)
    sampleStruct.(fieldName) = {tbl};
else
    sampleStruct.(fieldName) = tbl;
end
end

% =========================================================================
% Helper: build one summary row.
% =========================================================================
function row = makeSummaryRow(featureRow, alignedPPM, sampleIDs, rowIDs, originalPPMs, multiPeakSamples, skippedSamples)
row = table();
row.FeatureID = featureRow.FeatureID;
row.FeatureSource = featureRow.FeatureSource;
row.FeatureRegionType = featureRow.FeatureRegionType;
row.FeatureStartPPM = featureRow.FeatureStartPPM;
row.FeatureEndPPM = featureRow.FeatureEndPPM;
row.FeatureCenterPPM = featureRow.FeatureCenterPPM;
row.AlignedPPM = alignedPPM;
row.NumSamplesAligned = numel(sampleIDs);
row.SampleIDs = {sampleIDs};
row.RowIDs = {rowIDs};
row.OriginalPPMs = {originalPPMs};
row.MultiPeakSamples = {multiPeakSamples};
row.SkippedSamples = {skippedSamples};
end
