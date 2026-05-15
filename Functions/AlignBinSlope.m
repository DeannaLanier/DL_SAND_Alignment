function [dataOut, sigmoidMetaTable] = AlignBinSlope(dataIn, varargin)
% AlignBinSlope: bin-widening, slope-based alignment on tabular SAND peaks.
%
% Deanna Lanier - Updated October 1, 2025
% Updated: Nov 7, 2025 (add 'PPMColName' parameter)
%
% Stops when:
%   (a) any sample has >1 peak in-bin (overcrowded)  [hard stop]
%   (b) every sample contributes exactly one peak      [success]
%   (c) max bin size reached                          [incomplete]
%
% Slope:
%   MaxPossibleSlope = (#samples) / (InitialBinSize)
%   OverallSlope     = (unique samples included at STOP) / (STOP_width - InitialBinSize)
%   NormSlope        = OverallSlope / MaxPossibleSlope
%
% Peaks with amplitude < AmplitudeThreshold are ignored for alignment.
%
% Outputs per sample:
%   .Aligned_Table      = samples considered aligned (+ adjusted ppm in Aligned_PPM)
%   .Misaligned_Table   = samples considered misaligned (original ppm)
%   .All_Table          = all peaks with Aligned flag and Aligned_PPM
%   sigmoidMetaTable    = per-reference-bin diagnostics
%
% Parameters
%   'PeakTableField'   (char)  : name of field holding per-sample peak table (default 'Peaks')
%   'CenterSampleID'   (int)   : index of reference sample (default 1)
%   'PPMColName'       (char)  : name of ppm column to use. If omitted, auto-detects from
%                                {'Adj_Peak','adj_peak','ppm','freq_ppm'}
%   'InitialBinSize'   (double): ppm (default 0.001)
%   'MaxBinSize'       (double): ppm (default 0.02)
%   'BinIncrement'     (double): ppm (default 0.001)
%   'AmplitudeThreshold' (double): reject peaks below this amplitude (default 0.05)
%   'NormHighThreshold' (double): normalized slope ≥ this → perfectly aligned (default 0.70)
%   'NormLowThreshold'  (double): normalized slope < this → not aligned (default 0.30)
%   'Verbose'           (logical): print progress (default true)

p = inputParser;
addRequired(p, 'dataIn');

addParameter(p, 'PeakTableField', 'Peaks');    % struct containing per-sample peak table
addParameter(p, 'CenterSampleID', 1);          % reference sample (by index)

% NEW: allow caller to choose ppm column
addParameter(p, 'PPMColName', '', @(s) ischar(s) || isstring(s));

addParameter(p, 'InitialBinSize', 0.001);      % ppm
addParameter(p, 'MaxBinSize', 0.02);           % ppm
addParameter(p, 'BinIncrement', 0.001);        % ppm
addParameter(p, 'AmplitudeThreshold', 0.05);   % reject peaks below this amplitude

% Slope classification thresholds on normalized slope (0..1)
addParameter(p, 'NormHighThreshold', 0.70); % Anything above this value = perfectly aligned
addParameter(p, 'NormLowThreshold',  0.30); % Anything below this value not aligned

addParameter(p, 'Verbose', true);

parse(p, dataIn, varargin{:});

peakField        = p.Results.PeakTableField;
center_sample_id = p.Results.CenterSampleID;
ppmColNameParam  = char(p.Results.PPMColName);
initialBinSize   = p.Results.InitialBinSize;
maxBinSize       = p.Results.MaxBinSize;
binIncrement     = p.Results.BinIncrement;
ampThr           = p.Results.AmplitudeThreshold;
normHighThr      = p.Results.NormHighThreshold;
normLowThr       = p.Results.NormLowThreshold;
verbose          = p.Results.Verbose;

% ----------------------- bin widths -------------------------------------
nSteps    = max(2, round((maxBinSize - initialBinSize)/binIncrement) + 1);
binWidths = linspace(initialBinSize, maxBinSize, nSteps);
numWidths = numel(binWidths);

% ----------------------- data input -------------------------------------
if istable(dataIn)
    data = table2struct(dataIn);
    isTableInput = true;
else
    data = dataIn;
    isTableInput = false;
end
num_samples = numel(data);

% ----------------------- Reference peaks (ppm) --------------------------
refTable = data(center_sample_id).(peakField);

% Use requested ppm column if provided, otherwise auto-detect
if ~isempty(ppmColNameParam)
    adj_ppm_col = ensureVarExists(refTable, ppmColNameParam);
else
    adj_ppm_col = matchVar(refTable, {'Adj_Peak','adj_peak','ppm','freq_ppm'});
end
reference_ppm = refTable.(adj_ppm_col);

% Find amplitude column (optional)
amp_col = tryMatchVarOrEmpty(refTable, {'Amplitude','amplitude','amp','AMP','intensity','Intensity'});
hasAmp  = ~isempty(amp_col);

% ----------------------- Init per-sample bookkeeping --------------------
used = cell(num_samples, 1);
ppmValsPerSample = cell(num_samples,1);
ampValsPerSample = cell(num_samples,1);

for s = 1:num_samples
    tbl = data(s).(peakField);

    % Validate ppm column presence for each sample table
    if ~isempty(ppmColNameParam)
        ppm_col_this = ensureVarExists(tbl, ppmColNameParam);
    else
        % Keep the same detected column name across all samples when possible;
        % if missing in a specific sample, attempt auto-detect for that table.
        if ismember(adj_ppm_col, tbl.Properties.VariableNames)
            ppm_col_this = adj_ppm_col;
        else
            ppm_col_this = matchVar(tbl, {'Adj_Peak','adj_peak','ppm','freq_ppm'});
        end
    end

    % Create and/or clear alignment columns
    if ~ismember('Aligned', tbl.Properties.VariableNames)
        tbl.Aligned = zeros(height(tbl),1);
    else
        tbl.Aligned(:) = 0;
    end
    if ~ismember('Aligned_PPM', tbl.Properties.VariableNames)
        tbl.Aligned_PPM = tbl.(ppm_col_this);
    else
        tbl.Aligned_PPM = tbl.(ppm_col_this);
    end

    data(s).(peakField) = tbl;                 % write back
    used{s}             = false(height(tbl),1);
    ppmValsPerSample{s} = tbl.(ppm_col_this);

    if hasAmp && ismember(amp_col, tbl.Properties.VariableNames)
        ampValsPerSample{s} = tbl.(amp_col);
    else
        % no amp col or not present in this table: accept all as above threshold
        ampValsPerSample{s} = inf(height(tbl),1);
    end
end

% ----------------------- Slope ceiling (smallest bin) -------------------
maxPossibleSlope = num_samples / initialBinSize;

% ----------------------- Iterate reference peaks ------------------------
metaCells = cell(numel(reference_ppm),1);

for refID = 1:numel(reference_ppm)
    ref_ppm = reference_ppm(refID);

    totalCount      = zeros(numWidths,1);
    uniqueCount     = zeros(numWidths,1);
    overlapFlag     = false(numWidths,1);
    allPeaksPerBin  = cell(numWidths,1);
    matchIdxPerSmpl = cell(numWidths, num_samples);

    STOP_id = numWidths;
    stop_reason = 'none';

    %  Bin widening
    for w = 1:numWidths
        halfw  = binWidths(w)/2;
        rangeL = ref_ppm - halfw;
        rangeR = ref_ppm + halfw;

        sampleIDs = [];
        peaksPPM  = [];

        for s = 1:num_samples
            ppmVals = ppmValsPerSample{s};
            ampVals = ampValsPerSample{s};

            valid = ~used{s} & (ppmVals >= rangeL) & (ppmVals <= rangeR) & (ampVals >= ampThr);
            id    = find(valid);
            matchIdxPerSmpl{w,s} = id;

            if numel(id) == 1
                sampleIDs = [sampleIDs; s];             
                peaksPPM  = [peaksPPM; ppmVals(id)];    
            elseif numel(id) > 1
                overlapFlag(w) = true; % overcrowded: >1 peak from a single sample
            end
        end

        totalCount(w)      = numel(peaksPPM);
        uniqueCount(w)     = numel(unique(sampleIDs));
        allPeaksPerBin{w}  = peaksPPM;

        % Stop #1: overcrowded (hard stop)
        if overlapFlag(w)
            STOP_id     = w;
            stop_reason = 'overcrowded';
            break
        end

        % Stop #2: all samples contribute exactly one peak
        if uniqueCount(w) == num_samples
            counts = cellfun(@numel, matchIdxPerSmpl(w,:));
            if all(counts == 1)
                STOP_id     = w;
                stop_reason = 'all_included';
                break
            end
        end
    end

    % Stop #3: Incomplete
    if strcmp(stop_reason,'none')
        STOP_id     = numWidths;
        stop_reason = 'incomplete';
    end

    %  Slope computation 
    stopWidth     = binWidths(STOP_id);
    stopUnique    = uniqueCount(STOP_id);
    denom         = max(stopWidth - initialBinSize, eps); % width added since start
    overallSlope  = stopUnique / denom;                   % samples per ppm
    normSlope     = min(1, max(0, overallSlope / maxPossibleSlope));

    %  Region label based on stop # and slope
    if strcmp(stop_reason,'overcrowded')
        regionType = 'overcrowded';
    elseif strcmp(stop_reason,'all_included')
        if normSlope >= normHighThr
            regionType = 'aligned';
        elseif normSlope >= normLowThr
            regionType = 'moderate';
        else
            regionType = 'weak';
        end
    else
        regionType = 'incomplete';
    end

    % Commit alignment when valid 
    if strcmp(stop_reason,'all_included') && normSlope >= normLowThr
        binMidPPM = median(allPeaksPerBin{STOP_id});
        for s = 1:num_samples
            id = matchIdxPerSmpl{STOP_id, s};
            data(s).(peakField).Aligned(id)     = 1;
            data(s).(peakField).Aligned_PPM(id) = binMidPPM;
            used{s}(id) = true;
        end
    end

    % Metadata row 
    m = struct;
    m.RefPPM                  = ref_ppm;
    m.BinStart                = ref_ppm - stopWidth/2;
    m.BinEnd                  = ref_ppm + stopWidth/2;
    m.BinWidths               = binWidths(:);
    m.UniqueCount             = uniqueCount(:);
    m.Counts                  = totalCount(:);
    m.OverlapFlag             = overlapFlag(:);
    m.StopID                  = STOP_id;
    m.StopReason              = stop_reason;
    m.StopWidth_ppm           = stopWidth;
    m.OverallSlope_pk_per_ppm = overallSlope;
    m.MaxPossibleSlope        = maxPossibleSlope;
    m.NormSlope               = normSlope;
    m.NormHighThreshold       = normHighThr;
    m.NormLowThreshold        = normLowThr;
    m.RegionType              = regionType;
    m.AllPeaksAtStop          = allPeaksPerBin{STOP_id};
    metaCells{refID}          = m;

    if verbose
        fprintf('Ref %.6f | stop=%s @ %.5f ppm | uniq=%d | normSlope=%.2f | region=%s\n', ...
            ref_ppm, stop_reason, stopWidth, uniqueCount(STOP_id), normSlope, regionType);
    end
end

sigmoidMetaTable = struct2table([metaCells{:}]);

%  Split outputs per sample
for s = 1:num_samples
    tbl = data(s).(peakField);
    data(s).Aligned_Table    = tbl(tbl.Aligned == 1, :);
    data(s).Misaligned_Table = tbl(tbl.Aligned == 0, :);
    data(s).All_Table        = tbl;
end

if isTableInput
    dataOut = struct2table(data);
else
    dataOut = data;
end
end

% ----------------------- helpers ----------------------------------------
function name = ensureVarExists(tbl, varName)
% Assert a variable exists in table, else throw a descriptive error.
if ~ismember(varName, tbl.Properties.VariableNames)
    error('PPM column "%s" not found. Available columns: %s', ...
        varName, strjoin(tbl.Properties.VariableNames, ', '));
end
name = varName;
end

function name = matchVar(tbl, candidates)
for c = 1:numel(candidates)
    id = find(strcmpi(candidates{c}, tbl.Properties.VariableNames), 1);
    if ~isempty(id)
        name = tbl.Properties.VariableNames{id};
        return;
    end
end
error('Expected one of: %s', strjoin(candidates, ', '));
end

function name = tryMatchVarOrEmpty(tbl, candidates)
for c = 1:numel(candidates)
    id = find(strcmpi(candidates{c}, tbl.Properties.VariableNames), 1);
    if ~isempty(id)
        name = tbl.Properties.VariableNames{id};
        return;
    end
end
name = '';
end
