function [SandTable, processing] = ampFilter_DL(SandTable, varargin)
%{ Deanna Lanier (Adapted from Jojo's original amp_filter function)  

% ampFilter_DL
% Filters peaks by amplitude using a PER-SAMPLE threshold derived from the
% "negative ppm" region (freq_ppm < 0) within EACH sample (by default).
%
% Default per-sample threshold rule:
%   ampThresh(s) = max( amplitude for peaks where freq_ppm < 0 ) for sample s
% If a sample has no peaks in that ppm region, it falls back to 0 for that sample.
%
% Filtering rule (per sample):
%   keep if amplitude >= ampThresh(s)
%   remove if amplitude <  ampThresh(s)
%
% INPUT
%   SandTable : table with per-sample peak tables stored in a column (cell array)
%
% NAME-VALUE PAIRS (optional)
%   'Processing'      : existing processing struct to append into (default: struct())
%   'InputVar'        : input column (default: 'Con_RR')
%   'OutputVar'       : output column (default: 'Con_RR_AF')
%   'AmpThresh'       : numeric scalar (global override) OR numeric vector length nSamples
%                       If empty/not provided -> auto per-sample from freq_ppm < 0
%   'FreqPpmMaskFcn'  : function handle to define "noise region" rows (default: @(ppm) ppm < 0)
%   'MetadataStep'    : metadata step name (default: 'AT')
%
% OUTPUTS
%   SandTable  : updated SandTable with OutputVar column containing filtered peak tables
%   processing : appended under processing.metadata_AT (so it won’t overwrite RR metadata)
%       processing.metadata_AT.table.Amplitude_Removed : removed rows per sample (cell of tables)
%       processing.metadata_AT.info.AT                 : details (threshold(s), counts, etc.)

% -----------------------------
% Parse inputs
% -----------------------------
p = inputParser;
p.addRequired('SandTable', @(x) istable(x));
p.addParameter('Processing',      struct(), @(x) isstruct(x) || isempty(x));
p.addParameter('InputVar',        'Con_RR', @(x) ischar(x) || isstring(x));
p.addParameter('OutputVar',       'Con_RR_AF', @(x) ischar(x) || isstring(x));
p.addParameter('AmpThresh',       [], @(x) isempty(x) || (isnumeric(x) && all(isfinite(x(:)))));
p.addParameter('FreqPpmMaskFcn',  @(ppm) ppm < 0, @(f) isa(f,'function_handle'));
p.addParameter('MetadataStep',    'AT', @(x) ischar(x) || isstring(x));
p.parse(SandTable, varargin{:});

processing = p.Results.Processing;
if isempty(processing), processing = struct(); end

inputVar   = string(p.Results.InputVar);
outputVar  = string(p.Results.OutputVar);
ampThreshIn = p.Results.AmpThresh;
maskFcn    = p.Results.FreqPpmMaskFcn;
stepName   = string(p.Results.MetadataStep);

% -----------------------------
% Validate required columns
% -----------------------------
if ~ismember(inputVar, string(SandTable.Properties.VariableNames))
    error('SandTable is missing the input column "%s".', inputVar);
end

% -----------------------------
% Initialize output column
% -----------------------------
if ~ismember(outputVar, string(SandTable.Properties.VariableNames))
    SandTable.(outputVar) = SandTable.(inputVar);
end

nSamples = height(SandTable);

% -----------------------------
% Determine per-sample thresholds
% -----------------------------
autoUsed = false;
autoReason = strings(nSamples,1);
nCandidatesPerSample = zeros(nSamples,1);

% If user provided thresholds:
% - scalar => use same for all samples
% - vector length nSamples => per-sample override
if ~isempty(ampThreshIn)
    if isscalar(ampThreshIn)
        ampThreshVec = repmat(double(ampThreshIn), nSamples, 1);
    else
        ampThreshIn = ampThreshIn(:);
        if numel(ampThreshIn) ~= nSamples
            error('If AmpThresh is a vector, it must have length nSamples=%d.', nSamples);
        end
        ampThreshVec = double(ampThreshIn);
    end
else
    % Auto per-sample from freq_ppm mask region (default ppm < 0)
    autoUsed = true;
    ampThreshVec = zeros(nSamples,1); % default fallback if no candidates

    for s = 1:nSamples
        Ti = SandTable.(inputVar){s};
        if isempty(Ti) || ~istable(Ti)
            ampThreshVec(s) = 0;
            autoReason(s) = "Empty or non-table input; defaulted threshold to 0.";
            continue;
        end

        varNames = string(Ti.Properties.VariableNames);
        if ~ismember("freq_ppm", varNames)
            error('Sample %d: table missing "freq_ppm" column.', s);
        end
        if ~ismember("amplitude", varNames)
            error('Sample %d: table missing "amplitude" column.', s);
        end

        ppm = Ti.freq_ppm;
        a   = Ti.amplitude;

        candMask = maskFcn(ppm);
        candA = a(candMask & isfinite(a));

        nCandidatesPerSample(s) = numel(candA);

        if isempty(candA)
            ampThreshVec(s) = 0;
            autoReason(s) = "No peaks in freq_ppm mask region; defaulted threshold to 0.";
        else
            ampThreshVec(s) = max(candA);
            autoReason(s) = sprintf("Threshold = max(amplitude) in freq_ppm mask region = %.6g.", ampThreshVec(s));
        end
    end
end

% -----------------------------
% Prepare metadata containers
% -----------------------------
removedTables = cell(nSamples,1);
removedCounts = zeros(nSamples,1);
origCounts    = nan(nSamples,1);
keptCounts    = nan(nSamples,1);

% -----------------------------
% Main filtering loop (PER-SAMPLE thresholds)
% -----------------------------
for s = 1:nSamples

    Ti = SandTable.(inputVar){s};

    if isempty(Ti)
        warning('Sample %d: %s empty. Skipping.', s, inputVar);
        SandTable.(outputVar){s} = Ti;
        removedTables{s} = Ti;
        origCounts(s) = 0;
        keptCounts(s) = 0;
        removedCounts(s) = 0;
        continue;
    end

    if ~istable(Ti)
        error('Sample %d: %s entry is not a table.', s, inputVar);
    end

    varNames = string(Ti.Properties.VariableNames);
    if ~ismember("amplitude", varNames)
        error('Sample %d: table missing "amplitude" column.', s);
    end

    thr = ampThreshVec(s);
    a   = Ti.amplitude;

    keepMask = (a >= thr);

    Ti_kept    = Ti(keepMask, :);
    Ti_removed = Ti(~keepMask, :);

    SandTable.(outputVar){s} = Ti_kept;

    removedTables{s} = Ti_removed;
    origCounts(s)    = height(Ti);
    keptCounts(s)    = height(Ti_kept);
    removedCounts(s) = height(Ti_removed);
end

fprintf('Created %s using PER-SAMPLE amplitude thresholds (Input=%s).\n', outputVar, inputVar);

% -----------------------------
% Append metadata WITHOUT overwriting previous processing
% -----------------------------
if ~isfield(processing, 'metadata_AT') || ~isstruct(processing.metadata_AT)
    processing.metadata_AT = struct();
end

processing.metadata_AT.table = struct();
processing.metadata_AT.table.Amplitude_Removed = removedTables';

info = struct();
info.timestamp       = datetime('now');
info.step            = stepName;
info.inputVar        = inputVar;
info.outputVar       = outputVar;

% Per-sample thresholds (main change)
info.ampThreshPerSample = ampThreshVec;

% Convenience summaries
info.minThresh = min(ampThreshVec);
info.maxThresh = max(ampThreshVec);
info.medianThresh = median(ampThreshVec);

info.autoUsed        = autoUsed;
info.freqPpmMaskDesc = 'User-provided mask function handle (default @(ppm) ppm < 0).';

if autoUsed
    info.autoReasonPerSample = cellstr(autoReason);
    info.autoCandidatesPerSample = nCandidatesPerSample;
    info.totalAutoCandidates = sum(nCandidatesPerSample);
else
    info.autoReasonPerSample = {};
end

info.perSample = table((1:nSamples)', ampThreshVec, origCounts, removedCounts, keptCounts, ...
    'VariableNames', {'SampleIndex','AmpThresh','OriginalPeaks','RemovedPeaks','KeptPeaks'});

info.totalRemoved    = sum(removedCounts, 'omitnan');
info.totalOriginal   = sum(origCounts, 'omitnan');
info.note            = 'Peaks removed if amplitude < AmpThresh(sample). Threshold is per-sample (derived from freq_ppm mask region by default).';

processing.metadata_AT.info = struct();
processing.metadata_AT.info.(stepName) = info;

end
