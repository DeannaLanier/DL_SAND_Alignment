function [SandTable, processing] = removePeakRegionsFromSandTable(SandTable, regionsToRemove, varargin)
% removePeakRegionsFromSandTable
% Removes peaks whose freq_ppm fall within one or more ppm regions.
%
% INPUTS
%   SandTable         : table with per-sample peak tables stored in a column (cell array)
%   regionsToRemove   : Nx2 numeric array of [low high] ppm regions
%                       Example: [4.655 5.000; 2.90 3.05]
%
% NAME-VALUE PAIRS (optional)
%   'InputVar'        : column in SandTable containing peak tables (default: 'Con')
%   'OutputVar'       : column to write filtered peak tables (default: 'Con_RR')
%   'MetadataStep'    : metadata step name (default: 'RR')
%
% OUTPUTS
%   SandTable         : updated SandTable with OutputVar column containing filtered peak tables
%   processing        : struct with metadata:
%       processing.metadata.table.Regions_Removed : removed rows per sample (cell of tables)
%       processing.metadata.info.RR              : details (regions, counts, etc.)

% -----------------------------
% Parse inputs
% -----------------------------
p = inputParser;
p.addRequired('SandTable', @(x) istable(x));
p.addRequired('regionsToRemove', @(x) isnumeric(x) && size(x,2)==2 && ~isempty(x));
p.addParameter('InputVar',   'Con',   @(x) ischar(x) || isstring(x));
p.addParameter('OutputVar',  'Con_RR',@(x) ischar(x) || isstring(x));
p.addParameter('MetadataStep','RR',               @(x) ischar(x) || isstring(x));
p.parse(SandTable, regionsToRemove, varargin{:});

inputVar    = string(p.Results.InputVar);
outputVar   = string(p.Results.OutputVar);
stepName    = string(p.Results.MetadataStep);

regions = regionsToRemove;
% Normalize each row to [low high]
regions = sort(regions, 2);

% -----------------------------
% Initialize output column
% -----------------------------
if ~ismember(outputVar, string(SandTable.Properties.VariableNames))
    % If OutputVar doesn't exist, initialize from InputVar (if present)
    if ismember(inputVar, string(SandTable.Properties.VariableNames))
        SandTable.(outputVar) = SandTable.(inputVar);
    else
        error('SandTable is missing the input column "%s".', inputVar);
    end
end

nSamples = height(SandTable);

% -----------------------------
% Prepare metadata containers
% -----------------------------
removedTables = cell(nSamples,1);
removedCounts = zeros(nSamples,1);
origCounts    = nan(nSamples,1);
keptCounts    = nan(nSamples,1);

% -----------------------------
% Main loop
% -----------------------------
for s = 1:nSamples

    Ti = SandTable.(inputVar){s};

    if isempty(Ti)
        warning('Sample %d: %s is empty. Skipping.', s, inputVar);
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
    if ~ismember("freq_ppm", varNames)
        error('Sample %d: peak table missing "freq_ppm" column.', s);
    end

    ppm = Ti.freq_ppm;

    % Build a "remove if in ANY region" mask
    removeMask = false(size(ppm));
    for r = 1:size(regions,1)
        low  = regions(r,1);
        high = regions(r,2);
        removeMask = removeMask | (ppm >= low & ppm <= high);
    end

    keepMask = ~removeMask;

    Ti_kept    = Ti(keepMask, :);
    Ti_removed = Ti(removeMask, :);

    SandTable.(outputVar){s} = Ti_kept;

    removedTables{s} = Ti_removed;
    origCounts(s)    = height(Ti);
    keptCounts(s)    = height(Ti_kept);
    removedCounts(s) = height(Ti_removed);
end

% -----------------------------
% Build processing metadata
% -----------------------------
processing = struct();
processing.metadata = struct();

% Only the rows removed at this step
processing.metadata.table = struct();
processing.metadata.table.Regions_Removed = removedTables';

% Info/details about what was done
info = struct();
info.timestamp        = datetime('now');
info.step             = stepName;
info.inputVar         = inputVar;
info.outputVar        = outputVar;
info.regions_removed  = regions;              % Nx2 [low high]
info.perSample = table((1:nSamples)', origCounts, removedCounts, keptCounts, ...
    'VariableNames', {'SampleIndex','OriginalPeaks','RemovedPeaks','KeptPeaks'});
info.totalRemoved     = sum(removedCounts, 'omitnan');
info.totalOriginal    = sum(origCounts, 'omitnan');
info.note             = 'Peaks removed if freq_ppm falls within ANY specified region (inclusive bounds).';

processing.metadata.info = struct();
processing.metadata.info.(stepName) = info;

end
