function out = alignAndReconstructValidatedRidges( ...
    outVal, alignedData, reorderOut, varargin)
% alignAndReconstructValidatedRidges
%
% Take the *validated* ridges from validateRidgesInteractive,
% snap all included peaks in each ridge to a common ppm, and
% (optionally) reconstruct spectra using only those aligned peaks.
%
% INPUTS
%   outVal      : output struct from validateRidgesInteractive
%                 (must contain ridgePointsTable with .Accepted)
%   alignedData : table or struct containing peak table cell array
%                 (e.g., SandTable or DogStruct.FinalTable)
%   reorderOut  : struct from detectRidgesGlobal / reorderSamplesByGuidePeak
%                 (must contain .order with mapping reordered -> original)
%
% NAME–VALUE PAIRS
%   'PeakTableField' : name of field/variable with peak tables
%                      (default: 'Misaligned_Table')
%   'ppmCol'         : ppm column name inside each peak table
%                      (default: 'freq_ppm')
%   'AlignStatistic' : 'median' or 'mean' (default: 'median')
%   'NewFieldName'   : name for NEW validated+aligned peak field
%                      (default: 'ValidatedAligned_Table')
%
%   'ReconstructFcn' : [] (default) OR function handle / function name.
%                      Called as:
%                        recOut = ReconstructFcn(recStruct, ...
%                                   'PeakTableField', NewFieldName, ...
%                                   ReconstructArgs{:});
%   'ReconstructArgs': cell array of extra args passed to ReconstructFcn
%                      (default: {})
%
% OUTPUT struct "out" with:
%   .alignedDataUpdated : alignedData with ppm values modified in-place
%                         for validated peaks (AlignStatistic per ridge)
%   .validatedTables    : cell array of per-sample tables containing ONLY
%                         the validated & aligned peaks
%   .ridgeTargets       : table (RidgeID, TargetPPM)
%   .recOut             : output of ReconstructFcn (if provided), else []
%

% ---------------- Parse inputs ----------------
p = inputParser;
p.addParameter('PeakTableField', 'Misaligned_Table', @(s)ischar(s)||isstring(s));
p.addParameter('ppmCol',         'freq_ppm',      @(s)ischar(s)||isstring(s));
p.addParameter('AlignStatistic', 'median',           @(s)ischar(s)||isstring(s));
p.addParameter('NewFieldName',   'ValidatedAligned_Table', @(s)ischar(s)||isstring(s));
p.addParameter('ReconstructFcn', [],                 @(f) isempty(f) || isa(f,'function_handle') || ischar(f) || isstring(f));
p.addParameter('ReconstructArgs',{},                 @(c) iscell(c));

p.parse(varargin{:});
S = p.Results;

peakField   = char(S.PeakTableField);
ppmCol      = char(S.ppmCol);
alignStat   = lower(string(S.AlignStatistic));
newField    = char(S.NewFieldName);
reconFcn    = S.ReconstructFcn;
reconArgs   = S.ReconstructArgs;

% ---------------- Basic checks ----------------
if ~isfield(outVal, 'ridgePointsTable')
    error('outVal must contain ridgePointsTable (from validateRidgesInteractive).');
end
rp = outVal.ridgePointsTable;

if ~ismember('Accepted', rp.Properties.VariableNames)
    error('ridgePointsTable must have an "Accepted" logical column (run validateRidgesInteractive first).');
end
if ~ismember('RidgeID', rp.Properties.VariableNames) || ...
   ~ismember('SampleOrderId', rp.Properties.VariableNames) || ...
   ~ismember('PeakRow', rp.Properties.VariableNames) || ...
   ~ismember('PeakPPM', rp.Properties.VariableNames)
    error('ridgePointsTable must contain RidgeID, SampleOrderId, PeakRow, and PeakPPM columns.');
end

if ~isfield(reorderOut,'order')
    error('reorderOut must contain field "order" (reordered -> original sample index).');
end
order = reorderOut.order(:);

% Unpack alignedData peak tables
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

nSamples = numel(Tcell);
if numel(order) ~= nSamples
    warning('order length (%d) and number of peak tables (%d) differ; assuming first min(...) entries match.', numel(order), nSamples);
    nSamples = min(nSamples, numel(order));
    Tcell    = Tcell(1:nSamples);
    order    = order(1:nSamples);
end

% ---------------- Filter to accepted ridge points ----------------
accMask = rp.Accepted;
if ~any(accMask)
    warning('No accepted points in ridgePointsTable.Accepted == true. Nothing to align.');
end

rpAcc = rp(accMask, :);
if isempty(rpAcc)
    % Just return original alignedData if nothing to do
    out = struct();
    out.alignedDataUpdated = alignedData;
    out.validatedTables    = repmat({table()}, nSamples, 1);
    out.ridgeTargets       = table([],[], 'VariableNames',{'RidgeID','TargetPPM'});
    out.recOut             = [];
    return;
end

% ---------------- Compute target ppm per ridge ----------------
ridgeIDs = unique(rpAcc.RidgeID);
targetPPM = zeros(numel(ridgeIDs),1);

for i = 1:numel(ridgeIDs)
    rid = ridgeIDs(i);
    Id = (rpAcc.RidgeID == rid);
    y   = rpAcc.PeakPPM(Id);

    switch alignStat
        case "median"
            targetPPM(i) = median(y, 'omitnan');
        case "mean"
            targetPPM(i) = mean(y, 'omitnan');
        otherwise
            error('AlignStatistic "%s" not supported. Use "median" or "mean".', alignStat);
    end
end

ridgeTargets = table(ridgeIDs, targetPPM, ...
    'VariableNames', {'RidgeID','TargetPPM'});

% Map RidgeID -> targetPPM
% (we can just use lookups via ridgeIDs/targetPPM arrays)

% ---------------- Update ppm values & build validated-only tables ----------------
% Copy alignedData so we can modify in-place
alignedDataUpdated = alignedData;

% validatedTables: per-sample cell array of ONLY validated peaks
validatedTables = cell(nSamples, 1);

for i = 1:numel(ridgeIDs)
    rid = ridgeIDs(i);
    tppm = targetPPM(i);

    Id = (rpAcc.RidgeID == rid);
    sOrdVec = rpAcc.SampleOrderId(Id);
    rowVec  = rpAcc.PeakRow(Id);

    for j = 1:numel(sOrdVec)
        sOrd = sOrdVec(j);
        if sOrd < 1 || sOrd > numel(order)
            continue;
        end
        origId = order(sOrd);

        if origId < 1 || origId > numel(Tcell)
            continue;
        end

        Ti = Tcell{origId};
        if isempty(Ti) || ~ismember(ppmCol, Ti.Properties.VariableNames)
            continue;
        end

        thisRow = rowVec(j);
        if thisRow < 1 || thisRow > height(Ti)
            continue;
        end

        % --- 1) Update ppm in the main peak table copy ---
        if istable(alignedDataUpdated)
            Tmain = alignedDataUpdated.(peakField){origId};
            Tmain.(ppmCol)(thisRow) = tppm;
            alignedDataUpdated.(peakField){origId} = Tmain;
        else
            Tmain = alignedDataUpdated.(peakField){origId};
            Tmain.(ppmCol)(thisRow) = tppm;
            alignedDataUpdated.(peakField){origId} = Tmain;
        end

        % --- 2) Add this (modified) row to validated-only table ----
        TiUpdated = Tmain(thisRow, :);  % row AFTER ppm update

        if isempty(validatedTables{origId})
            validatedTables{origId} = TiUpdated;
        else
            validatedTables{origId} = [validatedTables{origId}; TiUpdated]; %#ok<AGROW>
        end
    end
end

% ---------------- Attach validatedTables into a new field ----------------
if istable(alignedDataUpdated)
    alignedDataUpdated.(newField) = validatedTables;
else
    alignedDataUpdated.(newField) = validatedTables;
end

% ---------------- Optional reconstruction ----------------
recOut = [];
if ~isempty(reconFcn)
    recStruct = alignedDataUpdated;  % structure/table with newField added

    if isa(reconFcn, 'function_handle')
        recOut = reconFcn(recStruct, 'PeakTableField', newField, reconArgs{:});
    else
        % function name in char/string form
        recOut = feval(reconFcn, recStruct, 'PeakTableField', newField, reconArgs{:});
    end
end

% ---------------- Package outputs ----------------
out = struct();
out.alignedDataUpdated = alignedDataUpdated;
out.validatedTables    = validatedTables;
out.ridgeTargets       = ridgeTargets;
out.recOut             = recOut;

end
