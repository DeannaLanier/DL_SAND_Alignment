
function out = updateSandTableWithValidatedRidges(SandTable, outVal, outAligned, varargin)
% updateSandTableWithValidatedRidges
%
% Build/update master SAND table after ridge validation.
%
% This function:
%   1) pulls accepted ridge points from outVal
%   2) creates a per-sample table of accepted ridges (NOT aligned)
%   3) creates a per-sample table of remaining peaks after removing accepted ridges
%   4) copies the aligned validated tables from outAligned into the master SAND table
%   5) optionally reconstructs/plots:
%         - accepted ridges only
%         - remaining peaks only
%
% INPUTS
%   SandTable   : master SAND table
%   outVal      : output from validateRidgesInteractive
%   outAligned  : output from alignAndReconstructValidatedRidges
%
% OPTIONAL NAME-VALUE PAIRS
%   'SourcePeakField'        : source peak-table field in SandTable
%                              (default 'Peaks_Condensed_AT')
%   'AlignedFieldName'       : output field name to store aligned validated peaks
%                              (default 'ValidatedAligned_Table')
%   'AcceptedFieldName'      : output field name to store accepted ridges only
%                              not aligned (default 'ValidatedRidgesOnly')
%   'RemainingFieldName'     : output field name to store peaks after removing
%                              accepted ridges (default 'Remaining_After_pHIT')
%   'ppmCol'                 : ppm column name for reconstruction
%                              (default 'freq_ppm')
%   'MakePlots'              : true/false (default true)
%   'AcceptedInterFactor'    : interfactor for accepted-ridges plot
%                              (default [])
%   'RemainingInterFactor'   : interfactor for remaining-peaks plot
%                              (default [])
%   'FigureVisible'          : 'on' or 'off' (default 'on')
%
% OUTPUT
%   out : struct with fields
%       .SandTableUpdated
%       .AcceptedRidgesOnly
%       .RemainingPeaksOnly
%       .acceptedIDs
%       .rp_keep
%       .reconAccepted   (if MakePlots)
%       .reconRemaining  (if MakePlots)
%       .figAccepted     (if MakePlots)
%       .figRemaining    (if MakePlots)

    p = inputParser;
    p.addRequired('SandTable', @(x) istable(x) || isstruct(x));
    p.addRequired('outVal', @isstruct);
    p.addRequired('outAligned', @isstruct);
    p.addParameter('SourcePeakField', 'Peaks_Condensed_AT', @(s)ischar(s)||isstring(s));
    p.addParameter('AlignedFieldName', 'ValidatedAligned_Table', @(s)ischar(s)||isstring(s));
    p.addParameter('AcceptedFieldName', 'ValidatedRidgesOnly', @(s)ischar(s)||isstring(s));
    p.addParameter('RemainingFieldName', 'Remaining_After_pHIT', @(s)ischar(s)||isstring(s));
    p.addParameter('ppmCol', 'freq_ppm', @(s)ischar(s)||isstring(s));
    p.addParameter('MakePlots', true, @(b)islogical(b)&&isscalar(b));
    p.addParameter('AcceptedInterFactor', [], @(v) isempty(v) || (isnumeric(v)&&isscalar(v)));
    p.addParameter('RemainingInterFactor', [], @(v) isempty(v) || (isnumeric(v)&&isscalar(v)));
    p.addParameter('FigureVisible', 'on', @(s)ischar(s)||isstring(s));
    p.parse(SandTable, outVal, outAligned, varargin{:});
    S = p.Results;

    srcField       = char(S.SourcePeakField);
    alignedField   = char(S.AlignedFieldName);
    acceptedField  = char(S.AcceptedFieldName);
    remainingField = char(S.RemainingFieldName);
    ppmCol         = char(S.ppmCol);

    mustHaveOutVal = {'acceptedRidgeIDs','ridgePointsTable','order'};
    for k = 1:numel(mustHaveOutVal)
        if ~isfield(outVal, mustHaveOutVal{k})
            error('outVal must contain field "%s".', mustHaveOutVal{k});
        end
    end

    acceptedIDs = outVal.acceptedRidgeIDs;
    rp          = outVal.ridgePointsTable;
    order       = outVal.order(:);

    if isempty(rp)
        rp_keep = rp;
    else
        if ~ismember('RidgeID', rp.Properties.VariableNames)
            error('outVal.ridgePointsTable must contain column "RidgeID".');
        end
        rp_keep = rp(ismember(rp.RidgeID, acceptedIDs), :);
    end

    Torig = getPeakTableCellArray(SandTable, srcField);
    nSamples = numel(Torig);

    AcceptedRidgesOnly = cell(nSamples,1);
    RemainingPeaksOnly = cell(nSamples,1);

    for origId = 1:nSamples
        Ti = Torig{origId};

        if isempty(Ti)
            AcceptedRidgesOnly{origId} = Ti;
            RemainingPeaksOnly{origId} = Ti;
            continue;
        end

        sOrd = find(order == origId, 1, 'first');

        if isempty(sOrd)
            AcceptedRidgesOnly{origId} = Ti([], :);
            RemainingPeaksOnly{origId} = Ti;
            continue;
        end

        if isempty(rp_keep)
            acceptedRows = [];
        else
            neededCols = {'SampleOrderId','PeakRow'};
            for c = 1:numel(neededCols)
                if ~ismember(neededCols{c}, rp_keep.Properties.VariableNames)
                    error('outVal.ridgePointsTable must contain column "%s".', neededCols{c});
                end
            end
            acceptedRows = rp_keep.PeakRow(rp_keep.SampleOrderId == sOrd);
        end

        acceptedRows = unique(acceptedRows);
        acceptedRows = acceptedRows(acceptedRows >= 1 & acceptedRows <= height(Ti));

        if isempty(acceptedRows)
            AcceptedRidgesOnly{origId} = Ti([], :);
        else
            AcceptedRidgesOnly{origId} = Ti(acceptedRows, :);
        end

        keepMask = true(height(Ti),1);
        keepMask(acceptedRows) = false;
        RemainingPeaksOnly{origId} = Ti(keepMask, :);
    end

    SandTableUpdated = SandTable;
    SandTableUpdated = setFieldLike(SandTableUpdated, acceptedField, AcceptedRidgesOnly);
    SandTableUpdated = setFieldLike(SandTableUpdated, remainingField, RemainingPeaksOnly);

    alignedCell = tryExtractAlignedCell(outAligned, alignedField, nSamples);
    if ~isempty(alignedCell)
        SandTableUpdated = setFieldLike(SandTableUpdated, alignedField, alignedCell);
    else
        warning('Could not find aligned validated tables in outAligned. Updated table will not include field "%s".', alignedField);
    end

    out = struct();
    out.SandTableUpdated   = SandTableUpdated;
    out.AcceptedRidgesOnly = AcceptedRidgesOnly;
    out.RemainingPeaksOnly = RemainingPeaksOnly;
    out.acceptedIDs        = acceptedIDs;
    out.rp_keep            = rp_keep;
    out.figAccepted        = [];
    out.figRemaining       = [];
    out.reconAccepted      = struct([]);
    out.reconRemaining     = struct([]);

    if S.MakePlots
        tmpData = SandTableUpdated;
        [r_acc.reconstructedPPMs, r_acc.reconstructed_ppm, r_acc.X] = ...
            reconstructSpectraFromPeaks(tmpData, ...
            'peaksColumn', acceptedField, ...
            'ppmColName', ppmCol);

        r_acc.X_reordered = r_acc.X(order, :);

        if isempty(S.AcceptedInterFactor)
            if any(r_acc.X_reordered(:))
                accIF = 0.02 * max(abs(r_acc.X_reordered(:)));
            else
                accIF = 1;
            end
        else
            accIF = S.AcceptedInterFactor;
        end

        figAccepted = figure('Color','w', 'Visible', S.FigureVisible);
        hold on;
        for i = size(r_acc.X_reordered,1):-1:1
            plot(r_acc.reconstructed_ppm, ...
                r_acc.X_reordered(i,:) + (i-1)*accIF, ...
                'LineWidth', 1.2);
        end       
        set(gca, 'XDir', 'reverse');
        xlabel('Chemical Shift (ppm)', 'FontSize', 12);
        ylabel('Signal Intensity (offset)', 'FontSize', 12);
        title('Accepted Ridges Only (Reordered)', 'FontSize', 14);
        set(gca, 'FontSize', 12);
        grid on;
        hold off;

        tmpData2 = SandTableUpdated;
        [r_rem.reconstructedPPMs, r_rem.reconstructed_ppm, r_rem.X] = ...
            reconstructSpectraFromPeaks(tmpData2, ...
            'peaksColumn', remainingField, ...
            'ppmColName', ppmCol);

        r_rem.X_reordered = r_rem.X(order, :);

        if isempty(S.RemainingInterFactor)
            if any(r_rem.X_reordered(:))
                remIF = 0.001 * max(abs(r_rem.X_reordered(:)));
            else
                remIF = 1;
            end
        else
            remIF = S.RemainingInterFactor;
        end

        figRemaining = figure('Color','w', 'Visible', S.FigureVisible);
        hold on;
        for i = size(r_rem.X_reordered,1):-1:1
            plot(r_rem.reconstructed_ppm, ...
                r_rem.X_reordered(i,:) + (i-1)*remIF, ...
                'LineWidth', 1.2);
        end
        set(gca, 'XDir', 'reverse');
        xlabel('Chemical Shift (ppm)', 'FontSize', 12);
        ylabel('Signal Intensity (offset)', 'FontSize', 12);
        title('Remaining Peaks After Removing Accepted Ridges (Reordered)', 'FontSize', 14);
        set(gca, 'FontSize', 12);
        grid on;
        hold off;

        out.reconAccepted = r_acc;
        out.reconRemaining = r_rem;
        out.figAccepted = figAccepted;
        out.figRemaining = figRemaining;
    end
end

function Tcell = getPeakTableCellArray(SandTable, srcField)
    if istable(SandTable)
        if ~ismember(srcField, SandTable.Properties.VariableNames)
            error('SandTable does not contain field/variable "%s".', srcField);
        end
        Tcell = SandTable.(srcField);
    elseif isstruct(SandTable)
        if ~isfield(SandTable, srcField)
            error('SandTable struct does not contain field "%s".', srcField);
        end
        Tcell = SandTable.(srcField);
    else
        error('SandTable must be a table or struct.');
    end

    if ~iscell(Tcell)
        error('Field "%s" must be a cell array of per-sample peak tables.', srcField);
    end
end

function obj = setFieldLike(obj, fieldName, value)
    if istable(obj)
        obj.(fieldName) = value;
    elseif isstruct(obj)
        obj.(fieldName) = value;
    else
        error('Unsupported object type.');
    end
end

function alignedCell = tryExtractAlignedCell(outAligned, desiredField, nSamples)
    alignedCell = [];

    if isfield(outAligned, desiredField)
        candidate = outAligned.(desiredField);
        if iscell(candidate) && numel(candidate) == nSamples
            alignedCell = candidate;
            return;
        end
    end

    commonNames = {'ValidatedAligned_Table','Aligned','AlignedTables','AlignedPeakTables'};
    for k = 1:numel(commonNames)
        nm = commonNames{k};
        if isfield(outAligned, nm)
            candidate = outAligned.(nm);
            if iscell(candidate) && numel(candidate) == nSamples
                alignedCell = candidate;
                return;
            end
        end
    end

    nestedNames = fieldnames(outAligned);
    for k = 1:numel(nestedNames)
        val = outAligned.(nestedNames{k});
        try
            if istable(val) && ismember(desiredField, val.Properties.VariableNames)
                candidate = val.(desiredField);
                if iscell(candidate) && numel(candidate) == nSamples
                    alignedCell = candidate;
                    return;
                end
            elseif isstruct(val) && isfield(val, desiredField)
                candidate = val.(desiredField);
                if iscell(candidate) && numel(candidate) == nSamples
                    alignedCell = candidate;
                    return;
                end
            end
        catch
        end
    end
end
