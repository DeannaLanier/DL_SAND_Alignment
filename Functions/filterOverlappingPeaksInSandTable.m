function Sandtable = filterOverlappingPeaksInSandTable(Sandtable, varargin)
% filterOverlappingPeaksInSandTable
%
% Removes near-overlapping peaks within each sample's peak table.
% For peaks within a specified ppm tolerance, the peak with the largest
% amplitude is retained and smaller peaks are moved to a removed-peaks table.
%
% Example:
%   Dog_11.Sandtable = filterOverlappingPeaksInSandTable(Dog_11.Sandtable);
%
%   Dog_11.Sandtable = filterOverlappingPeaksInSandTable(Dog_11.Sandtable, ...
%       'PeakTableField', 'Con_RR_AF', ...
%       'NewPeakTableField', 'Con_RR_AF2', ...
%       'PPMColName', 'freq_ppm', ...
%       'AmplitudeColName', 'amplitude', ...
%       'PPMTolerance', 0.002);

    p = inputParser;

    addRequired(p, 'Sandtable', @istable);

    addParameter(p, 'PeakTableField', 'Con_RR_AF', @ischar);
    addParameter(p, 'NewPeakTableField', 'Con_RR_AF2', @ischar);
    addParameter(p, 'RemovedPeakTableField', 'Con_RR_AF2_Removed', @ischar);
    addParameter(p, 'PPMColName', 'freq_ppm', @ischar);
    addParameter(p, 'AmplitudeColName', 'amplitude', @ischar);
    addParameter(p, 'PPMTolerance', 0.002, @isnumeric);
    addParameter(p, 'Verbose', true, @islogical);

    parse(p, Sandtable, varargin{:});

    peakField    = p.Results.PeakTableField;
    newField     = p.Results.NewPeakTableField;
    removedField = p.Results.RemovedPeakTableField;
    ppmCol       = p.Results.PPMColName;
    ampCol       = p.Results.AmplitudeColName;
    ppmTol       = p.Results.PPMTolerance;
    verbose      = p.Results.Verbose;

    if ~ismember(peakField, Sandtable.Properties.VariableNames)
        error('Peak table field "%s" was not found in Sandtable.', peakField);
    end

    nSamples = height(Sandtable);

    filteredTables = cell(nSamples, 1);
    removedTables  = cell(nSamples, 1);

    totalRemoved = 0;

    for i = 1:nSamples

        T = Sandtable.(peakField){i};

        if isempty(T) || height(T) == 0
            filteredTables{i} = T;
            removedTables{i}  = table();
            continue
        end

        if ~ismember(ppmCol, T.Properties.VariableNames)
            error('Sample %d is missing ppm column "%s".', i, ppmCol);
        end

        if ~ismember(ampCol, T.Properties.VariableNames)
            error('Sample %d is missing amplitude column "%s".', i, ampCol);
        end

        % Sort by ppm so nearby peaks are adjacent
        T.SamplePeakOriginalIndex = (1:height(T))';
        T = sortrows(T, ppmCol, 'ascend');

        ppmVals = T.(ppmCol);
        ampVals = abs(T.(ampCol));

        keepMask = true(height(T), 1);
        removedReason = strings(height(T), 1);
        retainedIndex = NaN(height(T), 1);

        % Build local clusters of peaks within ppmTol
        clusterStart = 1;

        while clusterStart <= height(T)

            clusterEnd = clusterStart;

            while clusterEnd < height(T) && ...
                    abs(ppmVals(clusterEnd + 1) - ppmVals(clusterEnd)) <= ppmTol
                clusterEnd = clusterEnd + 1;
            end

            clusterId = clusterStart:clusterEnd;

            if numel(clusterId) > 1

                % Keep largest amplitude peak in the cluster
                [~, localMaxId] = max(ampVals(clusterId));
                keepId = clusterId(localMaxId);

                removeId = setdiff(clusterId, keepId);

                keepMask(removeId) = false;
                removedReason(removeId) = "near_overlap_lower_amplitude";
                retainedIndex(removeId) = T.SamplePeakOriginalIndex(keepId);

            end

            clusterStart = clusterEnd + 1;

        end

        T_filtered = T(keepMask, :);
        T_removed  = T(~keepMask, :);

        if ~isempty(T_removed)
            T_removed.RemovedReason = removedReason(~keepMask);
            T_removed.RetainedOriginalPeakIndex = retainedIndex(~keepMask);
            T_removed.PPMToleranceUsed = repmat(ppmTol, height(T_removed), 1);
        end

        % Restore original-ish order
        T_filtered = sortrows(T_filtered, 'SamplePeakOriginalIndex');
        T_removed  = sortrows(T_removed, 'SamplePeakOriginalIndex');

        filteredTables{i} = T_filtered;
        removedTables{i}  = T_removed;

        totalRemoved = totalRemoved + height(T_removed);

        if verbose
            fprintf('Sample %d: kept %d / removed %d overlapping peaks\n', ...
                i, height(T_filtered), height(T_removed));
        end

    end

    Sandtable.(newField) = filteredTables;
    Sandtable.(removedField) = removedTables;

    if verbose
        fprintf('\nFinished filtering overlapping peaks.\n');
        fprintf('Total removed peaks: %d\n', totalRemoved);
        fprintf('Filtered peak table saved as: %s\n', newField);
        fprintf('Removed peaks saved as: %s\n', removedField);
    end

end