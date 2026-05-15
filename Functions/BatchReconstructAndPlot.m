function results = BatchReconstructAndPlot(alignedData, varargin)
% BatchReconstructAndPlot
% Reconstruct and plot spectra for selected peak-table fields in a SAND
% struct/table using reconstructSpectraFromPeaks. Returns a struct with all
% outputs (PPM axes, matrices, and figure handles) per field.
%
% Usage example:
%   results = BatchReconstructAndPlot(alignedData, ...
%       'PPMOverrides', struct( ...
%           'Peaks_Aligned',   'Aligned_ppm', ...
%           'Peaks_Condensed', 'freq_ppm' ...
%       ), ...
%       'ShowLegend', false, ...
%       'FigureVisibility', 'on');
%
% Requires: reconstructSpectraFromPeaks(data, 'ppmColName', ..., 'peaksColumn', ...)
%
% Deanna Lanier — Nov 7, 2025
% Updated — Dec 2, 2025: Only reconstruct/plot fields explicitly requested
% (via PeakFields or PPMOverrides).

p = inputParser;
addRequired(p, 'alignedData');

% Which fields to process.
% If left empty, and PPMOverrides is non-empty, we will use the
% fieldnames of PPMOverrides instead.
addParameter(p, 'PeakFields', {}, @(c) iscellstr(c) || isstring(c));

% ppm column overrides by field name, e.g. struct('All_Table','freq_ppm')
addParameter(p, 'PPMOverrides', struct(), @isstruct);

% Plotting options
addParameter(p, 'LineWidth', 1.0, @isscalar);
addParameter(p, 'ShowLegend', false, @islogical);
addParameter(p, 'FigureVisibility', 'on', @(s) ischar(s) || isstring(s)); % 'on' | 'off'
addParameter(p, 'TitlePrefix', '', @(s) ischar(s) || isstring(s));

% If your reconstruct function supports additional args (besides ppmColName
% and peaksColumn which this wrapper controls), pass them here:
addParameter(p, 'ReconstructArgs', {}, @iscell);

parse(p, alignedData, varargin{:});
PeakFieldsParam   = p.Results.PeakFields;
PPMOverrides      = p.Results.PPMOverrides;
LineWidth         = p.Results.LineWidth;
ShowLegend        = p.Results.ShowLegend;
FigureVisibility  = char(p.Results.FigureVisibility);
TitlePrefix       = char(p.Results.TitlePrefix);
ReconstructArgs   = p.Results.ReconstructArgs;

% ---------------- Normalize input to struct array ----------------
if istable(alignedData)
    dataS = table2struct(alignedData);
else
    dataS = alignedData;
end

% ---------------- Determine which fields to process --------------
if isempty(PeakFieldsParam)
    ovNames = fieldnames(PPMOverrides);
    if ~isempty(ovNames)
        % Only reconstruct/plot the fields that appear in PPMOverrides
        PeakFields = ovNames;
    else
        % Fallback default behavior if nothing specified
        PeakFields = {'Aligned_Table','Misaligned_Table','Peaks_Condensed','All_Table'};
    end
else
    PeakFields = cellstr(PeakFieldsParam);
end

% Helper: check a field exists and is a table for at least one sample
hasAnyFieldTable = @(fName) any(arrayfun(@(s) isfield(s, fName) && ...
    istable(s.(fName)) && ~isempty(s.(fName)), dataS));

% Prepare results struct
results = struct();
results.FieldsProcessed = {};
results.FieldResults    = struct();

for k = 1:numel(PeakFields)
    fName = char(PeakFields{k});

    % Determine ppm column override for this field (if any)
    ppmCol = '';
    if isfield(PPMOverrides, fName)
        ppmCol = PPMOverrides.(fName);
    end

    % Validate presence; skip if entirely missing
    if ~hasAnyFieldTable(fName)
        warnMsg = sprintf('Field "%s" not found or empty for all samples. Skipping.', fName);
        warning(warnMsg);
        results.FieldResults.(fName).status  = 'missing';
        results.FieldResults.(fName).message = warnMsg;
        continue
    end

    % Build arguments for reconstructSpectraFromPeaks
    try
        rcArgs = [{'peaksColumn', fName}, ReconstructArgs];
        if ~isempty(ppmCol)
            rcArgs = [{'ppmColName', ppmCol}, rcArgs];
        end

        [reconstructedPPMs, reconstructed_ppm, X] = ...
            reconstructSpectraFromPeaks(dataS, rcArgs{:});

        % Make a figure
        fig = figure('Visible', FigureVisibility, 'Color', 'w', ...
                     'Name', ['Reconstruction - ' fName]);
        hold on
        nS = size(X,1);
        for i = 1:nS
            plot(reconstructed_ppm, X(i,:), ...
                 'LineWidth', LineWidth, ...
                 'DisplayName', sprintf('Sample %d', i));
        end
        set(gca, 'XDir','reverse');
        xlabel('Chemical Shift (ppm)');
        ylabel('Signal Intensity');

        ttl = strtrim(sprintf('%s%s', TitlePrefix, fName));
        if isempty(ttl), ttl = fName; end
        title(ttl, 'Interpreter','none');

        if ShowLegend
            legend('show');
        end
        set(gca, 'FontName','Helvetica');

        % Store outputs
        fr = struct();
        fr.reconstructedPPMs = reconstructedPPMs;
        fr.reconstructed_ppm = reconstructed_ppm;
        fr.X                 = X;
        fr.figureHandle      = fig;
        if ~isempty(ppmCol)
            fr.ppmColUsed = ppmCol;
        else
            fr.ppmColUsed = '(auto-detected / default in reconstructSpectraFromPeaks)';
        end
        fr.nSamples          = size(X,1);
        fr.nPoints           = size(X,2);

        results.FieldResults.(fName) = fr;
        results.FieldsProcessed{end+1} = fName;

    catch ME
        warning('Reconstruction failed for field "%s": %s', fName, ME.message);
        results.FieldResults.(fName).status  = 'error';
        results.FieldResults.(fName).message = ME.message;
    end
end

end % function

% ----------------------------- helpers ----------------------------------
function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
