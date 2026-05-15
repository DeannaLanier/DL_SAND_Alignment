function validatedBins = validateIncompleteConsensusBins(candidateBins, ppm, X, varargin)
% validateIncompleteConsensusBins_v2
% Interactive visual validation of incomplete/aligned candidate regions.
%
% This version supports:
%   1) Multiple RegionType values at once
%   2) Direct use of AlignBinSlopeConsensus_v2 IncompleteCandidates
%   3) Adjusted incomplete columns:
%        AdjustedBinStart_Incomplete
%        AdjustedBinEnd_Incomplete
%        AdjustedRegionType_Incomplete
%   4) Manual ppm-range editing before promoting a region
%
% Controls:
%   y = keep/promote with current range
%   e = edit ppm range, then keep/promote
%   n = reject
%   s = skip/unsure
%   q = quit
%
% Example:
% validatedIncomplete = validateIncompleteConsensusBins_v2( ...
%     consensusOut.IncompleteCandidates, ...
%     pHIT.outUpdate.reconRemaining.reconstructed_ppm, ...
%     pHIT.outUpdate.reconRemaining.X, ...
%     'RegionTypes', {'incomplete_aligned_candidate','incomplete_moderate_candidate'}, ...
%     'WindowPad', 0.02, ...
%     'TitlePrefix', 'Adjusted Incomplete Candidate Validation');

p = inputParser;

addRequired(p, 'candidateBins');
addRequired(p, 'ppm');
addRequired(p, 'X');

addParameter(p, 'RegionTypes', {'incomplete_aligned_candidate','incomplete_moderate_candidate'});
addParameter(p, 'WindowPad', 0.02);
addParameter(p, 'LineWidth', 1.2);
addParameter(p, 'FontName', 'Times New Roman');
addParameter(p, 'FontSize', 14);
addParameter(p, 'TitlePrefix', 'Incomplete Candidate Validation');
addParameter(p, 'UseAdjustedIncompleteColumns', true);
addParameter(p, 'SortByPPM', true);
addParameter(p, 'PauseAfterEdit', 0.5);

parse(p, candidateBins, ppm, X, varargin{:});

regionTypes = string(p.Results.RegionTypes);
windowPad   = p.Results.WindowPad;

% Make sure ppm is a row vector
ppm = ppm(:)';

% Make sure X is samples x ppm
if size(X,2) ~= numel(ppm) && size(X,1) == numel(ppm)
    X = X';
end

if size(X,2) ~= numel(ppm)
    error('X must be samples x ppm, or ppm x samples. Number of ppm points does not match X.');
end

reviewBins = candidateBins;

% -------------------------------------------------------------------------
% Determine which columns to use for plotting/validation
% -------------------------------------------------------------------------
useAdjusted = p.Results.UseAdjustedIncompleteColumns && ...
    all(ismember({'AdjustedBinStart_Incomplete','AdjustedBinEnd_Incomplete'}, ...
    reviewBins.Properties.VariableNames));

if useAdjusted
    binStartCol = 'AdjustedBinStart_Incomplete';
    binEndCol   = 'AdjustedBinEnd_Incomplete';
else
    if all(ismember({'ConsensusBinStart','ConsensusBinEnd'}, reviewBins.Properties.VariableNames))
        binStartCol = 'ConsensusBinStart';
        binEndCol   = 'ConsensusBinEnd';
    elseif all(ismember({'BinStart','BinEnd'}, reviewBins.Properties.VariableNames))
        binStartCol = 'BinStart';
        binEndCol   = 'BinEnd';
    else
        error('Could not find bin range columns. Expected adjusted, consensus, or BinStart/BinEnd columns.');
    end
end

% Determine region type column
if useAdjusted && ismember('AdjustedRegionType_Incomplete', reviewBins.Properties.VariableNames)
    regionCol = 'AdjustedRegionType_Incomplete';
elseif ismember('RegionType', reviewBins.Properties.VariableNames)
    regionCol = 'RegionType';
else
    regionCol = '';
end

% Filter by selected region types
if ~isempty(regionCol)
    keepMask = ismember(string(reviewBins.(regionCol)), regionTypes);
    reviewBins = reviewBins(keepMask, :);
else
    warning('No region type column found. Reviewing all rows.');
end

if isempty(reviewBins)
    warning('No candidate bins matched the requested RegionTypes.');
    validatedBins = reviewBins;
    return;
end

% Add standardized validation columns
reviewBins.ValidationDecision = repmat("unreviewed", height(reviewBins), 1);
reviewBins.PromoteToAligned   = false(height(reviewBins), 1);

reviewBins.OriginalBinStart   = reviewBins.(binStartCol);
reviewBins.OriginalBinEnd     = reviewBins.(binEndCol);
reviewBins.ValidatedBinStart  = reviewBins.(binStartCol);
reviewBins.ValidatedBinEnd    = reviewBins.(binEndCol);
reviewBins.ValidatedBinCenter = mean([reviewBins.ValidatedBinStart, reviewBins.ValidatedBinEnd], 2, 'omitnan');
reviewBins.RangeWasEdited     = false(height(reviewBins), 1);
reviewBins.ReviewNotes        = strings(height(reviewBins), 1);

% Add standardized display region type
if ~isempty(regionCol)
    reviewBins.ValidationRegionType = string(reviewBins.(regionCol));
else
    reviewBins.ValidationRegionType = repmat("candidate", height(reviewBins), 1);
end

if p.Results.SortByPPM
    reviewBins = sortrows(reviewBins, 'ValidatedBinCenter', 'ascend');
end

nBins = height(reviewBins);

fprintf('\nInteractive validation started.\n');
fprintf('Reviewing %d candidate regions.\n', nBins);
fprintf('y = keep as-is | e = edit range and keep | n = reject | s = skip | q = quit\n\n');

figure;

for i = 1:nBins

    binStart = reviewBins.ValidatedBinStart(i);
    binEnd   = reviewBins.ValidatedBinEnd(i);

    rangeLow  = min(binStart, binEnd);
    rangeHigh = max(binStart, binEnd);

    xMin = rangeLow  - windowPad;
    xMax = rangeHigh + windowPad;

    id = ppm >= xMin & ppm <= xMax;

    if ~any(id)
        warning('No ppm values found in plotting range for region %d. Skipping.', i);
        reviewBins.ValidationDecision(i) = "no_ppm_points_in_range";
        continue;
    end

    clf;
    hold on;

    plot(ppm(id), real(X(:,id))', 'k', 'LineWidth', p.Results.LineWidth);

    xline(rangeLow, '--', 'Start', 'LineWidth', 1.2);
    xline(rangeHigh, '--', 'End', 'LineWidth', 1.2);
    xline(mean([rangeLow, rangeHigh]), '-', 'Center', 'LineWidth', 1.3);

    set(gca, 'XDir', 'reverse');
    xlabel('Chemical Shift (ppm)', 'FontName', p.Results.FontName, 'FontSize', p.Results.FontSize);
    ylabel('Intensity', 'FontName', p.Results.FontName, 'FontSize', p.Results.FontSize);

    title(sprintf('%s %d of %d | %.4f-%.4f ppm | %s', ...
        p.Results.TitlePrefix, i, nBins, rangeLow, rangeHigh, string(reviewBins.ValidationRegionType(i))), ...
        'FontName', p.Results.FontName, ...
        'FontSize', p.Results.FontSize + 2, ...
        'Interpreter', 'none');

    set(gca, 'FontName', p.Results.FontName, 'FontSize', p.Results.FontSize);
    box on;
    drawnow;

    fprintf('\nRegion %d/%d\n', i, nBins);
    fprintf('Type: %s\n', string(reviewBins.ValidationRegionType(i)));
    fprintf('Current range: %.5f to %.5f ppm\n', rangeLow, rangeHigh);

    if ismember('AdjustedNormSlope_Incomplete', reviewBins.Properties.VariableNames)
        fprintf('Adjusted NormSlope: %.4f\n', reviewBins.AdjustedNormSlope_Incomplete(i));
    end

    if ismember('BestCoverageUniqueCount', reviewBins.Properties.VariableNames)
        fprintf('Best coverage unique count: %d\n', reviewBins.BestCoverageUniqueCount(i));
    end

    decision = input('Decision? y = keep, e = edit range, n = reject, s = skip, q = quit: ', 's');
    decision = lower(strtrim(decision));

    if strcmp(decision, 'y')

        reviewBins.ValidationDecision(i) = "keep_as_is";
        reviewBins.PromoteToAligned(i) = true;

    elseif strcmp(decision, 'e')

        fprintf('\nEnter corrected ppm range for this feature.\n');
        newA = input('New ppm boundary 1: ');
        newB = input('New ppm boundary 2: ');

        newStart = min(newA, newB);
        newEnd   = max(newA, newB);

        reviewBins.ValidatedBinStart(i)  = newStart;
        reviewBins.ValidatedBinEnd(i)    = newEnd;
        reviewBins.ValidatedBinCenter(i) = mean([newStart, newEnd], 'omitnan');
        reviewBins.RangeWasEdited(i)     = true;
        reviewBins.ValidationDecision(i) = "keep_edited_range";
        reviewBins.PromoteToAligned(i)   = true;

        editId = ppm >= (newStart - windowPad) & ppm <= (newEnd + windowPad);

        clf;
        hold on;

        plot(ppm(editId), real(X(:,editId))', 'k', 'LineWidth', p.Results.LineWidth);
        xline(newStart, '--', 'Edited start', 'LineWidth', 1.5);
        xline(newEnd, '--', 'Edited end', 'LineWidth', 1.5);
        xline(mean([newStart,newEnd]), '-', 'Edited center', 'LineWidth', 1.5);

        set(gca, 'XDir', 'reverse');
        xlabel('Chemical Shift (ppm)', 'FontName', p.Results.FontName, 'FontSize', p.Results.FontSize);
        ylabel('Intensity', 'FontName', p.Results.FontName, 'FontSize', p.Results.FontSize);
        title(sprintf('Edited Region %d | %.4f-%.4f ppm', i, newStart, newEnd), ...
            'FontName', p.Results.FontName, ...
            'FontSize', p.Results.FontSize + 2);
        set(gca, 'FontName', p.Results.FontName, 'FontSize', p.Results.FontSize);
        box on;
        drawnow;

        pause(p.Results.PauseAfterEdit);

    elseif strcmp(decision, 'n')

        reviewBins.ValidationDecision(i) = "reject";
        reviewBins.PromoteToAligned(i) = false;

    elseif strcmp(decision, 's')

        reviewBins.ValidationDecision(i) = "skip";
        reviewBins.PromoteToAligned(i) = false;

    elseif strcmp(decision, 'q')

        fprintf('Validation stopped early.\n');
        break;

    else

        reviewBins.ValidationDecision(i) = "unreviewed_invalid_entry";
        reviewBins.PromoteToAligned(i) = false;

    end
end

validatedBins = reviewBins;

end
