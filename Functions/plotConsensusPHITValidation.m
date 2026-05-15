function figH = plotConsensusPHITValidation(R_aligned, R_pHIT, SandTable, varargin)
% plotConsensusPHITValidation
% Creates a 3-panel diagnostic figure for evaluating consensus alignment + pHIT.
%
% Panel 1: reconstructed spectra from consensus/AlignBinSlope aligned table
% Panel 2: reconstructed spectra from pHIT aligned peaks
% Panel 3: SAND dot plot of All_Table + pHIT aligned table peaks together
%
% AlignBinSlope/consensus-aligned peaks are colored blue.
% pHIT-aligned peaks are colored purple.
%
% Example:
% figH = plotConsensusPHITValidation( ...
%     R_aligned, ...
%     R_pHIT, ...
%     syntheticUrine.alignedData, ...
%     'AllTableName','All_Table', ...
%     'PHITTableName','pHITRidgesAlign', ...
%     'PPMColName','freq_ppm', ...
%     'AlignedPPMCol','Aligned_PPM', ...
%     'UseAlignedPPM',true, ...
%     'InterFactor',1000, ...
%     'XLim',[7.5 8.3]);
%
% Required inputs:
%   R_aligned.reconstructed_ppm, R_aligned.X
%   R_pHIT.reconstructed_ppm, R_pHIT.X
%   SandTable/data struct or table containing per-sample All_Table and pHIT tables

p = inputParser;
addRequired(p, 'R_aligned');
addRequired(p, 'R_pHIT');
addRequired(p, 'SandTable');

addParameter(p, 'AllTableName', 'All_Table');
addParameter(p, 'PHITTableName', 'pHITRidgesAlign');
addParameter(p, 'PPMColName', 'freq_ppm');
addParameter(p, 'AlignedPPMCol', 'Aligned_PPM');
addParameter(p, 'UseAlignedPPM', true);
addParameter(p, 'AlignedFlagCol', 'Aligned');
addParameter(p, 'AlignmentSourceCol', 'AlignmentSource');
addParameter(p, 'ConsensusSourceLabels', {'consensus','consensus_unique','consensus_incomplete_promoted','AlignBinSlope','aligned','manual_promoted_incomplete','incomplete_promoted','UniqueBins','IncompleteCandidates'});
addParameter(p, 'PHITSourceLabels', {'pHIT','pHIT_ridge','pHITRidgesAlign'});

addParameter(p, 'InterFactor', []);
addParameter(p, 'XLim', []);
addParameter(p, 'LineWidth', 1.1);
addParameter(p, 'DotSize', 18);
addParameter(p, 'FontName', 'Times New Roman');
addParameter(p, 'FontSize', 13);
addParameter(p, 'TitleFontSize', 15);
addParameter(p, 'FigureName', 'Consensus + pHIT Alignment Validation');
addParameter(p, 'ShowLegend', true);
addParameter(p, 'SampleOrder', []);

% Colors requested by user
addParameter(p, 'AlignBinSlopeColor', [0.0000 0.3000 0.8500]); % blue
addParameter(p, 'PHITColor', [0.5200 0.1600 0.7800]);          % purple
addParameter(p, 'OtherColor', [0.35 0.35 0.35]);

parse(p, R_aligned, R_pHIT, SandTable, varargin{:});
opts = p.Results;

% Convert table to struct for easier indexing
if istable(SandTable)
    dataStruct = table2struct(SandTable);
else
    dataStruct = SandTable;
end

nSamples = numel(dataStruct);
if isempty(opts.SampleOrder)
    sampleOrder = 1:nSamples;
else
    sampleOrder = opts.SampleOrder;
end

% Make X matrices samples x points
XA = orientX(R_aligned.X, R_aligned.reconstructed_ppm);
XP = orientX(R_pHIT.X, R_pHIT.reconstructed_ppm);

% Optional vertical stacking for spectra
[XA_plot, yOffsetsA] = applyStacking(XA, opts.InterFactor);
[XP_plot, yOffsetsP] = applyStacking(XP, opts.InterFactor);

% Build dot plot data
[dotPPM, dotY, dotType] = collectDotPlotData(dataStruct, sampleOrder, opts);

figH = figure('Color','w','Name',opts.FigureName);
t = tiledlayout(3,1,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>
ax = gobjects(3,1);

%% Panel 1: Consensus/AlignBinSlope aligned reconstructed spectra
ax(1) = nexttile;
hold on;
plot(R_aligned.reconstructed_ppm, real(XA_plot)', 'Color', opts.AlignBinSlopeColor, 'LineWidth', opts.LineWidth);
set(gca,'XDir','reverse');
title('Consensus / AlignBinSlope Aligned Peaks', 'Interpreter','none', ...
    'FontName',opts.FontName, 'FontSize',opts.TitleFontSize);
ylabel('Intensity', 'FontName',opts.FontName, 'FontSize',opts.FontSize);
set(gca,'FontName',opts.FontName,'FontSize',opts.FontSize,'Box','on');

%% Panel 2: pHIT aligned reconstructed spectra
ax(2) = nexttile;
hold on;
plot(R_pHIT.reconstructed_ppm, real(XP_plot)', 'Color', opts.PHITColor, 'LineWidth', opts.LineWidth);
set(gca,'XDir','reverse');
title('pHIT Aligned Peaks', 'Interpreter','none', ...
    'FontName',opts.FontName, 'FontSize',opts.TitleFontSize);
ylabel('Intensity', 'FontName',opts.FontName, 'FontSize',opts.FontSize);
set(gca,'FontName',opts.FontName,'FontSize',opts.FontSize,'Box','on');

%% Panel 3: Dot plot from All_Table + pHIT peaks
ax(3) = nexttile;
hold on;

isConsensus = dotType == "AlignBinSlope";
isPHIT      = dotType == "pHIT";
isOther     = dotType == "Other";

if any(isOther)
    scatter(dotPPM(isOther), dotY(isOther), opts.DotSize, ...
        'MarkerFaceColor', opts.OtherColor, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.35, ...
        'DisplayName','Other / unaligned');
end

if any(isConsensus)
    scatter(dotPPM(isConsensus), dotY(isConsensus), opts.DotSize, ...
        'MarkerFaceColor', opts.AlignBinSlopeColor, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.90, ...
        'DisplayName','AlignBinSlope consensus');
end

if any(isPHIT)
    scatter(dotPPM(isPHIT), dotY(isPHIT), opts.DotSize, ...
        'MarkerFaceColor', opts.PHITColor, ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.90, ...
        'DisplayName','pHIT aligned');
end

set(gca,'XDir','reverse','YDir','reverse');
title('SAND Peak Map: Consensus + pHIT Aligned Peaks', 'Interpreter','none', ...
    'FontName',opts.FontName, 'FontSize',opts.TitleFontSize);
xlabel('Chemical Shift (ppm)', 'FontName',opts.FontName, 'FontSize',opts.FontSize);
ylabel('Sample', 'FontName',opts.FontName, 'FontSize',opts.FontSize);
set(gca,'FontName',opts.FontName,'FontSize',opts.FontSize,'Box','on');

if opts.ShowLegend
    legend('Location','bestoutside','Interpreter','none');
end

% Apply x limits if requested
if ~isempty(opts.XLim)
    for a = 1:numel(ax)
        xlim(ax(a), opts.XLim);
    end
end

linkaxes(ax,'x');

end

%% Helper functions
function Xout = orientX(X, ppm)
    if size(X,2) == numel(ppm)
        Xout = X;
    elseif size(X,1) == numel(ppm)
        Xout = X';
    else
        Xout = X;
        warning('Could not confidently orient X matrix. Expected one dimension to match ppm length.');
    end
end

function [Xplot, offsets] = applyStacking(X, interFactor)
    n = size(X,1);
    if isempty(interFactor) || interFactor == 0
        offsets = zeros(n,1);
        Xplot = X;
    else
        offsets = (0:n-1)' .* interFactor;
        Xplot = X + offsets;
    end
end

function [dotPPM, dotY, dotType] = collectDotPlotData(dataStruct, sampleOrder, opts)
    dotPPM = [];
    dotY = [];
    dotType = strings(0,1);

    consensusLabels = string(opts.ConsensusSourceLabels);
    phitLabels = string(opts.PHITSourceLabels);

    for ii = 1:numel(sampleOrder)
        s = sampleOrder(ii);
        sampleY = ii;

        % ---- All_Table peaks ----
        if isfield(dataStruct(s), opts.AllTableName)
            tbl = dataStruct(s).(opts.AllTableName);

            if istable(tbl) && height(tbl) > 0 && ismember(opts.PPMColName, tbl.Properties.VariableNames)
                ppmVals = getPlotPPM(tbl, opts);
                typeVals = classifyAllTablePeaks(tbl, opts, consensusLabels, phitLabels);

                dotPPM = [dotPPM; ppmVals(:)]; %#ok<AGROW>
                dotY = [dotY; repmat(sampleY, numel(ppmVals), 1)]; %#ok<AGROW>
                dotType = [dotType; typeVals(:)]; %#ok<AGROW>
            end
        end

        % ---- pHIT table peaks, if stored separately ----
        if isfield(dataStruct(s), opts.PHITTableName)
            tblP = dataStruct(s).(opts.PHITTableName);

            if istable(tblP) && height(tblP) > 0 && ismember(opts.PPMColName, tblP.Properties.VariableNames)
                ppmValsP = getPlotPPM(tblP, opts);

                dotPPM = [dotPPM; ppmValsP(:)]; %#ok<AGROW>
                dotY = [dotY; repmat(sampleY, numel(ppmValsP), 1)]; %#ok<AGROW>
                dotType = [dotType; repmat("pHIT", numel(ppmValsP), 1)]; %#ok<AGROW>
            end
        end
    end
end

function ppmVals = getPlotPPM(tbl, opts)
    if opts.UseAlignedPPM && ismember(opts.AlignedPPMCol, tbl.Properties.VariableNames)
        ppmVals = tbl.(opts.AlignedPPMCol);
        missingMask = isnan(ppmVals);
        if any(missingMask)
            ppmVals(missingMask) = tbl.(opts.PPMColName)(missingMask);
        end
    else
        ppmVals = tbl.(opts.PPMColName);
    end
end

function typeVals = classifyAllTablePeaks(tbl, opts, consensusLabels, phitLabels)
    n = height(tbl);
    typeVals = repmat("Other", n, 1);

    % 1) Use AlignmentSource when available.
    hasSource = ismember(opts.AlignmentSourceCol, tbl.Properties.VariableNames);
    if hasSource
        sourceVals = string(tbl.(opts.AlignmentSourceCol));

        % Exact known labels
        typeVals(ismember(sourceVals, consensusLabels)) = "AlignBinSlope";
        typeVals(ismember(sourceVals, phitLabels)) = "pHIT";

        % More forgiving text matching for labels created by different versions
        sourceLower = lower(sourceVals);
        consensusMask = contains(sourceLower, "consensus") | ...
                        contains(sourceLower, "alignbinslope") | ...
                        contains(sourceLower, "manual") | ...
                        contains(sourceLower, "incomplete");
        phitMask = contains(sourceLower, "phit");

        typeVals(consensusMask) = "AlignBinSlope";
        typeVals(phitMask) = "pHIT";
    end

    % 2) Regardless of source label, use the Aligned flag as a fallback.
    % This fixes cases where AlignmentSource exists but does not exactly match
    % the expected labels. pHIT labels remain pHIT if already identified above.
    if ismember(opts.AlignedFlagCol, tbl.Properties.VariableNames)
        alignedVals = tbl.(opts.AlignedFlagCol);

        if islogical(alignedVals)
            alignedMask = alignedVals;
        elseif isnumeric(alignedVals)
            alignedMask = alignedVals == 1;
        else
            alignedMask = ismember(lower(string(alignedVals)), ["1","true","yes","aligned"]);
        end

        typeVals(alignedMask & typeVals ~= "pHIT") = "AlignBinSlope";
    end
end
