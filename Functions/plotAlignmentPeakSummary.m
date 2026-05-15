function summaryStats = plotAlignmentPeakSummary(dataIn, varargin)
% plotAlignmentPeakSummary
% Creates summary plots comparing the number and amplitude distribution of:
%   1) AlignBinSlope/consensus aligned peaks
%   2) pHIT aligned peaks
%   3) remaining misaligned peaks
%
% The function assumes each row/sample contains peak tables such as:
%   All_Table, Aligned_Table, Misaligned_Table, pHITRidgesAlign
%
% It preserves the original data and only reads from the supplied tables.
%
% Example:
% summaryStats = plotAlignmentPeakSummary(syntheticUrine.alignedData, ...
%     'AlignedTableName','Aligned_Table', ...
%     'MisalignedTableName','Misaligned_Table', ...
%     'PHITTableName','pHITRidgesAlign', ...
%     'AmplitudeCol','amplitude', ...
%     'UseAbsoluteAmplitude',true);

p = inputParser;
addRequired(p, 'dataIn');

addParameter(p, 'AlignedTableName', 'Aligned_Table');
addParameter(p, 'MisalignedTableName', 'Misaligned_Table');
addParameter(p, 'PHITTableName', 'pHITRidgesAlign');
addParameter(p, 'AmplitudeCol', 'amplitude');
addParameter(p, 'UseAbsoluteAmplitude', true);
addParameter(p, 'FontName', 'Times New Roman');
addParameter(p, 'FontSize', 14);
addParameter(p, 'MakeLogAmplitudePlot', true);
addParameter(p, 'FigureTitle', 'Peak Alignment Summary');
addParameter(p, 'Verbose', true);

parse(p, dataIn, varargin{:});

alignedName    = p.Results.AlignedTableName;
misalignedName = p.Results.MisalignedTableName;
phitName       = p.Results.PHITTableName;
ampCol         = p.Results.AmplitudeCol;
useAbsAmp      = p.Results.UseAbsoluteAmplitude;
fontName       = p.Results.FontName;
fontSize       = p.Results.FontSize;
makeLogPlot    = p.Results.MakeLogAmplitudePlot;

% Convert table row structure to struct if needed
if istable(dataIn)
    dataStruct = table2struct(dataIn);
else
    dataStruct = dataIn;
end

nSamples = numel(dataStruct);

alignedAmp = [];
misalignedAmp = [];
phitAmp = [];

alignedCountBySample = zeros(nSamples,1);
misalignedCountBySample = zeros(nSamples,1);
phitCountBySample = zeros(nSamples,1);

for s = 1:nSamples

    % ---------- Aligned peaks ----------
    if isfield(dataStruct(s), alignedName) && istable(dataStruct(s).(alignedName))
        T = dataStruct(s).(alignedName);
        alignedCountBySample(s) = height(T);
        alignedAmp = [alignedAmp; getAmpVector(T, ampCol, useAbsAmp)]; %#ok<AGROW>
    end

    % ---------- Misaligned peaks ----------
    if isfield(dataStruct(s), misalignedName) && istable(dataStruct(s).(misalignedName))
        T = dataStruct(s).(misalignedName);
        misalignedCountBySample(s) = height(T);
        misalignedAmp = [misalignedAmp; getAmpVector(T, ampCol, useAbsAmp)]; %#ok<AGROW>
    end

    % ---------- pHIT peaks ----------
    if isfield(dataStruct(s), phitName) && istable(dataStruct(s).(phitName))
        T = dataStruct(s).(phitName);
        phitCountBySample(s) = height(T);
        phitAmp = [phitAmp; getAmpVector(T, ampCol, useAbsAmp)]; %#ok<AGROW>
    end
end

% Remove NaNs
alignedAmp = alignedAmp(~isnan(alignedAmp));
misalignedAmp = misalignedAmp(~isnan(misalignedAmp));
phitAmp = phitAmp(~isnan(phitAmp));

countLabels = categorical({'Expanding_Window','pHIT','Misaligned'});
countLabels = reordercats(countLabels, {'Expanding_Window','pHIT','Misaligned'});
peakCounts = [numel(alignedAmp), numel(phitAmp), numel(misalignedAmp)];

% ---------- Summary statistics table ----------
summaryStats = table;
summaryStats.Group = ["Expanding_Window"; "pHIT"; "Misaligned"];
summaryStats.TotalPeaks = peakCounts(:);
summaryStats.MeanAmplitude = [mean(alignedAmp,'omitnan'); mean(phitAmp,'omitnan'); mean(misalignedAmp,'omitnan')];
summaryStats.MedianAmplitude = [median(alignedAmp,'omitnan'); median(phitAmp,'omitnan'); median(misalignedAmp,'omitnan')];
summaryStats.Q1Amplitude = [prctile(alignedAmp,25); prctile(phitAmp,25); prctile(misalignedAmp,25)];
summaryStats.Q3Amplitude = [prctile(alignedAmp,75); prctile(phitAmp,75); prctile(misalignedAmp,75)];
summaryStats.IQRAmplitude = summaryStats.Q3Amplitude - summaryStats.Q1Amplitude;
summaryStats.MinAmplitude = [min(alignedAmp,[],'omitnan'); min(phitAmp,[],'omitnan'); min(misalignedAmp,[],'omitnan')];
summaryStats.MaxAmplitude = [max(alignedAmp,[],'omitnan'); max(phitAmp,[],'omitnan'); max(misalignedAmp,[],'omitnan')];
summaryStats.MeanCountPerSample = [mean(alignedCountBySample,'omitnan'); mean(phitCountBySample,'omitnan'); mean(misalignedCountBySample,'omitnan')];
summaryStats.MedianCountPerSample = [median(alignedCountBySample,'omitnan'); median(phitCountBySample,'omitnan'); median(misalignedCountBySample,'omitnan')];

% ---------- Long-form amplitude table for plotting ----------
ampVals = [alignedAmp; phitAmp; misalignedAmp];
groups = [repmat("Expanding_Window", numel(alignedAmp), 1); ...
          repmat("pHIT", numel(phitAmp), 1); ...
          repmat("Misaligned", numel(misalignedAmp), 1)];

ampPlotTable = table(categorical(groups), ampVals, ...
    'VariableNames', {'Group','Amplitude'});
ampPlotTable.Group = reordercats(ampPlotTable.Group, ...
    {'Expanding_Window','pHIT','Misaligned'});

% ---------- Plot ----------
figure('Color','w');
t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
title(t, p.Results.FigureTitle, 'FontName', fontName, 'FontSize', fontSize + 4, 'FontWeight','bold');

% Top: peak counts
nexttile;
b = bar(countLabels, peakCounts, 'FaceColor','flat');
b.CData = [0.1 0.35 0.85; 0.45 0.15 0.75; 0.55 0.55 0.55];
ylabel('Number of peaks', 'FontName', fontName, 'FontSize', fontSize);
title('Peak count by alignment category', 'FontName', fontName, 'FontSize', fontSize + 2);
set(gca, 'FontName', fontName, 'FontSize', fontSize, 'Box','off');
grid on;

for i = 1:numel(peakCounts)
    text(i, peakCounts(i), sprintf(' %d', peakCounts(i)), ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', ...
        'FontName',fontName, ...
        'FontSize',fontSize - 1);
end

% Bottom: amplitude distributions
nexttile;
boxchart(ampPlotTable.Group, ampPlotTable.Amplitude);
ylabel(getAmplitudeLabel(useAbsAmp, makeLogPlot), 'FontName', fontName, 'FontSize', fontSize);
title('Amplitude distribution by alignment category', 'FontName', fontName, 'FontSize', fontSize + 2);
set(gca, 'FontName', fontName, 'FontSize', fontSize, 'Box','off');
grid on;

if makeLogPlot
    set(gca, 'YScale', 'log');
end

% Add median labels
hold on;
medVals = [summaryStats.MedianAmplitude(1), summaryStats.MedianAmplitude(2), summaryStats.MedianAmplitude(3)];
for i = 1:numel(medVals)
    if ~isnan(medVals(i)) && medVals(i) > 0
        text(i, medVals(i), sprintf(' median %.3g', medVals(i)), ...
            'HorizontalAlignment','left', ...
            'VerticalAlignment','bottom', ...
            'FontName',fontName, ...
            'FontSize',fontSize - 2);
    end
end


if p.Results.Verbose
    disp(summaryStats);
end

end

function amp = getAmpVector(T, ampCol, useAbsAmp)
    if isempty(T) || ~istable(T) || ~ismember(ampCol, T.Properties.VariableNames)
        amp = [];
        return;
    end

    amp = T.(ampCol);
    amp = amp(:);

    if useAbsAmp
        amp = abs(amp);
    end
end

function labelText = getAmplitudeLabel(useAbsAmp, makeLogPlot)
    if useAbsAmp
        labelText = 'Absolute amplitude';
    else
        labelText = 'Amplitude';
    end

    if makeLogPlot
        labelText = [labelText ' (log scale)'];
    end
end
