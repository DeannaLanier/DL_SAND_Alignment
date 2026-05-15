function dogOut = runIterativeAlignBinSlopeConsensus(dogIn, varargin)
% runTwoPassAlignBinSlopeConsensus
%
% Purpose:
%   Runs the expanding-window consensus alignment twice:
%
%   PASS 1:
%       1. Reconstruct original Con_RR_AF spectra
%       2. Run AlignBinSlopeConsensus on Con_RR_AF
%       3. Interactively validate incomplete candidates
%       4. Build a temporary aligned/misaligned table
%
%   PASS 2:
%       5. Reconstruct the PASS 1 Misaligned_Table
%       6. Run AlignBinSlopeConsensus on PASS 1 Misaligned_Table
%       7. Interactively validate incomplete candidates
%
%   FINAL:
%       8. Build final aligned/misaligned/all tables using BOTH pass 1 and pass 2
%       9. Reconstruct final outputs
%
% Important:
%   This function is intentionally simple. It does NOT automatically loop through
%   many passes. It only runs consensus twice, with validation between passes.
%
% Example:
%   Dog_1 = runTwoPassAlignBinSlopeConsensus(Dog_1, ...
%       'DogID', 1, ...
%       'FigureVisibility', 'on');

%% Parse inputs
p = inputParser;
p.FunctionName = 'runTwoPassAlignBinSlopeConsensus';

addRequired(p, 'dogIn', @isstruct);

addParameter(p, 'DogID', [], @(x) isempty(x) || isnumeric(x) || ischar(x) || isstring(x));

addParameter(p, 'OriginalPeakTableField', 'Con_RR_AF', @(x) ischar(x) || isstring(x));
addParameter(p, 'PPMColName', 'freq_ppm', @(x) ischar(x) || isstring(x));

addParameter(p, 'AlignedTableName', 'Aligned_Table', @(x) ischar(x) || isstring(x));
addParameter(p, 'MisalignedTableName', 'Misaligned_Table', @(x) ischar(x) || isstring(x));
addParameter(p, 'AllTableName', 'All_Table', @(x) ischar(x) || isstring(x));
addParameter(p, 'AlignedPPMColName', 'Aligned_PPM', @(x) ischar(x) || isstring(x));

addParameter(p, 'RegionTypes', {'incomplete_aligned_candidate','incomplete_moderate_candidate'}, @iscell);
addParameter(p, 'WindowPad', 0.02, @(x) isnumeric(x) && isscalar(x));

addParameter(p, 'ShowLegend', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'FigureVisibility', 'on', @(x) ischar(x) || isstring(x));
addParameter(p, 'Verbose', true, @(x) islogical(x) || isnumeric(x));

% Pass 1 consensus parameters
pass1Default = struct( ...
    'InitialBinSize', 0.001, ...
    'MaxBinSize', 0.02, ...
    'BinIncrement', 0.0001, ...
    'NormLowThreshold', 0.10, ...
    'RegionTypesToKeep', {{'aligned','moderate'}}, ...
    'MinUniqueCount', 5, ...
    'MergeTolerancePPM', 0.002);

% Pass 2 consensus parameters
pass2Default = struct( ...
    'InitialBinSize', 0.001, ...
    'MaxBinSize', 0.02, ...
    'BinIncrement', 0.0001, ...
    'NormLowThreshold', 0.10, ...
    'RegionTypesToKeep', {{'aligned','moderate'}}, ...
    'MinUniqueCount', 5, ...
    'MergeTolerancePPM', 0.002);

addParameter(p, 'Pass1Params', pass1Default, @isstruct);
addParameter(p, 'Pass2Params', pass2Default, @isstruct);

parse(p, dogIn, varargin{:});
opt = p.Results;

dogOut = dogIn;

origField       = char(opt.OriginalPeakTableField);
ppmCol          = char(opt.PPMColName);
alignedField    = char(opt.AlignedTableName);
misalignedField = char(opt.MisalignedTableName);
allField        = char(opt.AllTableName);
alignedPPMCol   = char(opt.AlignedPPMColName);

if isfield(dogOut, 'Sandtable')
    sandTable = dogOut.Sandtable;
elseif isfield(dogOut, 'SandTable')
    sandTable = dogOut.SandTable;
else
    error('dogIn must contain either dogIn.Sandtable or dogIn.SandTable.');
end

if isempty(opt.DogID)
    dogLabel = 'Dog';
else
    dogLabel = sprintf('Dog_%s', string(opt.DogID));
end

fprintf('\n\n========================================\n');
fprintf('STARTING TWO-PASS ALIGNMENT: %s\n', dogLabel);
fprintf('========================================\n');

%% ------------------------------------------------------------------------
% PASS 1: reconstruct original, run consensus, validate incomplete
%% ------------------------------------------------------------------------

fprintf('\nPASS 1: Reconstructing original field: %s\n', origField);

[dogOut.reconstructedPPMs, dogOut.reconstructed_ppm, dogOut.X] = reconstructSpectraFromPeaks( ...
    sandTable, ...
    'ppmColName', ppmCol, ...
    'peaksColumn', origField);

fprintf('PASS 1: Running AlignBinSlopeConsensus on %s\n', origField);

dogOut.pass1.consensusOut = runConsensusLocal( ...
    sandTable, origField, ppmCol, opt.Pass1Params, opt.Verbose);

fprintf('\n\n========================================\n');
fprintf('PASS 1 VALIDATION: %s\n', dogLabel);
fprintf('========================================\n');

dogOut.pass1.validatedIncomplete = validateIncompleteConsensusBins( ...
    dogOut.pass1.consensusOut.IncompleteCandidates, ...
    dogOut.reconstructed_ppm, ...
    dogOut.X, ...
    'RegionTypes', opt.RegionTypes, ...
    'WindowPad', opt.WindowPad, ...
    'TitlePrefix', sprintf('%s Pass 1 Incomplete Candidate Validation', dogLabel));

fprintf('\nPASS 1 validation function returned.\n');
fprintf('If you finished validating, press Enter to build PASS 1 temporary aligned/misaligned tables.\n');
pause;

fprintf('\nPASS 1: Building temporary aligned/misaligned tables.\n');

[dogOut.pass1.alignedData, ...
 dogOut.pass1.alignmentSummary, ...
 dogOut.pass1.featureTable] = buildAlignedSandTablesFromConsensus( ...
    sandTable, ...
    dogOut.pass1.consensusOut.UniqueBins, ...
    dogOut.pass1.validatedIncomplete, ...
    'PeakTableField', origField, ...
    'PPMColName', ppmCol, ...
    'AlignedTableName', alignedField, ...
    'MisalignedTableName', misalignedField, ...
    'AllTableName', allField, ...
    'AlignStatistic', 'median', ...
    'MultiplePeakRule', 'closest', ...
    'Verbose', opt.Verbose);

fprintf('\nPASS 1: Reconstructing temporary PASS 1 results.\n');

dogOut.pass1.results = BatchReconstructAndPlot( ...
    dogOut.pass1.alignedData, ...
    'PPMOverrides', struct( ...
        alignedField,    alignedPPMCol, ...
        misalignedField, ppmCol, ...
        allField,        alignedPPMCol, ...
        origField,       ppmCol), ...
    'ShowLegend', opt.ShowLegend, ...
    'FigureVisibility', opt.FigureVisibility);

%% ------------------------------------------------------------------------
% PASS 2: run consensus on PASS 1 Misaligned_Table, validate incomplete
%% ------------------------------------------------------------------------

fprintf('\n\n========================================\n');
fprintf('PASS 2: Running consensus on PASS 1 %s\n', misalignedField);
fprintf('========================================\n');

dogOut.pass2.consensusOut = runConsensusLocal( ...
    dogOut.pass1.alignedData, misalignedField, ppmCol, opt.Pass2Params, opt.Verbose);

fprintf('\nPASS 2: Reconstructing PASS 1 misaligned spectra for validation.\n');

[dogOut.pass2.reconstructedPPMs, ...
 dogOut.pass2.reconstructed_ppm, ...
 dogOut.pass2.X] = reconstructSpectraFromPeaks( ...
    dogOut.pass1.alignedData, ...
    'ppmColName', ppmCol, ...
    'peaksColumn', misalignedField);

fprintf('\n\n========================================\n');
fprintf('PASS 2 VALIDATION: %s\n', dogLabel);
fprintf('========================================\n');

dogOut.pass2.validatedIncomplete = validateIncompleteConsensusBins( ...
    dogOut.pass2.consensusOut.IncompleteCandidates, ...
    dogOut.pass2.reconstructed_ppm, ...
    dogOut.pass2.X, ...
    'RegionTypes', opt.RegionTypes, ...
    'WindowPad', opt.WindowPad, ...
    'TitlePrefix', sprintf('%s Pass 2 Misaligned Candidate Validation', dogLabel));

fprintf('\nPASS 2 validation function returned.\n');
fprintf('If you finished validating, press Enter to build final aligned/misaligned tables.\n');
pause;

%% ------------------------------------------------------------------------
% FINAL: combine pass 1 + pass 2 consensus/validated candidates, then build once
%% ------------------------------------------------------------------------

fprintf('\n\n========================================\n');
fprintf('FINAL BUILD: %s\n', dogLabel);
fprintf('========================================\n');

% Build pass 2 aligned data from the pass 1 misaligned table only
[dogOut.pass2.alignedData, ...
 dogOut.pass2.alignmentSummary, ...
 dogOut.pass2.featureTable] = buildAlignedSandTablesFromConsensus( ...
    dogOut.pass1.alignedData, ...
    dogOut.pass2.consensusOut.UniqueBins, ...
    dogOut.pass2.validatedIncomplete, ...
    'PeakTableField', misalignedField, ...
    'PPMColName', ppmCol, ...
    'AlignedTableName', alignedField, ...
    'MisalignedTableName', misalignedField, ...
    'AllTableName', allField, ...
    'AlignStatistic', 'median', ...
    'MultiplePeakRule', 'closest', ...
    'Verbose', opt.Verbose);

% Merge pass 2 aligned peaks into pass 1 cumulative output
dogOut.alignedData = mergeTwoPassAlignedData( ...
    dogOut.pass1.alignedData, ...
    dogOut.pass2.alignedData, ...
    alignedField, misalignedField, allField, ppmCol, alignedPPMCol);

fprintf('\nFINAL: Reconstructing final alignedData.\n');

dogOut.Final.results = BatchReconstructAndPlot( ...
    dogOut.alignedData, ...
    'PPMOverrides', struct( ...
        alignedField,    alignedPPMCol, ...
        misalignedField, ppmCol, ...
        allField,        alignedPPMCol, ...
        origField,       ppmCol), ...
    'ShowLegend', opt.ShowLegend, ...
    'FigureVisibility', opt.FigureVisibility);

fprintf('\n\n========================================\n');
fprintf('FINISHED TWO-PASS ALIGNMENT: %s\n', dogLabel);
fprintf('========================================\n');

end

%% ========================================================================
function consensusOut = runConsensusLocal(inputData, peakField, ppmCol, params, verboseFlag)

args = { ...
    'PeakTableField', peakField, ...
    'PPMColName', ppmCol, ...
    'InitialBinSize', params.InitialBinSize, ...
    'MaxBinSize', params.MaxBinSize, ...
    'BinIncrement', params.BinIncrement, ...
    'NormLowThreshold', params.NormLowThreshold, ...
    'RegionTypesToKeep', params.RegionTypesToKeep, ...
    'MinUniqueCount', params.MinUniqueCount, ...
    'MergeTolerancePPM', params.MergeTolerancePPM, ...
    'Verbose', verboseFlag};

if isfield(params, 'AmplitudeThreshold') && ~isempty(params.AmplitudeThreshold)
    args = [args, {'AmplitudeThreshold', params.AmplitudeThreshold}]; %#ok<AGROW>
end

consensusOut = AlignBinSlopeConsensus(inputData, args{:});

end

%% ========================================================================
function mergedData = mergeTwoPassAlignedData(pass1Data, pass2Data, alignedField, misalignedField, allField, ppmCol, alignedPPMCol)

mergedData = pass1Data;
nSamples = height(pass1Data);

for i = 1:nSamples

    pass1Aligned    = pass1Data.(alignedField){i};
    pass1Misaligned = pass1Data.(misalignedField){i};

    pass2Aligned    = pass2Data.(alignedField){i};
    pass2Misaligned = pass2Data.(misalignedField){i};

    if isempty(pass2Aligned) || height(pass2Aligned) == 0
        mergedAligned = pass1Aligned;
    else
        if ~ismember(alignedPPMCol, pass2Aligned.Properties.VariableNames)
            pass2Aligned.(alignedPPMCol) = pass2Aligned.(ppmCol);
        end

        if ~ismember('SecondPassAligned', pass2Aligned.Properties.VariableNames)
            pass2Aligned.SecondPassAligned = true(height(pass2Aligned), 1);
        end

        mergedAligned = harmonizedVertcat(pass1Aligned, pass2Aligned);
    end

    % The remaining misaligned table after pass 2 should become final misaligned.
    finalMisaligned = pass2Misaligned;

    if ~isempty(finalMisaligned) && height(finalMisaligned) > 0
        if ~ismember(alignedPPMCol, finalMisaligned.Properties.VariableNames)
            finalMisaligned.(alignedPPMCol) = finalMisaligned.(ppmCol);
        end
    end

    finalAll = harmonizedVertcat(mergedAligned, finalMisaligned);

    mergedData.(alignedField){i}    = mergedAligned;
    mergedData.(misalignedField){i} = finalMisaligned;
    mergedData.(allField){i}        = finalAll;
end

end

%% ========================================================================
function out = harmonizedVertcat(A, B)

if isempty(A) || height(A) == 0
    out = B;
    return
end

if isempty(B) || height(B) == 0
    out = A;
    return
end

varsA = A.Properties.VariableNames;
varsB = B.Properties.VariableNames;
allVars = unique([varsA varsB], 'stable');

A = addMissingVars(A, allVars);
B = addMissingVars(B, allVars);

A = A(:, allVars);
B = B(:, allVars);

out = [A; B];

end

%% ========================================================================
function T = addMissingVars(T, allVars)

for i = 1:numel(allVars)
    v = allVars{i};
    if ~ismember(v, T.Properties.VariableNames)
        T.(v) = NaN(height(T), 1);
    end
end

end
