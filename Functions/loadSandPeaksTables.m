function [peakStruct, peakTable, peakArray] = loadSandPeaksTables(inputDir, outputFormat)
% 
% 
%{
    Deanna Lanier 6.28.25 (Updated 4.28.25)

%
% 
%
% Usage:
% [peakStruct, peakTable, peakArray] = loadAndProcessPeaks('/path/to/data');
% [~, peakTable, peakArray] = loadAndProcessPeaks('/path/to/data', {'table', 'array'});
%
% Inputs:
%   inputDir     - Path to folder containing peak CSV files
%   outputFormat - (Optional) Cell array specifying which outputs to return:
%                  Options: 'struct', 'table', 'array'
%                  Default: {'struct', 'table', 'array'}
%
% Outputs:
%   peakStruct - Struct output of the raw peak data
%   peakTable  - Table version of the peak data with extra fields
%   peakArray  - Cell array where each element is a peak table
%}
    if nargin < 2
        outputFormat = {'struct', 'table', 'array'};
    end
    
    if ischar(outputFormat) || isstring(outputFormat)
        outputFormat = cellstr(outputFormat);
    end

    peakFiles = dir(fullfile(inputDir, '*.csv'));
    peakFiles = localNaturalSortDir(peakFiles);

    peakStruct = struct('Peaks', {}, 'Filename', {});

    for i = 1:numel(peakFiles)
        thisFile = fullfile(peakFiles(i).folder, peakFiles(i).name);
        peakStruct(i).Peaks = readtable(thisFile);
        peakStruct(i).Filename = peakFiles(i).name;
    end
    
    % Initialize outputs
    peakTable = [];
    peakArray = [];

    % If 'table' or 'array' is requested, convert struct to table
    if any(strcmpi(outputFormat, 'table')) || any(strcmpi(outputFormat, 'array'))
        peakTable = struct2table(peakStruct);

        % Add Run_ID by extracting numeric part of Filename
        peakTable.Run_ID = cellfun(@(x) str2double(regexp(x, '\d+', 'match', 'once')), peakTable.Filename);

        % Reorder columns
        peakTable = peakTable(:, [{'Run_ID'}, setdiff(peakTable.Properties.VariableNames, {'Run_ID'}, 'stable')]); 
        peakTable = movevars(peakTable, 'Filename', 'After', 'Peaks');
    end

    % If 'array' is requested, extract the Peaks field into a cell array
    if any(strcmpi(outputFormat, 'array'))
        numPeaks = height(peakTable);
        peakArray = cell(1, numPeaks);
        for i = 1:numPeaks
            peakArray{i} = peakTable.Peaks{i};
        end
    end

    % Clean up unrequested outputs
    if ~any(strcmpi(outputFormat, 'struct'))
        peakStruct = [];
    end
    if ~any(strcmpi(outputFormat, 'table'))
        peakTable = [];
    end
    if ~any(strcmpi(outputFormat, 'array'))
        peakArray = [];
    end
end

function files = localNaturalSortDir(files)


    if isempty(files)
        return;
    end

    names = {files.name};
    [~, id] = sort(localNaturalSortKeys(names));
    files = files(id);
end

function keys = localNaturalSortKeys(names)
% Numeric substrings are zero-padded so file2 sorts before file10.

    names = cellstr(names(:));
    keys = strings(size(names));

    for k = 1:numel(names)
        name = lower(names{k});
        parts = regexp(name, '\d+|\D+', 'match');

        key = strings(1, numel(parts));
        for p = 1:numel(parts)
            token = parts{p};
            if ~isempty(regexp(token, '^\d+$', 'once'))
                % Pad numeric pieces to preserve numeric ordering.
                key(p) = sprintf('%020d', str2double(token));
            else
                key(p) = string(token);
            end
        end

        keys(k) = strjoin(key, '');
    end
end

