function [minMaxArray, mergedPeakCell] = sandReconstruction(peakTables, FR, nmrParams)

%%% Author : Edison Lab
    % Code is written by Zarif H. -> updated 11/7/2025 Deanna L.
    % FR binning is modified from the concept developed by Chris E.
    % Signal reconstruction is modified from Leandro P.
    
%%% Dependencies
    % decayModelSAND.m
    
%%% Input
    % peakTables: cell array of SAND tables
    % FR: structure  containing XN and ppm
        % FR.XN is the spectral matrix, ideally normalized
        % FR.ppm is the ppm matrix
    % nmrParams: structure containing B0, SW, and TD
        % B0: field strength in Mhz
        % SW: sweep width in Hz
        % TD: # of acquired time points in the FID

%%% Output
    % mergedPeakCell: reconstructed SAND tables
        % combines multiple underlying peaks modeling the same Lorentzian
        
%%% Additional Comments:
    % automatically assigns max(1, number of cores - 2) for parallel processing

    %%%%%% Procuring the spectral and ppm matrices %%%%%%
    XN = FR.XN;
    reconstructedPPMs = FR.ppm;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% Binning the FR %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Finding the min/max of peaks in the matrix
    mins = islocalmin(XN,2); % Finding the local minimum of a matrix/function (may be numerically/2nd derivative)
                             % Finding minima row-wise
    
    %%% Finding the range of the ppms for the peaks in each spectrum 

    offset = 0.1;

    for i = 1:size(XN,1)                                        % iterating through each spectrum
        lower = min(peakTables{i}.freq_ppm) - offset;           % get the lowest ppm value of the spectrum
        upper = max(peakTables{i}.freq_ppm) + offset;           % get the highest ppm value of the spectrum
        lowerIdx = matchPPMs(lower, reconstructedPPMs(i,:));    % get the index of the lowest ppm value of the spectrum 
        upperIdx = matchPPMs(upper, reconstructedPPMs(i,:));    % get the index of the highest ppm value of the spectrum 
        mins(i, [lowerIdx, upperIdx]) = 1;
    end

    %%% Peak pick each spectrum 
    minMaxArray = struct;

    % Finding the peak maxima of the spectral set
    for i = 1:size(XN,1)
        [minMaxArray(i).maxima, minMaxArray(i).maxIDs] = findpeaks(XN(i,:)); % finds the local maxima/peak and their indices
        ppm = reconstructedPPMs(i, :);                                       
        minMaxArray(i).ppmMax = ppm(minMaxArray(i).maxIDs);                  % finds the corresponding ppm values
    end

    % Finding the peak minima of the spectral set
    for i = 1:size(XN,1)                                   % iterating through each spectrum                                              

        minArrayLower = [];                                % lower bounds of the current bin
        minArrayUpper = [];                                % upper bounds of the current bin

        for j = 1:length(minMaxArray(i).maxIDs)            % iterating through each peak (maximum) of a spectrum
            maxIdx = minMaxArray(i).maxIDs(j);

            %Locate the nearest minima surrounding the peak
            min1Idx = 0;
            min2Idx =0;

            %%% The following loops identify the minimum points sorrounding a maxima
            % a. backward tracking
            for l = maxIdx:-1:1
                if mins(i,l) == 1
                    min1Idx = l;
                    break;                                 % break the backward tracking after the first hit;
                    % no need to scan any further for lower minima for that peak
                end
            end
%Updated - DL 11/07/2025
            if min1Idx ~= 0
                minArrayLower = [minArrayLower min1Idx];
            end

            % b. forward tracking
            for l = maxIdx:size(XN, 2)                     
                if mins(i,l) == 1
                    min2Idx = l;
                    break;                                 % break the forward tracking after the first hit; 
                                                           % no need to scan any further for upper minima for that peak
                end
            end

            if min2Idx ~= 0
                minArrayUpper = [minArrayUpper min2Idx];
            end

        end

        minIDs = union(minArrayLower, minArrayUpper);      % building a common minimum ppm indices vector
        %%% Populating the minMaxArray

        

        % Comments: 
        % minArrayLower and minArrayUpper are the same apart from the first and the last indices respectively
        % So, union is not necessarily needed. Suppose, the upper limit is used. Then track the first index of the lower limit.
        % Add that as the first element of the combined array, and you should be golden.

   %Updated - DL 11/07/2025

        % minMaxArray(i).minIDsLower = minArrayLower;
        % minMaxArray(i).minIDsUpper = minArrayUpper;
        % minMaxArray(i).minIDs = minIDs;    
        % ppms = reconstructedPPMs(i, :);
        % minMaxArray(i).ppmMin = ppms(minIDs);

        % Remove invalids and duplicates and enforce bounds
        minIDs = unique(minIDs);
        minIDs = minIDs(isfinite(minIDs) & minIDs >= 1 & ...
            minIDs <= size(XN,2) & minIDs == floor(minIDs));

        % Populate structure fields
        minMaxArray(i).minIDsLower = minArrayLower;
        minMaxArray(i).minIDsUpper = minArrayUpper;
        minMaxArray(i).minIDs = minIDs;

        ppms = reconstructedPPMs(i, :);
        minMaxArray(i).ppmMin = ppms(minIDs);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%% Combining peaks for an entire experiment (parallel processing) %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% Parameters for parallel processing
    if isempty(gcp('nocreate'))
        numCores = feature('numcores');
        numWorkers = max(1, numCores - 2);
        parpool('local', numWorkers);
    end


    %%% Information to process the SAND .csvs 
    mergedPeakCell = cell(1,length(peakTables));                                    % empty cell array to be poplulated with the merged signals

    %tic;
    
    %%% NMR parameters
    operating_frequency_mhz = nmrParams.B0;
    sw_hz = nmrParams.SW;
    dw = 1 / sw_hz;
    td_points = nmrParams.TD;
    fid_time = td_points / 2 * dw;          % acquisition time

    t = linspace(0, fid_time, td_points);   % FID time-axis


    parfor d = 1:length(peakTables)                                                 % iterating through all the spectral peakTablesset
        binPPMs = minMaxArray(d).ppmMin;                                            % copying the ppm values of the spectrum minima 
        rowsTable = length(minMaxArray(d).ppmMin) - 1;
        peakTable = table(zeros(rowsTable,1), zeros(rowsTable,1), zeros(rowsTable,1), zeros(rowsTable,1), ...
            'VariableNames', {'index', 'freq_ppm', 'decay_hz', 'amplitude'});

        for p = 1:(length(minMaxArray(d).ppmMin) - 1)                               % iterating through the bins in one spectrum

            boundaryR = binPPMs(p);                                                 % boundary on the right;
            boundaryL = binPPMs(p+1);                                               % boundary on the left
            binMask = boundaryR <= peakTables{d}.freq_ppm & ...
                peakTables{d}.freq_ppm < boundaryL;                                 % logical mask to extract ppms within this boundary
            mergingIndices = peakTables{d}.index;                                   % peak indices of SAND .csvs 
            mergingIndices = mergingIndices(binMask);                               % peak indices of SAND .csvs within the bin


            %%% Initialization %%%
            sumFID = 0;
            peakAmplitude = 0;
            decayArray = [];    % to store the decay values 
                                % maximum is necessary for decay modeling

            if ~isempty(mergingIndices)

                %%%%%% Computing the merged signals %%%%%%
                for m = 1:length(mergingIndices)                                                    % iterating through the bin indices from the SAND .csvs

                    %%% Getting the parameters from SAND .csvs %%%
                    freq_ppm = peakTables{d}.freq_ppm(mergingIndices(m));
                    decay_hz = peakTables{d}.decay_hz(mergingIndices(m));
                    amplitude = peakTables{d}.amplitude(mergingIndices(m));

                    %%% Computing the FID and merged Amplitude %%%
                    freq_hz = freq_ppm * operating_frequency_mhz;                                   % converting frequency from ppm to Hz
                    decay_constant = 1 / decay_hz;                                                  % calculating the decay constant
                    fid = amplitude * exp(-t / decay_constant) .* (cos(2 * pi * freq_hz * t) ...
                        + 1i * sin(2 * pi * freq_hz * t));                                          % calculating the FID for each signal

                    sumFID = sumFID + fid;                                                          % adding up the signals
                    peakAmplitude = peakAmplitude + amplitude; % adding up the amplitude
                    decayArray = [decayArray, decay_hz];                                            % storing the decay values

                end

                %%% Computing the frequency from the merged signals %%%
                composite_spectrum = fftshift(fft(sumFID));                                 % Fourier transformation
                X = real(composite_spectrum);                                               % absorption signal

                reconstructedHz = linspace (-sw_hz, sw_hz, length(composite_spectrum));     % sweep-width in Hz

                [~, peakIdx] = max(X);                                                      % finding the index of the peak
                peakPPM = reconstructedHz(peakIdx) ./ operating_frequency_mhz;              % finding the corresponding peak frequency (ppm)

                %%% Decay modeling %%%
                maxDecay = max(decayArray);                                                 % max decay for initial value

                try 
                    lambda = decayModelSAND(composite_spectrum, reconstructedHz, maxDecay); % the modeled decay value
                catch
                    lambda = maxDecay;                                                      % returns maxDecay if fitting fails                                                        
                end

                %%% Populating the table with parameters %%%
                peakTable.freq_ppm(end-p+1) = peakPPM;
                peakTable.decay_hz(end-p+1) = lambda;
                peakTable.amplitude(end-p+1) = peakAmplitude;

            else

                %%% Flagging the ppm and amplitude where there is no peak%%% 
                %%% From the SAND table in the peak-picked bin; Needed to be removed %%%

                peakTable.freq_ppm(end-p+1) = Inf; 
                peakTable.decay_hz(end-p+1) = Inf;
                peakTable.amplitude(end-p+1) = Inf;  

            end

        end

        peakTable(peakTable.freq_ppm == Inf & peakTable.decay_hz == Inf & peakTable.amplitude == Inf, :) = [];
        peakTable.index = (1:height(peakTable))';
        %mergedPeakCell{1,d} = peakTable;            % copying the new peak table to the cell array
        mergedPeakCell{d} = peakTable;               % copying the new peak table to the cell array 

    end

    %toc;

    delete(gcp('nocreate'));
end