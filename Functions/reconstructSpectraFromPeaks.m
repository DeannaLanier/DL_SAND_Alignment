function [reconstructedPPMs, reconstructed_ppm, X] = reconstructSpectraFromPeaks(data, varargin)
%{
    Deanna Lanier 6.28.2025 (Updated 10.29.2025; updated 11.07.2025; Updated 4.29.2026)

    Reconstruct spectra from time-domain peak data.

    Supported inputs:
      • SAND Table: rows are samples, one column contains a peak table per row
      • SAND struct / struct array
      • Cell array of tables (one per sample)
      • Single table (one sample)
      • Numeric array (Nx3) -> [ppm, amplitude, decay_hz] (format selectable)

    Optional Parameters:
      'ppmColName'            - Preferred ppm name (default 'freq_ppm')
      'operatingFrequencyMHz' - Spectrometer frequency in MHz (default 600.13282159)
      'spectralWidthHz'       - Spectral width in Hz (default 12019.2307692308)
      'tdPoints'              - Total data points (default 65536)
      'lineBroadening'        - Line broadening in Hz (default 0.2)
      'NumericArrayFormat'    - 'ppm_amp_decay' (default) or 'amp_decay_ppm'
      'peaksColumn'           - SAND table/struct field to pull peak tables from
%}

    % ---------------- Parse inputs ----------------
    p = inputParser;
    addRequired(p, 'data');
    addParameter(p, 'ppmColName', 'freq_ppm', @(s)ischar(s)||isstring(s));
    addParameter(p, 'operatingFrequencyMHz', 600.13282159, @(v)isnumeric(v)&&isscalar(v));
    addParameter(p, 'spectralWidthHz', 12019.2307692308, @(v)isnumeric(v)&&isscalar(v));
    addParameter(p, 'tdPoints', 65536, @(v)isnumeric(v)&&isscalar(v));
    addParameter(p, 'lineBroadening', 0.2, @(v)isnumeric(v)&&isscalar(v));
    addParameter(p, 'NumericArrayFormat', 'ppm_amp_decay', @(s)ischar(s)||isstring(s));
    addParameter(p, 'peaksColumn', 'Peaks', @(s)ischar(s)||isstring(s));
    addParameter(p, 'ChunkSize', 500, @(v)isnumeric(v)&&isscalar(v)&&v>=1);
    parse(p, data, varargin{:});

    ppmColName              = char(p.Results.ppmColName);
    operating_frequency_mhz = p.Results.operatingFrequencyMHz;
    sw_hz                   = p.Results.spectralWidthHz;
    td_points               = round(p.Results.tdPoints);
    lineBroadening          = p.Results.lineBroadening;
    numericFmt              = char(p.Results.NumericArrayFormat);
    peaksColumn             = char(p.Results.peaksColumn);
    chunkSize               = round(p.Results.ChunkSize);

    % ---------------- Axes computed once ----------------
    dw                = 1 / sw_hz;
    fid_time          = (td_points / 2) * dw;
    t                 = linspace(0, fid_time, td_points);
    reconstructed_ppm = linspace(-sw_hz, sw_hz, td_points) ./ operating_frequency_mhz;

    % line-broadening
    lbDecay = pi * max(lineBroadening, 0);

    % ---------------- Convert input once ----------------
    peakTables = normalizeInputToPeakTables(data, peaksColumn, numericFmt);
    nSamples   = numel(peakTables);

    X = zeros(nSamples, td_points);

    % ----------------  reconstruction  ----------------
    for i = 1:nSamples
        [ppm, amp, decay_hz] = getPeakVectors(peakTables{i}, ppmColName, numericFmt);
        nPeaks = numel(ppm);

        if nPeaks == 0
            continue;
        end

        ppm      = ppm(:);
        amp      = amp(:);
        decay_hz = max(decay_hz(:), eps);
        freq_hz  = ppm .* operating_frequency_mhz;

        sumFID = complex(zeros(1, td_points));

        % This is the same idea as the simple loop, but it computes many
        % peak FIDs at once. Smaller chunks keep memory controlled.
        for firstPeak = 1:chunkSize:nPeaks
            id = firstPeak:min(firstPeak + chunkSize - 1, nPeaks);
            exponent = (-decay_hz(id) - lbDecay + 1i .* 2 .* pi .* freq_hz(id)) .* t;
            sumFID = sumFID + sum(amp(id) .* exp(exponent), 1);
        end

        X(i,:) = real(fftshift(fft(sumFID)));
    end

    reconstructedPPMs = repmat(reconstructed_ppm, nSamples, 1);

    % ---------------- Helper functions ----------------

    function cells = normalizeInputToPeakTables(d, peakCol, numericFormat)
        % Return a cell array where each cell is one sample's peak table or numeric array.

        if istable(d)
            varNames = d.Properties.VariableNames;

            % SAND master table: one row per sample, one column stores peak tables.
            if ismember(peakCol, varNames)
                cells = d.(peakCol);
                if ~iscell(cells)
                    cells = num2cell(cells);
                end
                cells = cells(:);
                return;
            end

            % Otherwise treat it as one already-extracted peak table.
            cells = {d};
            return;
        end

        if isstruct(d)
            if numel(d) > 1
                cells = cell(numel(d), 1);
                for k = 1:numel(d)
                    cells{k} = getStructPeakTable(d(k), peakCol);
                end
            else
                cells = {getStructPeakTable(d, peakCol)};
            end
            return;
        end

        if iscell(d)
            cells = d(:);
            for k = 1:numel(cells)
                if isstruct(cells{k})
                    cells{k} = getStructPeakTable(cells{k}, peakCol);
                elseif ~(istable(cells{k}) || isnumeric(cells{k}))
                    error('Unsupported cell element at index %d. Expected table, numeric array, or struct.', k);
                end
            end
            return;
        end

        if isnumeric(d)
            cells = {d};
            return;
        end

        error('Unsupported input type for data.');
    end

    function tbl = getStructPeakTable(s, peakCol)
        if ~isfield(s, peakCol)
            error('Field "%s" not found in struct. Use ''peaksColumn'' to choose the correct field.', peakCol);
        end
        tbl = s.(peakCol);
        if ~(istable(tbl) || isnumeric(tbl))
            error('Struct field "%s" must be a peak table or numeric peak array.', peakCol);
        end
    end

    function [ppm, amp, decay_hz] = getPeakVectors(tbl, ppmName, numericFormat)
        % Extract only the three vectors needed for reconstruction.
        % No alias scanning and no table copying.

        if istable(tbl)
            if ~ismember(ppmName, tbl.Properties.VariableNames)
                error('PPM column "%s" not found. Use ''ppmColName'' to choose the correct ppm/aligned ppm column.', ppmName);
            end
            if ~ismember('amplitude', tbl.Properties.VariableNames)
                error('Required column "amplitude" not found.');
            end
            if ~ismember('decay_hz', tbl.Properties.VariableNames)
                error('Required column "decay_hz" not found.');
            end

            ppm      = tbl.(ppmName);
            amp      = tbl.amplitude;
            decay_hz = tbl.decay_hz;
            return;
        end

        if isnumeric(tbl)
            if size(tbl,2) < 3
                error('Numeric peak array must have at least 3 columns.');
            end
            switch lower(numericFormat)
                case 'ppm_amp_decay'
                    ppm = tbl(:,1); amp = tbl(:,2); decay_hz = tbl(:,3);
                case 'amp_decay_ppm'
                    amp = tbl(:,1); decay_hz = tbl(:,2); ppm = tbl(:,3);
                otherwise
                    error('Unknown NumericArrayFormat "%s".', numericFormat);
            end
            return;
        end

        error('Peak data must be a table or numeric array.');
    end
end
