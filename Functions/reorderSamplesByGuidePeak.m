function out = reorderSamplesByGuidePeak(ppm, X, ppmRangeGuide, NP, varargin)
% {
%   Deanna Lanier 1.26.2026 
%   
%       reorderSamplesByGuidePeak
%       Reorder spectra based on a guide peak region using Sicong's original approach.
%
%    INPUTS
%       ppm           - 1 x N vector of chemical shifts (ppm)
%       X             - S x N matrix of spectra (rows = samples)
%       ppmRangeGuide - 1 x 2 vector [ppm_right ppm_left] defining the guide region
%       NP            - number of guide peaks
%
%   OPTIONAL INPUTS
%       'InterFactor' - vertical spacing between stacked spectra (default: 1e4)
%       'MakePlots'   - logical, whether to make stacked plots (default: true)
%
%    OUTPUT 
%       out.order          - sample order (indices into original rows of X)
%       out.X_reordered    - X reordered according to guide peak
%       out.OrderMat       - [origIndex, reorderedIndex]
%       out.locV           - peak locations (after reorder, as in Sicong code)
%       out.pks1           - peak heights (after reorder)
%       out.ppmRegion      - ppm vector for the guide region
%       out.idRegion      - [idStart idEnd] of the guide region in ppm
%       out.X_stack        - stacked reordered spectra (for plotting)
%
% NOTES
%   - Requires ReorderAlign_FindPeaks on the path (in the Edison Lab toolbox?).

    %% ---------------- Parse inputs ----------------
    p = inputParser;
    p.addParameter('InterFactor', 10000, @(v)isnumeric(v)&&isscalar(v)&&v>0);
    p.addParameter('MakePlots',   true,   @(v)islogical(v)&&isscalar(v));    p.parse(varargin{:});
    S = p.Results;

    %% ---------------- Style parameters ----------------
    fontName   = 'Times New Roman';
    tickSize   = 14;
    labelSize  = 18;
    titleSize  = 20;

    %% ---------------- Basic checks ----------------
    ppm = ppm(:)';
    [nSamples, nPoints] = size(X);

    if numel(ppm) ~= nPoints
        error('ppm length must match number of columns in X.');
    end

    ppmRangeGuide = ppmRangeGuide(:)';

    %% ---------------- Initialization ----------------
    X_reo   = X;
    sample = nSamples;
    OrderMat = (1:sample)';

    %% ---------------- Guide region ----------------
    [~, idL] = min(abs(ppm - ppmRangeGuide(1)));
    [~, idR] = min(abs(ppm - ppmRangeGuide(2)));

    id1 = min(idL, idR);
    id2 = max(idL, idR);

    ppmRegion = ppm(id1:id2);
    XRs = X_reo(:, id1:id2);

    %% ---------------- Find guide peaks ----------------
    [locV, pks1] = ReorderAlign_FindPeaks(XRs, 'guide', NP);

    XRs_withKey = [mean(locV, 2), XRs];
    [~, order] = sortrows(XRs_withKey);

    X_reordered = X_reo(order, :);
    locV_reo    = locV(order, :);
    pks1_reo    = pks1(order, :);

    OrderMat = [OrderMat, order];

    %% ---------------- Build stacked spectra ----------------
    X_stack      = X_reordered + repmat(S.InterFactor*(0:sample-1)', 1, size(X_reordered,2));
    X_stack_orig = X           + repmat(S.InterFactor*(0:sample-1)', 1, size(X,2));

    %% ---------------- Plotting ----------------
    if S.MakePlots

        fig = figure('Color','w');

        % --- Axes 1: Experimental order ---
        ax(1) = subplot(2,1,1);
        plot(ppm, X_stack_orig', 'LineWidth', 1.5);
        set(ax(1), 'XDir','reverse', ...
                   'FontName',fontName, ...
                   'FontSize',tickSize);
        title('Stacked spectra (experimental order)', ...
              'FontName',fontName,'FontSize',titleSize);
        xlabel('ppm','FontName',fontName,'FontSize',labelSize);
        ylabel('Intensity (offset)','FontName',fontName,'FontSize',labelSize);

        % --- Axes 2: Reordered ---
        ax(2) = subplot(2,1,2);
        plot(ppm, X_stack', 'LineWidth', 1.5);
        set(ax(2), 'XDir','reverse', ...
                   'FontName',fontName, ...
                   'FontSize',tickSize);
        title('Stacked spectra (reordered by guide peak)', ...
              'FontName',fontName,'FontSize',titleSize);
        xlabel('ppm','FontName',fontName,'FontSize',labelSize);
        ylabel('Intensity (offset)','FontName',fontName,'FontSize',labelSize);

        % --- Link axes ---
        linkaxes(ax, 'xy');

        % --- Guide peak zoom figure ---
        fig2 = figure('Color','w');
        ax2 = axes(fig2);

        plot(ppmRegion, X_stack(:, id1:id2)', 'LineWidth', 1.5);
        hold on;
        plot(ppmRegion(locV_reo), ...
             (0:sample-1)'*S.InterFactor + pks1_reo, '*r', 'MarkerSize', 6);
        hold off;

        set(ax2, 'XDir','reverse', ...
                 'FontName',fontName, ...
                 'FontSize',tickSize);

        title('Guide peak used for reordering', ...
              'FontName',fontName,'FontSize',titleSize);
        xlabel('ppm','FontName',fontName,'FontSize',labelSize);
        ylabel('Intensity (offset)','FontName',fontName,'FontSize',labelSize);

        linkaxes(ax2, 'xy'); %#ok<*LAXES>
    end

    %% ---------------- Outputs ----------------
    out = struct();
    out.order         = order;
    out.X_reordered   = X_reordered;
    out.OrderMat      = OrderMat;
    out.locV          = locV_reo;
    out.pks1          = pks1_reo;
    out.ppmRegion     = ppmRegion;
    out.idRegion     = [id1 id2];
    out.X_stack       = X_stack;
    out.X_stack_orig  = X_stack_orig;   
    out.ppm           = ppm;
    out.InterFactor   = S.InterFactor;   


end
