function decay = decayModelSAND(compositeSpectrum, frequencyHz, decayVal)
    %%% Author: Edison Lab
        % Code written by Zarif H.
        % Idea developed by Zarif H.
        
    %%% Modeling using the Lorentzian function
    %%% Defining a Lorentzian %%%
    %      h: intensity
    %      f: frequency; random variable
    %     f0: center frequency
    % omega: 0.5 of the Full Width Half Maximum (FWHM)

    % Lorentzian = @(A, lambda, f0, f) (A .* lambda) ./ (lambda.^2+(2*pi*(f-f0)) .^ 2);
    
    if decayVal == 0
        
        decay = decayVal;
    
    else

        Lorentzian = @(h, omega, f0, f) (h .* omega) ./ (omega.^2+ (f-f0).^2);

        %Y = real(composite_spectrum);
        %X = freq_hz;
        realFT = real(compositeSpectrum);

        model = fittype(Lorentzian, 'coefficients', ...
            {'h', 'omega', 'f0'}, 'independent', 'f');

        omega0 = decayVal./(2 .* pi);
        [peakIntensity, maxIdx] = max(realFT);
        frequency0 = frequencyHz(maxIdx);


        initVals = [peakIntensity, omega0, frequency0];                                    % initial values
        %params = fit(X(:), Y(:), model, 'StartPoint', initVals);                          % parameter estimation
        params = fit(frequencyHz(:), realFT(:), model, 'StartPoint', initVals);            % parameter estimation


        %h = params.h;
        %f0 = params.f0;
        omega = params.omega;

        %%% Calculating lambda from omega %%%
        FWHM = omega .* 2;
        decay = FWHM .* pi;
        
    end
    
end