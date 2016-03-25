%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Detrmine the threshold intensity using N.Otsu method
% ( IEEE Trans. Syst. Man and Cybernetics,
%   vol. SMC-9, NO. 1, JANUARY 1979 )

function [thr,mu_b] = MS_S2D_GetOtsuThresholds(counts, num_thr)
    thr = 0.;
    mu_b = [0 0];
    probabilities = double(counts)/double(sum(counts));
%   disp(['probabilities=' num2str(probabilities)]);

    % Determine the best threshold
    best_threshold = 0;
    best_eta = 0;   
    n = numel(probabilities);
    for k=1:(numel(probabilities)-1)
        omega_0 = sum(probabilities(1    :k ));
        omega_1 = sum(probabilities((k+1):n ));
        % Get mu_0
        mu_0     = dot(probabilities(1    :k), double([1    :k]       ))   /omega_0;
        mu_1     = dot(probabilities((k+1):n), double([(k+1):n]       ))   /omega_1;
        mu_T     = dot(probabilities(1    :n), double([1    :n]       ));
        sigma2_0 = dot(probabilities(1    :k), double([1    :k] - mu_0).^2)/omega_0;
        sigma2_1 = dot(probabilities((k+1):n), double([(k+1):n] - mu_1).^2)/omega_1; 
        sigma2_T = dot(probabilities(1    :n), double([1    :n] - mu_T).^2);
        % Get eta
        sigma2_W = omega_0*sigma2_0        + omega_1*sigma2_1;
        sigma2_B = omega_0*(mu_0 - mu_T)^2 + omega_1*(mu_1 - mu_T)^2;
        eta = sigma2_B/sigma2_T;
        if  best_eta < eta
            best_eta = eta;
            best_threshold = k; 
            mu_b     = [mu_0 mu_1];
        end
    end
    thr = [ best_threshold ];
%   disp(['Final threshold=' num2str(thr)]);


