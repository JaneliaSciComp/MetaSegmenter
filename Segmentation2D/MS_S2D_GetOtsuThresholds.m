%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Detrmine the threshold intensity using N.Otsu method
% ( IEEE Trans. Syst. Man and Cybernetics,
%   vol. SMC-9, NO. 1, JANUARY 1979 )

function thr = MS_S2D_GetOtsuThresholds(Igr, max_intensity, num_thr)
    thr = 0.;
    debug = 0;
    if num_thr == 1
        threshold = 0;
        max_Igr = uint8(max(max(Igr)));
        if max_intensity > max_Igr
           max_intensity = max_Igr;
        end
%       disp(['max_intensity=' num2str(max_intensity) ' min_intenisty=' num2str(min(min(Igr)))]);

        % Determine the probabilities
        counts        = zeros(1, max_intensity + 1);
        probabilities = zeros(1, max_intensity + 1);
        for i=0:max_intensity
%           disp(['i=' num2str(i) ' sum(sum(Igr == i))=' num2str(sum(sum(Igr == i)))]);
            counts(i+1) = sum(sum(Igr == i));
        end
%       disp(['counts=' num2str(counts)]);
        probabilities = double(counts)/double(sum(counts));
%       disp(['probabilities=' num2str(probabilities)]);

        % Determine the best threshold
        best_threshold = 0;
        best_eta = 0;   
        for k=1:(numel(probabilities)-1)
            omega_0 = sum(probabilities(1    :k                   ));
            omega_1 = sum(probabilities((k+1):numel(probabilities)));
            % Get mu_0
            mu_0 = 0;
            for i=1:k
                mu_0 = mu_0 + double(i)*probabilities(i);
            end
            mu_0 = mu_0/omega_0;
            % Get mu_1
            mu_1 = 0;
            for i=(k+1):numel(probabilities)
                mu_1 = mu_1 + double(i)*probabilities(i);
            end
            mu_1 = mu_1/omega_1;
            % Get mu_T
            mu_T = 0;
            for i=1:numel(probabilities)
                mu_T = mu_T + double(i)*probabilities(i);
            end
            % Get sigma2_0
            sigma2_0 = 0;
            for i=1:k
                sigma2_0 = sigma2_0 + (double(i) - mu_0)^2*probabilities(i);
            end
            sigma2_0 = sigma2_0/omega_0;
            % Get sigma2_1
            sigma2_1 = 0;
            for i=(k+1):numel(probabilities)
                sigma2_1 = sigma2_1 + (double(i) - mu_1)^2*probabilities(i);
            end
            sigma2_1 = sigma2_1/omega_1;
            % Get sigma2_T
            sigma2_T = 0;
            for i=1:numel(probabilities)
                sigma2_T = sigma2_T + (double(i) - mu_T)^2*probabilities(i);
            end
            % Get eta
            sigma2_W = omega_0*sigma2_0        + omega_1*sigma2_1;
            sigma2_B = omega_0*(mu_0 - mu_T)^2 + omega_1*(mu_1 - mu_T)^2;
            eta = sigma2_B/sigma2_T;
%       disp(['    k=' num2str(k) ' sigma2_W=' num2str(sigma2_W) ' sigma2_B=' num2str(sigma2_B) ' eta=' num2str(eta) ' best_eta='  num2str(best_eta) ' best_threshold=' num2str(best_threshold)]); 
            if  best_eta < eta
                best_eta = eta;
                best_threshold = k; 
            end
        end
        thr = [ best_threshold ];
    elseif num_thr == 2
        threshold1 = 0;
        threshold2 = 0
        max_Igr = uint8(max(max(Igr)));
        if max_intensity > max_Igr
           max_intensity = max_Igr;
        end
%       disp(['max_intensity=' num2str(max_intensity) ' min_intenisty=' num2str(min(min(Igr)))]);

        % Determine the probabilities
        counts        = zeros(1, max_intensity + 1);
        probabilities = zeros(1, max_intensity + 1);
        for i=0:max_intensity
%           disp(['i=' num2str(i) ' sum(sum(Igr == i))=' num2str(sum(sum(Igr == i)))]);
            counts(i+1) = sum(sum(Igr == i));
        end
        disp(['counts=' num2str(counts)]);
        probabilities = double(counts)/double(sum(counts));
        disp(['probabilities=' num2str(probabilities)]);

        % Determine the best threshold
        best_threshold1 = 0;
        best_threshold2 = 0;
        best_eta = 0;  
        for k1=1:(numel(probabilities)-2)
            for k2=(k1+1):(numel(probabilities)-1)
                omega_0 = sum(probabilities(1     :k1                 ));
                omega_1 = sum(probabilities((k1+1):k2                 ));
                omega_2 = sum(probabilities((k2+1):numel(probabilities)));
            end

            % Get mu_0
            mu_0 = 0;
            for i=1:k1
                mu_0 = mu_0 + double(i)*probabilities(i);
            end
            mu_0 = mu_0/omega_0;
            % Get mu_1
            mu_1 = 0;
            for i=(k1+1):k2                       
                mu_1 = mu_1 + double(i)*probabilities(i);
            end
            mu_1 = mu_1/omega_1;
            % Get mu_2
            mu_2 = 0;
            for i=(k2+1):numel(probabilities)
                mu_2 = mu_2 + double(i)*probabilities(i);
            end
            mu_2 = mu_2/omega_2;
            % Get mu_T
            mu_T = 0;
            for i=1:numel(probabilities)
                mu_T = mu_T + double(i)*probabilities(i);
            end

            % Get sigma2_0
            sigma2_0 = 0;
            for i=1:k1
                sigma2_0 = sigma2_0 + (double(i) - mu_0)^2*probabilities(i);
            end
            sigma2_0 = sigma2_0/omega_0;
            % Get sigma2_1
            sigma2_1 = 0;
            for i=(k1+1):k2                        
                sigma2_1 = sigma2_1 + (double(i) - mu_1)^2*probabilities(i);
            end
            sigma2_1 = sigma2_1/omega_1;
            % Get sigma2_2
            sigma2_2 = 0;
            for i=(k2+1):numel(probabilities)
                sigma2_2 = sigma2_2 + (double(i) - mu_2)^2*probabilities(i);
            end
            sigma2_2 = sigma2_2/omega_2;

            % Get sigma2_T
            sigma2_T = 0;
            for i=1:numel(probabilities)
                sigma2_T = sigma2_T + (double(i) - mu_T)^2*probabilities(i);
            end

            % Get eta
            sigma2_I = omega_0*sigma2_0        + omega_1*sigma2_1        + omega_2*sigma2_2;
            sigma2_E = omega_0*(mu_0 - mu_T)^2 + omega_1*(mu_1 - mu_T)^2 + omega_2*(mu_2 - mu_T)^2;

            disp(['sigma2_I=' num2str(sigma2_I) ' sigma2_E=' num2str(sigma2_E) ...
                 ' sigma2_I+sigma2_E=' num2str(sigma2_I+sigma2_E) ' sigma2_T=' num2str(sigma2_T)]);
            % GD
            sigma2_I = max(max(omega_0*sigma2_0, omega_1*sigma2_1), omega_2*sigma2_2);
            sigma2_E = min(min(omega_0*(mu_0 - mu_T)^2, omega_1*(mu_1 - mu_T)^2), omega_2*(mu_2 - mu_T)^2);
%           eta = sigma2_B/sigma2_T;
            eta = 0.;
            if sigma2_I > 0
                eta = sigma2_E/sigma2_I;
                eta = sigma2_E;
            end
%       disp(['    k=' num2str(k) ' sigma2_W=' num2str(sigma2_W) ' sigma2_B=' num2str(sigma2_B) ' eta=' num2str(eta) ' best_eta='  num2str(best_eta) ' best_threshold=' num2str(best_threshold)]);
            if  best_eta < eta
                best_eta = eta;
                best_threshold1 = k1;
                best_threshold2 = k2;
                best_mu_0 = mu_0;
                best_mu_1 = mu_1;
                best_mu_2 = mu_2;
                best_mu_T = mu_T;
                best_sig2_0 = sigma2_0;
                best_sig2_1 = sigma2_1;
                best_sig2_2 = sigma2_2;
                best_sig2_T = sigma2_T;
                best_sig2_I = sigma2_I;
                best_sig2_E = sigma2_E;
            end
            threshold_correction = 0;
        end
        if debug > 0
            disp(['best_k1=' num2str(best_threshold1) ' best_k2=' num2str(best_threshold2)]);
            disp(['best_mu_0=' num2str(best_mu_0) ' best_mu_1=' num2str(best_mu_1) ' best_mu_2=' num2str(best_mu_2) ...
                  ' best_mu_T=' num2str(best_mu_T)]);
            disp(['best_sig2_0=' num2str(best_sig2_0) ' best_sig2_1=' num2str(best_sig2_1) ...
                  ' best_sig2_2=' num2str(best_sig2_2) ' best_sig2_T=' num2str(best_sig2_T)]);
        end
        thr = [best_threshold1 best_threshold2];
    end
%   disp(['Final threshold=' num2str(thr)]);


