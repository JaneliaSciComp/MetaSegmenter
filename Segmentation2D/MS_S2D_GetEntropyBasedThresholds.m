%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%
% Compute the threshold using
% 1) Shannon entropy (a == 0, q == 1) 
% 2) Renyi's entropy (a != 1, q == 1) (Sahoo et al, Pattern Recognition, Vol. 30, No. 1, pp. 71-84, 1997)
% 3) Tsallis entropy (a == 1; q != 1) (de Albuquerque et al, Pattern Recognition Letters 25 (2004) 1059.1065)

function thr_best = MS_S2D_GetEntropyBasedThresholds(counts, a, q)
    probabilities = double(counts)/double(sum(counts));
    thr_best = 0;
    n = numel(probabilities);
    H_best = 0;
    for k=1:(numel(probabilities)-1)
        H  = 0.;
        Pb = sum(probabilities(   1 :k));
        Pw = sum(probabilities((k+1):n));
        if a == 1 && q == 1
            % Shannon entropy
            H  = -Pb*log(Pb) - Pw*log(Pw);
        elseif a ~= 1 && q == 1
            % Renyi's entropy
            Hb = log(sum((probabilities(   1 :k)/Pb)^a))/(1 - a);
            Hw = log(sum((probabilities((k+1):n)/Pb)^a))/(1 - a);
            H = Hb + Hw; 
        elseif a == 1 && q ~= 1
            % Tsallis entropy
            Sb = (1 - sum((probabilities(   1 :k)/Pb)^q))/(q - 1);
            Sw = (1 - sum((probabilities((k+1):n)/Pw)^q))/(q - 1);
            H  = Sb + Sw + (1-q)*Sb*SW;
        else
            disp('At lest one of parameters a, q must be == 1');
        end
        if H < H_best
            thr_best = k;
            H_best   = H;
        end
    end


