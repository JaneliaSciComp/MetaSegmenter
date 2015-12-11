function [BWdist, ind_best] = MS_S2D_BWDistance(Ibw1, Ibw2, dist_type, verbose)    
    % Compute a distance between two binary images 
    CC1 = bwconncomp(Ibw1); 
    CC2 = bwconncomp(Ibw2);
    Nc1 = CC1.NumObjects;
    Nc2 = CC2.NumObjects;
    % Compute the vectors of component sizes:
    CSV1 = [];
    
    for k=1:numel(CC1.PixelIdxList)
        CSV1 = [CSV1 numel(CC1.PixelIdxList{k})];
    end
    CSV1 = sort(CSV1, 'descend');
    CSV2 = [];
    for k=1:numel(CC2.PixelIdxList)
        CSV2 = [CSV2 numel(CC2.PixelIdxList{k})];
    end
    CSV2 = sort(CSV2, 'descend');

    Nc = min(numel(CSV1), numel(CSV2));
    if verbose > 0
        disp(['Nc=', num2str(Nc)]);
        disp(['CSV1=' num2str(CSV1(1:Nc))]);
        disp(['CSV2=' num2str(CSV2(1:Nc))]);
    end
    BWdist = 0.;
    if dist_type == 1   
        for i=1:Nc
            BWdist = BWdist + (CSV1(i) - CSV2(i))^2;
        end
        BWdist = sqrt(BWdist);
        ind_best = Nc;
    else
        ind1 = find_num_good_regions(CSV1, Nc1);
        ind2 = find_num_good_regions(CSV2, Nc2);
        BWdist = 0;  
        ind_max = max(ind1, ind2);
        Nc = min(Nc, ind_max); 
        for i=1:Nc
            BWdist = BWdist + (CSV1(i) - CSV2(i))^2;
        end
        BWdist = sqrt(BWdist);
    end
%   disp(['Nc=' num2str(Nc) ' ind=' num2str(min(ind1,ind2)), ' BWdist=' num2str(BWdist)]);

% -----------------------------------------------------------------

function ind = find_num_good_regions(CSV, Nc)
    max_dsize = 0;
    max_ratio = 0.;
    ind = Nc;
    for i=Nc:-1:2
        dsize = CSV(i-1) - CSV(i);
        if i < Nc 
            ratio = double(dsize)/double(CSV(i));           
            if max_ratio < ratio
                ind = i;
                max_ratio = ratio;
            end
%           disp(['    i=' num2str(i) ' max_ratio=' num2str(max_ratio)]);
        end
    end



