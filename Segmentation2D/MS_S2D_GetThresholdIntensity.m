%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: compute a threshold intensity for each region
% and then a matrix of threshold intensities       
% using a smooth interpolation between regions
% (so that every pixel will get its own threshold)
%        

function [M_thr1, M_thr2, M_thr] = MS_S2D_GetThresholdIntensity(Igr, subsections, options)
    
    % Compute threshold intrensities for different regions
    nx = int32(options.nx);
    ny = int32(options.ny);
    dx = int32(options.dx);
    dy = int32(options.dy);
    
    thresholds1 = zeros(ny, nx);
    thresholds2 = zeros(ny, nx);
    thresholds  = zeros(ny, nx);
    y_center    = int32(zeros(1, ny));
    x_center    = int32(zeros(1, nx));
    disp(' ');
    disp(['size(Igr)=' num2str(size(Igr)) ' num_subsections=' num2str(nx * ny) ]);
    max_Igr = max(max(Igr));
    if max_Igr > 255
        Igr = mat2gray(double(Igr)/max_Igr);
    end
    disp(['In MS_S2D_GetThresholdIntensity: max(Igr)=' num2str(max(max(Igr)))]);

    size1 = int32(size(Igr,1));
    size2 = int32(size(Igr,2));

    fid = '';
    if length(options.outThr) > 0
        fid = fopen(options.outThr, 'w');
    end
    for k=1:numel(subsections)
        iy = 1 + floor(double(k-1)/double(nx));
        jx = k - nx*(iy -1);
        y_center(iy) = int32(round(mean(subsections(k).ypixels)));
        x_center(jx) = int32(round(mean(subsections(k).xpixels)));
        [ymin ymax] = get_bounds(iy, ny, y_center(1), size1);
        [xmin xmax] = get_bounds(jx, nx, x_center(1), size2);
        Igr_reg = Igr(subsections(k).ypixels,subsections(k).xpixels);
        [thr1, thr2, thr] = get_threshold_intensity_for_one_subsection([ymin ymax], ...
                                [xmin xmax], Igr_reg, options, fid);
        thresholds1(iy, jx) = thr1;
        thresholds2(iy, jx) = thr2;
        thresholds( iy, jx) = thr;
    end
    
    if fid ~= ''
        fclose(fid);
    end

    M_thr1 = do_not_interpolate_intensity_thresholds(Igr, nx, ny, ...
                 x_center, y_center, thresholds1);
    M_thr2 = do_not_interpolate_intensity_thresholds(Igr, nx, ny, ...
                 x_center, y_center, thresholds2);
    M_thr  = do_not_interpolate_intensity_thresholds(Igr, nx, ny, ...
                 x_center, y_center, thresholds);

% -- --------------------------------------------------------------------------

function M_thr = do_not_interpolate_intensity_thresholds(Igr, nx, ny, ...
                 x_center, y_center, thresholds_reg)
    size1 = int32(size(Igr,1));
    size2 = int32(size(Igr,2));
    M_thr = zeros(size(Igr));

    for iy=1:ny
        [ymin, ymax] = get_bounds(iy, ny, y_center(1), size1);
        for jx=1:nx
            [xmin, xmax] = get_bounds(jx, nx, x_center(1), size2);
            M_thr(ymin:ymax, xmin:xmax) = thresholds_reg(int32(iy),int32(jx));
        end
    end

% -- --------------------------------------------------------------------------

function [ind_min, ind_max] = get_bounds(i, n, center, my_size)
    if i < n
        ind_min = 1 + (i -1)*center*2;
        ind_max =      i    *center*2;
    else
        ind_min = 1 + (i -1)*center*2;
        ind_max = my_size;
    end
    
% -- --------------------------------------------------------------------------

function M_thr = interpolate_intensity_thresholds_matlab(Igr, nx, ny, ...
                 dx, dy, x_center, y_center, thresholds_reg, method)
    size1 = int32(size(Igr,1));
    size2 = int32(size(Igr,2));
    M_thr = zeros(size1, size2);
    disp([' size(M_thr)=' num2str(size(M_thr))]);
    XI    = zeros(size1, size2);
    YI    = zeros(size1, size2);
    for i=1:size1
        XI(i,:) = double(1:size2);  
    end
    for j=1:size2
        YI(:,j) = double(1:size1)';
    end
    cxsize = numel(x_center);
    cysize = numel(y_center);
    for i=1:cysize
        X(i,:) = double(x_center);
        X(i,1) = 1;
        X(i,numel(x_center)) = size2;
    end
    for j=1:cxsize
        Y(:,j) = double(y_center)';
        Y(1,j) = 1;
        Y(numel(y_center),j) = size1;
    end
    M_thr = interp2(X, Y, double(thresholds_reg), XI, YI, method);

% -- --------------------------------------------------------------------------

% Compute a matrix of threshold intrensities by a smooth interpolation

function M_thr = interpolate_intensity_thresholds_ms(Igr, nx, ny, ...
                 x_center, y_center, thresholds_reg)
    size1 = int32(size(Igr,1));
    size2 = int32(size(Igr,2));
    M_thr = zeros(size(Igr));
   
    M_thr(1:int32(y_center(1))     ,1:int32(x_center(1))     ) = thresholds_reg(1 ,       1 );
    M_thr(int32(y_center(ny)):size1,1:int32(x_center(1))     ) = thresholds_reg(int32(ny),1 );
    M_thr(1:int32(y_center(1))     ,int32(x_center(nx)):size2) = thresholds_reg(1 ,       int32(nx));
    M_thr(int32(y_center(ny)):size1,int32(x_center(nx)):size2) = thresholds_reg(int32(ny),int32(nx));

    if ny > 1
        % Linear interpolation
        for i=(y_center(1)+1):(y_center(ny)-1)
            iy =  max(find(y_center < i));
            for j=1:x_center(1)       
                M_thr(i,j) = thresholds_reg(iy  ,1)  ...
                           +(thresholds_reg(iy+1,1)-thresholds_reg(iy,1))...
                           *double(i             -y_center(iy)) ...
                           /double(y_center(iy+1)-y_center(iy));
            end
            for j=x_center(nx):size2
                M_thr(i,j) = thresholds_reg(iy  ,nx)  ...
                           +(thresholds_reg(iy+1,nx)-thresholds_reg(iy,nx))...
                           *double(i             -y_center(iy)) ...
                           /double(y_center(iy+1)-y_center(iy));
            end
        end
    end
    if nx > 1
        % Linear interpolation
        for j=(x_center(1)+1):(x_center(nx)-1)
            jx = max(find(x_center < j));
            for i=1:y_center(1)  
                M_thr(i,j) = thresholds_reg(1,jx  )  ...
                           +(thresholds_reg(1,jx+1)-thresholds_reg(1,jx))...
                           *double(j             -x_center(jx)) ...
                           /double(x_center(jx+1)-x_center(jx));
            end
            for i=y_center(ny):size1
                M_thr(i,j) = thresholds_reg(ny,jx  )  ...
                           +(thresholds_reg(ny,jx+1)-thresholds_reg(ny,jx))...
                           *double(j             -x_center(jx)) ...
                           /double(x_center(jx+1)-x_center(jx));
            end
        end
    end
    if nx > 1 && ny > 1
        % Smooth nonlinear interpolation
        for i=(y_center(1)+1):(y_center(ny)-1)
            iy = max(find(y_center < i));
            for j=(x_center(1)+1):(x_center(nx)-1)
                jx = max(find(x_center < j));
                if     i==y_center(iy) && j==x_center(jx)
                    M_thr(i,j) = thresholds_reg(iy,jx);
                elseif i==y_center(iy+1) && j==x_center(jx)
                    M_thr(i,j) = thresholds_reg(iy+1,jx);
                elseif i==y_center(iy) && j==x_center(jx+1)
                    M_thr(i,j) = thresholds_reg(iy,jx+1);
                elseif i==y_center(iy+1) && j==x_center(jx+1)
                    M_thr(i,j) = thresholds_reg(iy+1,jx+1);
                else
                    % Assign a weighted mean if thresholds at adjacent centers,
                    % where the weights are the reciprocal distances
                    w = zeros(1,4);
                    w(1) = 1./sqrt(double(i-y_center(iy  ))^2+double(j-x_center(jx  ))^2);
                    w(2) = 1./sqrt(double(i-y_center(iy+1))^2+double(j-x_center(jx  ))^2);
                    w(3) = 1./sqrt(double(i-y_center(iy  ))^2+double(j-x_center(jx+1))^2);
                    w(4) = 1./sqrt(double(i-y_center(iy+1))^2+double(j-x_center(jx+1))^2);
                    thr = [thresholds_reg(iy,jx)   thresholds_reg(iy+1,jx) ...
                           thresholds_reg(iy,jx+1) thresholds_reg(iy+1,jx+1)];
                    M_thr(i,j) = dot(w,thr)/sum(w);
                end
            end
        end
    end

% -- --------------------------------------------------------------------------

function [threshold, dist1, dist2] = ...
    get_optimal_threshold_using_bisection(Igr, threshold1, threshold2, options)
    Ibw1 = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold1/255.), 0);
    Ibw2 = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold2/255.), 0);

    [dist, ind] = MS_S2D_BWDistance(Ibw1, Ibw2, options.verbose);
    dist1 = dist;
    dist2 = dist;
    disp(['dist=' num2str(dist)]);
    threshold1_best = threshold1;
    threshold2_best = threshold2;
    while (threshold1_best - threshold2_best) > 1
        threshold = round((threshold1_best + threshold2_best)/2);
        Ibw = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold/255.), 0);
        [dist1, ind] = MS_S2D_BWDistance(Ibw1, Ibw, 0);
        [dist2, ind] = MS_S2D_BWDistance(Ibw2, Ibw, 0);
        if dist1 < dist2
%           Ibw1 = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold/255.), 0);
            threshold1_best = threshold;
        else
%           Ibw2 = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold/255.), 0);
            threshold2_best = threshold;
        end
    end

    % Computing the final values
    threshold = threshold1_best;
    Ibw = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold/255.), 0);
    [dist1, ind] = MS_S2D_BWDistance(Ibw1, Ibw, 0);
    [dist2, ind] = MS_S2D_BWDistance(Ibw2, Ibw, 0);

% -- --------------------------------------------------------------------------

function [threshold] = get_optimal_threshold(Igr, threshold1, threshold2, options)
    Ibw1 = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold1/255.), 0);
    Ibw2 = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold2/255.), 0);
    Ibw_prev = Ibw1;
    [dist12, ind12] = MS_S2D_BWDistance(Ibw_prev, Ibw2, 0);
    disp(' ');
    disp(['dist12=' num2str(dist12)]);

    threshold      = threshold1;
    threshold_best = threshold1;
    disp(['... before while loop threshold1=', num2str(threshold1) ' threshold2=' num2str(threshold2) ]);
    while (threshold > threshold2)
        threshold = threshold -1;
        Ibw = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold/255.), 0);
        [dist , ind ] = MS_S2D_BWDistance(Ibw_prev, Ibw, 0);
        [dist1, ind1] = MS_S2D_BWDistance(Ibw, Ibw1, 0);
        [dist2, ind2] = MS_S2D_BWDistance(Ibw, Ibw2, 0);
%       disp(['threshold=' num2str(threshold) ' dist=' num2str(dist)]);
        threshold_best = threshold+1;
        if (dist1 > dist2)
            disp(['... breaking at dist=', num2str(dist)]);
            break;
        end
        Ibw_prev = Ibw;
    end

    % Computing the final values
    disp(' ');
    disp(['threshold_best=' num2str(threshold_best) ' dist_best=' num2str(dist) ' ind_best=' num2str(ind) ' dist2=' num2str(dist2)]);
    threshold = threshold_best;
    Ibw = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold_best/255.), 0);

% -- --------------------------------------------------------------------------

% Detrmine the threshold intensity using N.Otsu's method
% ( IEEE Trans. Syst. Man and Cybernetics,  
%   vol. SMC-9, NO. 1, JANUARY 1979 )
%

function [threshold1, threshold2, threshold] = ...
          get_threshold_intensity_for_one_subsection(ybounds, xbounds, Igr, ...
                                                     options, fid);
    disp(' ');
    disp(['In get_threshold_intensity_for_one_subsection: max(Igr)=' num2str(max(max(Igr)))]);
    if options.thr > 0
        threshold = options.thr;
    elseif options.thr2 > 0
        threshold = options.thr2;
    else
        if max(max(Igr)) <= 1
            Igr = round(double(Igr) * 255.); 
        end
        disp(['    max(Igr)=' num2str(max(max(Igr))) ' min(Igr)=' num2str(min(min(Igr)))]);
        disp(' ');
        threshold1 = 0;
   
        threshold_Otsu =  MS_S2D_GetOtsuThresholds(Igr, 255, 1);
        thresh = multithresh(Igr, 2);
        threshold1 = min(threshold_Otsu, thresh(2));
        threshold_Otsu2 = MS_S2D_GetOtsuThresholds(Igr, threshold_Otsu, 1);
        threshold2 = min(threshold_Otsu2, thresh(1));
        disp(' ');
        disp(['threshold1=', num2str(threshold1) ' threshold2=' num2str(threshold2)]);

        % Compute the distance between BW images
        Ibw1 = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold1/255.), 0);
        Ibw2 = MS_S2D_AddBoundaryPadding(im2bw(Igr/255., threshold2/255.), 0);

        [dist, ind] = MS_S2D_BWDistance(Ibw1, Ibw2, options.verbose);
        [threshold, dist1, dist2] = get_optimal_threshold_using_bisection(Igr, threshold1, threshold2, options);
%       [threshold] = get_optimal_threshold(Igr, threshold1, threshold2, options);
        disp(' ');
        disp(['optimal_threshold= '  num2str(threshold)]);
        if length(options.outThr > 0)
            fprintf(fid, 'y_min,max= %d %d x_min,max= %d %d Otsu1_thr= %d Otsu2_thr= %d optimal_thr= %d\n', ...
                    ybounds, xbounds, threshold1, threshold2, threshold);
            disp(['y_min,max= ' num2str(ybounds) ' x_min,max= ' num2str(xbounds) ...
                  ' Otsu1_thr= ' num2str(threshold1) ' Otsu2_thr= '  num2str(threshold2) ...
                  ' dist= '     num2str(dist)]);
        end
    end

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_GetThresholdIntensity(grayScaleImage [,desired_fraction_of_black ])');
    return

% -----------------------------------------------------------------------------

function cum_hist = get_cumulative_histogram(Hist, Hsum)
    cum_hist = zeros(1,numel(Hist));
    my_sum = 0;
    for i=2:numel(Hist)
        my_sum = my_sum + Hist(i-1);
        cum_hist(i) = my_sum/Hsum;
    end
    % Adjusting the last bin
    my_sum = my_sum + Hist(i);
    cum_hist(numel(Hist)) = my_sum/Hsum;

% -----------------------------------------------------------------------------

function threshold_intensity = get_threshold_intensity(cum_hist, fracBlack,...
                                                       edges, verbose)
    threshold_intensity = 0;
    num_bins = numel(cum_hist);
    i = 1;
    while i<= num_bins-1
        i = i + 1;
%       disp(['i=' num2str(i) ' edges(i)=' num2str(edges(i)) ...
%             ' cum_hist(i)=' num2str(cum_hist(i)) ' fracBlack=' num2str(fracBlack)]);
        if double(cum_hist(i)) >= double(fracBlack)
            deltaI = double(edges(i) - edges(i-1));
            deltaH = cum_hist(i) - cum_hist(i-1);
            threshold_intensity = edges(i-1) + deltaI*(fracBlack - cum_hist(i-1))/deltaH;
            break
        end
    end

% -----------------------------------------------------------------------------

function fracBlack = get_fracBlack(cum_hist, threshold_intensity, edges, verbose)
    fracBlack = cum_hist(max(find(edges < threshold_intensity)));


