%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: compute a threshold intensity for each region
% and then a matrix of threshold intensities       
% using a smooth interpolation between regions
% (so that every pixel will get its own threshold)
%        

function [M_thr, M_thr2] = MS_S2D_GetThresholdIntensity(Igr, fracBlack_id, subsections, options)
    % fracBlack_id == 1 for detection of neural boundaries and =2 for dark structures
    % Handle inputs
    num_edges  = 255; 

    % Compute threshold intrensities for different regions
    nx = int32(options.nx);
    ny = int32(options.ny);
    dx = int32(options.dx);
    dy = int32(options.dy);
    
    thresholds_reg = zeros(ny, nx);
    thresholds_reg2= zeros(ny, nx);
    y_center = int32(zeros(1, ny));
    x_center = int32(zeros(1, nx));
    disp(['fracBlack_id=' num2str(fracBlack_id) ' size(Igr)=' num2str(size(Igr)) ' num_subsections=' num2str(nx * ny) ]);
    for k=1:numel(subsections)
        iy = 1 + floor(double(k-1)/double(nx));
        jx = k - nx*(iy -1);
        y_center(iy) = int32(round(mean(subsections(k).ypixels)));
        x_center(jx) = int32(round(mean(subsections(k).xpixels)));
        Igr_reg = Igr(subsections(k).ypixels,subsections(k).xpixels);
        if fracBlack_id == 1 
            if length(options.membPr) == 0 
                [thr, thr2] = get_threshold_intensity3_for_one_subsection(iy, jx, Igr_reg,...  
                              fracBlack_id, num_edges, options);
            else
                thr         = get_threshold_intensity2_for_one_subsection(iy, jx, Igr_reg,...
                              fracBlack_id, num_edges, options);
                thr2 = thr;
            end
        else % fracBlack_id == 2 
            if length(options.mitoPr) == 0
                [thr, thr2] = get_threshold_intensity3_for_one_subsection(iy, jx, Igr_reg,...
                              fracBlack_id, num_edges, options);
            else
                thr2        = get_threshold_intensity2_for_one_subsection(iy, jx, Igr_reg,...
                                  fracBlack_id, num_edges, options);
                thr = thr2;
            end
        end
        thresholds_reg(iy, jx)  = thr;
        thresholds_reg2(iy, jx) = thr2;
    end

    method = 'linear'; % alternatives: 'ms' 'linear  'nearest' 'spline'  'cubic'
%   method = 'nearest';
%   method = 'spline';
    if strcmp(method, 'ms') | (nx == 1 & ny == 1)
        M_thr = interpolate_intensity_thresholds_ms(Igr, nx, ny, ...
                     x_center, y_center, thresholds_reg);
        M_thr2= interpolate_intensity_thresholds_ms(Igr, nx, ny, ...
                     x_center, y_center, thresholds_reg2);
    else
        M_thr = interpolate_intensity_thresholds_matlab(Igr, nx, ny, ...
                     dx, dy, x_center, y_center, thresholds_reg, method);
        M_thr2= interpolate_intensity_thresholds_matlab(Igr, nx, ny, ...
                     dx, dy, x_center, y_center, thresholds_reg2, method);
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

function display_histogram2(Igr, Hist, edges, A_B, mu_B, sig_B,...
                           A_W, mu_W, sig_W, my_title, options)
%   clear; clc;
    gauss_B = zeros(size(edges));
    gauss_W = zeros(size(edges));
    gauss_S = zeros(size(edges)); % sum of 2 gaussians
    for i=1:numel(edges)
        gauss_B(i) = A_B*exp(-((edges(i)-mu_B)/sig_B)^2);
        gauss_W(i) = A_W*exp(-((edges(i)-mu_W)/sig_W)^2);
        gauss_S(i) = gauss_B(i) + gauss_W(i);
    end

    disp(['numel(edges)=' num2str(numel(edges)) ' numel(Hist)=' num2str(numel(Hist))]);
    figure
%   imhist(Igr);
    bar(edges, Hist);
    hold on;
    plot(edges,gauss_B,'color','r', 'LineWidth',3);
    hold on;
    plot(edges,gauss_W,'color','r', 'LineWidth',3);
    hold on;
    plot(edges,gauss_S,'color','green', 'LineWidth',3);
    title(my_title);
    waitforbuttonpress;
    if options.closeAll
        close all;
    end

% -- --------------------------------------------------------------------------

% Detrmine the threshold intensity using N.Otsu's method
% ( IEEE Trans. Syst. Man and Cybernetics,  
%   vol. SMC-9, NO. 1, JANUARY 1979 )
%

function [threshold] = get_threshold_intensity2_for_one_subsection(iy, jx, Igr, ...
                                       fracBlack_id, num_edges, options);
    if fracBlack_id == 1 && options.thr > 0
        threshold = options.thr;
    elseif fracBlack_id == 2 && options.thr2 > 0
        threshold = options.thr2;
    else
        Igr = round(double(Igr) * 255.); 
        threshold = MS_S2D_GetOtsuThresholds(Igr, 1);
    end


% -- --------------------------------------------------------------------------

function [threshold, threshold2] = get_threshold_intensity3_for_one_subsection(i, j, Igr, ...
                                       fracBlack_id, num_edges, options);
    Igr = round(double(Igr) * 255.);
    thresholds = MS_S2D_GetOtsuThresholds(Igr, 2);

    if fracBlack_id == 1 && options.thr > 0
        thresholds(1) = options.thr;
    elseif fracBlack_id == 2 && options.thr2 > 0
        thresholds(2) = options.thr2;
    end

    threshold  = thresholds(1);
    threshold2 = thresholds(2);

    disp(['threshold=' num2str(thresholds(1)) ' threshold2=' num2str(thresholds(2))]);
    disp(['Final threshold=' num2str(threshold)]);

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_GetThresholdIntensity(grayScaleImage [,desired_fraction_of_black [,num_edges]])');
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

% -----------------------------------------------------------------------------

function edges = compute_edges(num_edges)
    edges = zeros(1, num_edges);
    step = 1/double(num_edges-1);
    edges(1) = 0.;
    edges(2) = step/2.;
    edges(num_edges  ) = 1.;
    edges(num_edges-1) = 1. - step/2.;
    for i=3:(num_edges-2)
        edges(i) = edges(i-1) + step;
    end

