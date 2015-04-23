%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: compute a threshold intensity for each region
% and then a matrix of threshold intensities       
% using a smooth interpolation between regions
% (so that every pixel will get its own threshold)
%        

function M_thr = MS_S2D_GetThresholdIntensity(Igr, fracBlack_id, subsections, options)
    % fracBlack_id == 1 for detection of neural boundaries and =2 for dark structures
    % Handle inputs
    num_edges  = 257;

    % Compute threshold intrensities for different regions
    nx = int32(options.nx);
    ny = int32(options.ny);
    
    thresholds_reg = zeros(ny, nx);
    y_center = int32(zeros(1, ny));
    x_center = int32(zeros(1, nx));
    for k=1:numel(subsections)
        iy = 1 + floor(double(k-1)/double(nx));
        jx = k - nx*(iy -1);
        y_center(iy) = int32(round(mean(subsections(k).ypixels)));
        x_center(jx) = int32(round(mean(subsections(k).xpixels)));
        Igr_reg = Igr(subsections(k).ypixels,subsections(k).xpixels,1);
        thresholds_reg(iy, jx) = get_threshold_intensity_for_one_subsection(iy, jx, Igr_reg,...  
                                          fracBlack_id, num_edges, options);
    end

    % Compute a matrix of threshold intrensities by a smooth interpolation
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

function display_histogram(Igr,Hist, edges,A_B,mu_B,sig_B,A_G,mu_G,sig_G,...
                           A_W,mu_W,sig_W,my_title,options)
%   clear; clc;
    gauss_B = zeros(size(edges));
    gauss_W = zeros(size(edges));
    gauss_S = zeros(size(edges)); % sum of 2 gaussians
    for i=1:numel(edges)
        gauss_B(i) = A_B*exp(-((edges(i)-mu_B)/sig_B)^2);
        gauss_G(i) = A_G*exp(-((edges(i)-mu_G)/sig_G)^2);
        gauss_W(i) = A_W*exp(-((edges(i)-mu_W)/sig_W)^2);
        gauss_S(i) = gauss_B(i) + gauss_G(i) + gauss_W(i);
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
    plot(edges,gauss_G,'color','r', 'LineWidth',3);
    hold on;
    plot(edges,gauss_S,'color','green', 'LineWidth',3);
    title(my_title);
    waitforbuttonpress;
    if options.closeAll
        close all;
    end

% -- --------------------------------------------------------------------------

function [ sig ] = find_right_width(ind, Hist, edges, A, mu)
    sig = 0.001;
    i = ind;
    while i < numel(edges)
        i = i+1;
        if Hist(i) > 0 && Hist(i) <= A/exp(1)
            sig = abs(edges(i) - mu);  % sig_W
            break;
        end
    end

% -- --------------------------------------------------------------------------

function [ sig ] = find_left_width(ind, Hist, edges, A, mu)
    sig = 0.03;
    i = ind;
    while i  > 1          
        i = i-1;    
        if Hist(i) > 0 && Hist(i) <= A/exp(1)
            sig = abs(edges(i) - mu);  % sig_W
            break;
        end
    end

% -- --------------------------------------------------------------------------

function[A_B, A_G, A_W, mu_B, mu_G, mu_W, sig_B, sig_G, sig_W, ind_B, ...
         ind_G, ind_W] = get_initial_guess(Hist, edges, options)
    N     = numel(edges);
    N2    = round(double(N)*0.4);    
    A_B   = max(Hist(1:N2));
    A_W   = max(Hist(N2:(N-1))); % exclude the last bucket 
    ind_B = min(find(Hist == A_B));
    ind_W = max(find(Hist == A_W));
    mu_B  = edges(ind_B);
    mu_W  = edges(ind_W);
%   disp(['Initial A_B, A_W, mu_B, mu_W=' num2str(A_B) ' '  num2str(A_W) ' ' num2str(mu_B)  ' ' num2str(mu_W)]);
    if ind_B == 1
        sig_B = 0.03;
    else
        sig_B = find_left_width(ind_B, Hist, edges, A_B, mu_B);
    end
    if ind_W == numel(edges)
        sig_W = 0.03;
    else
        sig_W = find_right_width(ind_W, Hist, edges, A_W, mu_W);
    end
    gauss_B = zeros(size(edges));
    gauss_W = zeros(size(edges));
    Hist1   = zeros(size(edges)); % sum of 2 gaussians
    for i=1:numel(edges)
        gauss_B(i) = A_B*exp(-((edges(i)-mu_B)/sig_B)^2);
        gauss_W(i) = A_W*exp(-((edges(i)-mu_W)/sig_W)^2);
        Hist1(i)   = Hist(i) - gauss_B(i) - gauss_W(i);
    end
    Hist1(numel(edges)) = 0.;
    Hist1(Hist1 < 0) = 0.;

    if abs(ind_W - ind_B) < 2
        if ind_B > 1
            ind_B = ind_B-1;
        end
        if ind_W < numel(edges) 
            ind_W = ind_W + 1;
        end
    end
    f_G   = fit(edges(ind_B:ind_W)',Hist1(ind_B:ind_W)','gauss1');
    A_G   = f_G.a1;
    mu_G  = f_G.b1;
    sig_G = f_G.c1;
    ind_G = min(find(edges > mu_B));
     
% -- --------------------------------------------------------------------------

function [A_B, A_W, mu_B, mu_W, sig_B, sig_W, ind_B, ind_W, Hist1] = update_edges(Hist, edges, A_G, mu_G, sig_G, ind_G)
    gauss_G = zeros(size(edges));
    Hist1   = zeros(size(edges)); % sum of 2 gaussians
    for i=1:numel(edges)
        gauss_G(i) = A_G*exp(-((edges(i)-mu_G)/sig_G)^2);
        Hist1(i)   = Hist(i) - gauss_G(i);
        if Hist1(i) < 0
            Hist1(i) = 0.;
        end
    end 

    if ind_G - 1 > 10
        f_B   = fit(edges(1:ind_G)',Hist1(1:ind_G)','gauss1');
        A_B   = f_B.a1;
        mu_B  = f_B.b1;
        sig_B = f_B.c1;
        ind_B = min(find(edges > mu_B));
    else
        ind_B = round((ind_G+1)/2);
        A_B   = Hist1(ind_B);                     
        mu_B  = edges(ind_B);                
        sig_B = sig_G;                         
    end

    if numel(edges) - ind_G > 10
        f_W   = fit(edges(ind_G:numel(edges))',Hist1(ind_G:numel(edges))','gauss1');
        A_W   = f_W.a1;
        mu_W  = f_W.b1;
        sig_W = f_W.c1;
        edges1= edges(ind_G:numel(edges));
        ind_W = ind_G + max(find(edges1 < mu_W));
    else
        ind_W = round((ind_G+numel(edges))/2);
        A_W   = Hist1(ind_W);
        mu_W  = edges(ind_W);
        sig_W = sig_G/2;
    end

% -- --------------------------------------------------------------------------

function [A_G, mu_G, sig_G, ind_G, Hist1] = update_center(Hist, edges, A_B, A_W, mu_B, mu_W, sig_B, sig_W, ind_B, ind_W)
    gauss_B = zeros(size(edges));
    gauss_W = zeros(size(edges));
    Hist1   = zeros(size(edges)); % sum of 2 gaussians
    for i=1:numel(edges)
        gauss_B(i) = A_B*exp(-((edges(i)-mu_B)/sig_B)^2);
        gauss_W(i) = A_W*exp(-((edges(i)-mu_W)/sig_W)^2);
        Hist1(i)   = Hist(i) - gauss_B(i) - gauss_W(i);
        if Hist1(i) < 0
            Hist1(i) = 0.;
        end
    end

    if ind_W - ind_B > 10
        f_G   = fit(edges(ind_B:ind_W)',Hist1(ind_B:ind_W)','gauss1');
        A_G   = f_G.a1;
        mu_G  = f_G.b1;
        sig_G = f_G.c1;
        ind_G = min(find(edges > mu_G));
    else
        ind_G = round((ind_B+ind_W)/2);
        A_G   = Hist1(ind_G);                 
        mu_G  = edges(ind_G);
        sig_G = (sig_B+sig_W)/2;
    end

% -- --------------------------------------------------------------------------

function [threshold, threshold2, A_B, mu_B, sig_B, A_G, mu_G, sig_G, A_W, mu_W, sig_W] = ...
          compute_threshold_intensity(i, j, Igr, Hist, edges, options)
    display_iterations = 0;
    method = 2; % may be = 1 or 2

    [A_B, A_G, A_W, mu_B, mu_G, mu_W, sig_B, sig_G, sig_W, ind_B, ind_G, ind_W] = ...
        get_initial_guess(Hist, edges, options);
    if options.hist & (options.dispOn | (i == options.sx & j == options.sy))
        display_histogram(Igr, Hist, edges, A_B, mu_B, sig_B, ...
                          A_G, mu_G, sig_G, A_W, mu_W, sig_W, ...
                          'Initial guess', options);
    end
    if options.verbose
        disp(['Initial guess [A_B, A_G, A_W, mu_B, mu_G, mu_W, sig_B, sig_G, sig_W=' ...
              num2str([A_B, A_G, A_W, mu_B, mu_G, mu_W, sig_B, sig_G, sig_W])]);
    end
    NUM_ITER = 3;
    disp(['NUM_ITER=' num2str(NUM_ITER)]);
    for it=1:NUM_ITER
        [A_G, mu_G, sig_G, ind_G, Hist1] = update_center(Hist, edges, A_B, A_W, mu_B, mu_W, sig_B, sig_W, ind_B, ind_W);
        disp(['  A_G=' num2str(A_G) ' mu_G=' num2str(mu_G) ' sig_G=' num2str(sig_G)]);
        if display_iterations | (options.hist & (i == options.sx & j == options.sy))
            display_histogram(Igr, Hist1, edges, A_B, mu_B, sig_B, ...
                              A_G, mu_G, sig_G, A_W, mu_W, sig_W, ...
                              ['Iter ' num2str(it) ' center'], options);
        end
        [A_B, A_W, mu_B, mu_W, sig_B, sig_W, ind_B, ind_W, Hist1] = update_edges(Hist, edges, A_G, mu_G, sig_G, ind_G);
        disp(['  A_B=' num2str(A_B) ' mu_B=' num2str(mu_B) ' sig_B=' num2str(sig_B)]);
        disp(['  A_W=' num2str(A_W) ' mu_W=' num2str(mu_W) ' sig_W=' num2str(sig_W)]);
        if display_iterations | (options.hist & (i == options.sx & j == options.sy))
            display_histogram(Igr, Hist1, edges, A_B, mu_B, sig_B, ...
                              A_G, mu_G, sig_G, A_W, mu_W, sig_W, ...
                              ['Iter ' num2str(it) ' edges'], options);
        end
    end
    % Threshold is a position of minimum of the smoothed histogram between 2 maxima
    if method == 1 % use minimum of hist between 2 maxima
        smoothed_Hist = smooth(Hist);
        smoothed_Hist_short = smoothed_Hist(ind_B:ind_W);
        min_val = min(smoothed_Hist_short);
        ind_min = min(find(smoothed_Hist_short == min_val));
        threshold = edges(ind_B + ind_min);
        threshold2 = edges(round((ind_B + ind_min + ind_W)/2));
    else % use approximations of hist by gaussian distributions
        gauss_B = zeros(size(edges));
        gauss_G = zeros(size(edges));
        gauss_W = zeros(size(edges));
        abs_diff= 10000*ones(size(edges));
        % Determine threshold for neural cell boundaries detection
        for i=ind_B:ind_G
            gauss_B(i) = A_B*exp(-((edges(i)-mu_B)/sig_B)^2);
            gauss_G(i) = A_G*exp(-((edges(i)-mu_G)/sig_G)^2);
            abs_diff(i)= abs(gauss_B(i) - gauss_G(i));
        end
        min_val = min(abs_diff);
        ind_min = min(find(abs_diff == min_val));
        threshold = edges(ind_min);
        % Determine threshold for dark structures detection
        abs_diff= 10000*ones(size(edges));
        for i=ind_G:ind_W
            gauss_W(i) = A_W*exp(-((edges(i)-mu_W)/sig_W)^2);
            gauss_G(i) = A_G*exp(-((edges(i)-mu_G)/sig_G)^2);
            abs_diff(i)= abs(gauss_W(i) - gauss_G(i));
        end
        min_val = min(abs_diff);
        ind_min = max(find(abs_diff == min_val));
        threshold2 = edges(ind_min);
    end
    if options.verbose
        disp(['method=' num2str(method) ' threshold=' num2str(threshold) ...
              ' threshold2=' num2str(threshold2)]);
        disp(' ' );
    end

% -- --------------------------------------------------------------------------

function threshold = get_threshold_intensity_for_one_subsection(i, j, Igr, ...
    fracBlack_id, num_edges, options);
    global Hist;
    global edges;
    Hist  = imhist(Igr,num_edges)';      % counts
    edges = compute_edges(num_edges); % using a custom function for edges
    disp(['options.hist=' num2str(options.hist)]);
    if options.fracBlack == 0 | options.fracBlack2 == 0
        [auto_threshold, auto_threshold2, A_B, mu_B, sig_B, A_G, mu_G, sig_G, A_W, mu_W, sig_W] = ...
            compute_threshold_intensity(i, j, Igr, Hist, edges, options);
            disp(['  A_B=' num2str(A_B) ' mu_B=' num2str(mu_B) ' sig_B=' num2str(sig_B)]);
            disp(['  A_G=' num2str(A_G) ' mu_G=' num2str(mu_G) ' sig_G=' num2str(sig_G)]);
            disp(['  A_W=' num2str(A_W) ' mu_W=' num2str(mu_W) ' sig_W=' num2str(sig_W)]);
            disp([' (automatic) threshold='  num2str(auto_threshold)]);
            disp([' (automatic) threshold2=' num2str(auto_threshold2)]);
            disp([' mu_B + sig_B=' num2str(mu_B + sig_B)]);
            disp([' mu_B + 2*sig_B=' num2str(mu_B + 2*sig_B)]);
            disp([' mu_W - sig_W=' num2str(mu_W - sig_W)]);
            disp([' mu_W - 2*sig_W=' num2str(mu_W - 2*sig_W)]);

        if options.fracBlack == 0
            threshold = max(auto_threshold, mu_B + sig_B);
            if options.verbose
                disp(['i=' num2str(i) ' j=' num2str(j) ' (automatic) threshold=' num2str(threshold)]);
            end
        end
        if options.fracBlack2 == 0
            threshold2 = auto_threshold2;
            threshold2 = max(auto_threshold2, mu_W - sig_W);
            if options.verbose
                disp(['i=' num2str(i) ' j=' num2str(j) ' (automatic) threshold2=' num2str(threshold2)]);
            end
        end
    else
        if options.verbose
            disp(['i=' num2str(i) ' j=' num2str(j) ' threshold=' num2str(threshold) 'threshold2=' num2str(threshold2)]);
        end
    end
    if options.hist & ((options.nx == 1 & options.ny == 1) | (i == options.sx & j == options.sy))
        display_histogram(Igr, Hist, edges, A_B, mu_B, sig_B, ...
                          A_G, mu_G, sig_G, A_W, mu_W, sig_W, ...
                          'Histogram of intensity - final', options);
    end
    Hsum  = sum(Hist');
    cum_hist = get_cumulative_histogram(Hist', Hsum);

    % Black-white image
    if  fracBlack_id == 1 & options.fracBlack > 0
        threshold = get_threshold_intensity(cum_hist, options.fracBlack, edges, options.verbose);
    elseif fracBlack_id == 2 & options.fracBlack2 > 0
        threshold = get_threshold_intensity(cum_hist, options.fracBlack2, edges, options.verbose);
    end
    if options.verbose
        disp(' ');
        disp(['threshold=' num2str(threshold)]);
        disp(' ');
    end

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

