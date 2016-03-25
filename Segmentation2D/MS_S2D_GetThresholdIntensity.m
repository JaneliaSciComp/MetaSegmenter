%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: compute a threshold intensity for each region
% and then a matrix of threshold intensities       
% using a smooth interpolation between regions
% (so that every pixel will get its own threshold)
%        

function [M_thr, M] = MS_S2D_GetThresholdIntensity(Igr, options)
    max_Igr = max(max(Igr));
    disp(['max_Igr=' num2str(max(max(Igr)))]);
    Igr1 = mat2gray(double(Igr)/double(max_Igr));

    if options.debug
        if exist('Segmentation2D.log', 'file') == 2
            delete('Segmentation2D.log');
        end
        diary('Segmentation2D.log'); % store a session log
        diary on;
    end

    disp(' ');
    disp('Computing initial guess to the thresholds matrix...');
    im_size = size(Igr1);
    if round(options.ny) >= 1 && round(options.nx) >= 1
        ny = round(options.ny);
        nx = round(options.nx);
    else
        ny = round(size(Igr,1)/options.minWin);
        nx = round(size(Igr,2)/options.minWin);
    end
    [M_thr,M_thr_upp,M_thr_low, M] = get_initial_guess_for_thresholds_matrix(Igr1, ...
                                  ny, nx, options);
    disp('Thresholds matrix for initial guess:');
    M
    if (options.dispOn2 > 0 && options.ny > 0 && options.nx > 0) || ...
        options.showIter == 2
        visualize_initial_guess(Igr1, M_thr, M_thr_upp, M_thr_low, options);
    end

    if round(options.nx) == 0 || round(options.ny) == 0
        % Use recursive approach to adaptive thresholding
        disp('Start refinement of the thresholds matrix ...');
        [M_thr, M] = refine_thresholds_matrix(Igr1, ny, nx, M_thr, ...
                                         M_thr_upp, M_thr_low, options);
        diary off;
    end

% -- --------------------------------------------------------------------------

function visualize_initial_guess(Igr1, M_thr, M_thr1, M_thr2, options)
    if options.debug == 2
        % Visualize a binary map for upper threshold
        Ibw1  = get_BW(Igr1, M_thr1, 1, options);
        MS_S2D_ShowImage(Ibw1, 'Binary map at upper threshold ', options);
        Ibw1 = produce_watershed_binary_map(Ibw1, 0, options);
        MS_S2D_ShowImage(Ibw1, 'Binary map at upper threshold after watershed)', options);
    end

    if options.debug == 2
        % Visualize a binary map for lower threshold
        Ibw2 = get_BW(Igr1, M_thr2, 1, options);
        MS_S2D_ShowImage(Ibw2, 'Binary map at lower threshold ', options);
        Ibw2 = produce_watershed_binary_map(Ibw2, 0, options);
        MS_S2D_ShowImage(Ibw2, 'Binary map at lower threshold after watershed)', options);
    end

    % Visualize a binary map for initial guess
    Ibw  = get_BW(Igr1, M_thr, 1, options);
    MS_S2D_ShowImageAndBoxes(Ibw, ['Initial guess (after bisection)'], '', ...
            [], [], [], [], options);
    Ibw = produce_watershed_binary_map(Ibw, 0, options);
    MS_S2D_ShowImage(Ibw, 'Initial guess (after bisection and watershed)', options);
    
% -- --------------------------------------------------------------------------

function Iws  = produce_watershed_binary_map(Ibw, remove_small_fragments, options)
    D = bwdist(Ibw);
    L = watershed(D);
    Iws = ones(size(D));
    % BW image with thin boundaries:
    Iws(L == 0) = 0;
    Iws = imfill(Iws, 'holes');
    if remove_small_fragments
        Iws = bwareaopen(Iws, double(options.minSize));
    end

% -- --------------------------------------------------------------------------

function counts = get_counts(vector, max_intensity)
    counts = zeros(1, max_intensity + 1);
    for i=0:max_intensity
        counts(i+1) = sum(round(vector) == i);
    end

% -- --------------------------------------------------------------------------

function counts = get_counts2(Igr, max_intensity)
    counts = zeros(1, max_intensity + 1);
    for i=0:max_intensity
        counts(i+1) = sum(sum(round(Igr) == i));
    end

% -- --------------------------------------------------------------------------

function [M_thr, M_thr1, M_thr2, M] = ...
          get_initial_guess_for_thresholds_matrix(Igr1, ny, nx, options)
    % Compute threshold intrensities for different subsections
    im_size = size(Igr1);
    dx = double(im_size(2))/double(nx);
    dy = double(im_size(1))/double(ny);
    M_thr  = zeros(im_size);
    M_thr1 = zeros(im_size);
    M_thr2 = zeros(im_size);
    M = zeros(ny, nx);

    thresholds  = zeros(ny, nx);
    disp(' ');
    disp(['size(Igr1)=' num2str(size(Igr1)) ' nx,ny=' num2str([nx ny]) ' num_subsections=' num2str(nx * ny) ]);

    fid = '';
    if length(options.outThr) > 0
        fid = fopen(options.outThr, 'w');
    end

    num_sect = nx * ny;
    for s = 1:num_sect
        iy = int32(mod(s, ny));
        if iy == 0
            iy = ny;
        end
        if s > ny
            jx = int32((s - iy)/ny)+1;
        else
            jx = 1;
        end
        disp([' s= ' num2str(s) '/' num2str(num_sect) ' iy=', num2str(iy) ' jx=' num2str(jx)]);


        [ymin, ymax] = get_bounds(iy, ny, dy, im_size(1));
        [xmin, xmax] = get_bounds(jx, nx, dx, im_size(2));

        if 0 % compute thresholds based on the entire image
            counts    = get_counts2(round(Igr1*255), 255);
            if options.a == 0 && options.q == 0
                [thr1,mu] = MS_S2D_GetOtsuThresholds(counts, 1);
            else
                thr1      = MS_S2D_GetEntropyBasedThresholds(counts, 1);
            end
            counts2   = get_counts2(round(Igr1*255), thr1);
            if options.a == 0 && options.q == 0
                [thr2,mu] = MS_S2D_GetOtsuThresholds(counts2, 1);
            else
                thr2      = MS_S2D_GetEntropyBasedThresholds(counts2, 1);
            end
        else % compute thresholds based on each individual window (this works better)
            counts    = get_counts2(round(Igr1(ymin:ymax, xmin:xmax)*255), 255);
            if options.a == 0 && options.q == 0
                [thr1,mu] = MS_S2D_GetOtsuThresholds(counts, 1);
            else
                thr1      = MS_S2D_GetEntropyBasedThresholds(counts, 1);
            end
            counts2   = get_counts2(round(Igr1(ymin:ymax, xmin:xmax)*255), thr1);
            if options.a == 0 && options.q == 0
                [thr2,mu] = MS_S2D_GetOtsuThresholds(counts2, 1);
            else
                thr2      = MS_S2D_GetEntropyBasedThresholds(counts2, 1);
            end
        end
        thr       = multithresh(Igr1, 2)*255;
        thr_upp   = round(max(thr1,  thr(2)));
        thr_low   = round(min(thr2, thr(1)));
 
        if length(options.membPr) > 0 && length(options.mitoPr) > 0
            thr_upp   = round(thr1);
            thr_low   = round(thr2);
        end

        if options.debug > 0      
            disp(['thr1=' num2str(thr1) ' thr2=' num2str(thr2) ...
                  ' thr=' num2str(round(thr))]);
        end

        Igr1b = Igr1([ymin:ymax], [xmin:xmax]);
        disp(' ');
        adjust2 = 1;
        if options.nx > 0 && options.ny > 0
            adjust2 = 0;
        end
        thr_bs = get_threshold_using_bisection(Igr1b, thr_upp, thr_low, ...
                                               adjust2, options);
        disp([' Ended at thr=' num2str(thr_bs)]);
        disp(' ');

        x1b = [1 size(Igr1b, 2)];
        y1b = [1 size(Igr1b, 1)];
        
        M_thr( ymin:ymax, xmin:xmax)= thr_bs;
        M_thr1(ymin:ymax, xmin:xmax)= thr_upp;
        M_thr2(ymin:ymax, xmin:xmax)= thr_low;
        M(iy, jx) = thr_bs;
    end

    if fid ~= ''
        fclose(fid);
    end

% -- --------------------------------------------------------------------------

function [ind_min, ind_max] = get_bounds(i, n, d, im_size)
    ind_min = round(1 + (i-1)*d);
    if i < n
        ind_max = round(i*d);
    else
        ind_max = im_size;
    end
    
% -- --------------------------------------------------------------------------

function Ibw = get_BW(Igr, M_thr, refine_BW, options)
    Ibw  = im2bw(Igr, 1);
    Ibw(Igr > M_thr/255.) = logical(1);  
    Ibw  = bwareaopen(Ibw, double(options.minSize));
    if refine_BW
%       Ibw = MS_S2D_AddBoundaryPadding(Ibw, 0);
        Ibw = imfill(Ibw, 'holes');
    end

% -- --------------------------------------------------------------------------

function thr_bs = get_threshold_using_bisection(Igr1, thr1, thr2, adjust2, options)
    disp('-----------------------------------------------------------------');
    disp(['Start bisection from thr=[' num2str([thr1 thr2]) ']']);

    Ibw1  = im2bw(Igr1, thr1/255.);
    Ibw2  = im2bw(Igr1, thr2/255.);

    x1b = [1 size(Igr1, 2)];
    y1b = [1 size(Igr1, 1)];

    dist = BWDistance(Ibw1, Ibw2, [] , [] , options);
    
    thr1_best = thr1;
    thr2_best = thr2;
    while (thr1_best - thr2_best) > 1
        thr = round((thr1_best + thr2_best)/2);
        Ibw = im2bw(Igr1, thr/255.);

        dist1 = BWDistance(Ibw1, Ibw, [] , [] , options);
        dist2 = BWDistance(Ibw2, Ibw, [] , [] , options);
        % NOTE: don't update Ibw1 and Ibw2!!! This will lead to wrong results
        if dist1 < dist2
%           Ibw1 = Ibw;       % updating Ibw1
            thr1_best = thr;
        else
            if adjust2
                Ibw2 = Ibw;       % updating Ibw1
            end
            thr2_best = thr;
        end
    end
    thr_bs = thr1_best;

% -- --------------------------------------------------------------------------

% Compute the background mito probability using bisection method
function mitoPrBg = compute_background_mitoPr(mitoPr, options)
    mitoPrBg = 0.5;
    my_strel = strel('disk', 1);
    max_mitoPr = max(max(mitoPr));
    if max_mitoPr > mitoPrBg
        Ibw = im2bw(mat2gray(mitoPr), mitoPrBg);
        CC = bwconncomp(Ibw);
        disp(['num regions before dilation/erosion=' num2str(CC.NumObjects)]);
%       MS_S2D_ShowImage(Ibw, 'Mito probabilities above background ', options);
        Ibw = imdilate(Ibw, my_strel);              
        Ibw = imerode( Ibw, my_strel);              
        Ibw1 = bwareaopen(Ibw, options.minSize);
        CC = bwconncomp(Ibw1);
        num_mito = CC.NumObjects;
        disp(['num mito regions=' num2str(num_mito)]);
        thr_max = mitoPrBg;
        thr_min = 0.1;
        while (thr_max-thr_min)*255 > 1
            thr = (thr_max + thr_min)/2;
            Ibw = im2bw(mat2gray(mitoPr), thr);
            Ibw = imdilate(Ibw, my_strel);           
            Ibw = imerode( Ibw, my_strel);           
            Ibw = bwareaopen(Ibw, options.minSize);
            CC  = bwconncomp(Ibw);
            nm  = CC.NumObjects;
            if nm > num_mito
                thr_min = thr;
            else
                thr_max = thr;
            end
        end
        mitoPrBg = thr_max;
        disp(['mitoPrBg=' num2str(mitoPrBg)]);
%       MS_S2D_ShowImage(Ibw1, 'Mito probabilities above background after dilation/erosion', options);
    end

% -- --------------------------------------------------------------------------

% Initial refinement - for ny x nx windows 
% Subsequent refinements - dividing each window into 4, 16, 64 etc.
%
function [M_thr, M] = refine_thresholds_matrix(Igr1, ny, nx, ...
                 M_thr, M_thr_upp, M_thr_low, options)
    im_size = size(Igr1);
   
    minGBF = 0.7;
    maxBIR = 0.3;

    subWin = [1:size(Igr1,1); 1:size(Igr1,2)];
    if numel(options.subWin) > 0
        limits = [ cellfun(@str2num,strsplit(options.subWin,','))];
        subWin = [limits(2):limits(4); limits(1):limits(3)];
    end
 
    membPr = [];
    if length(options.membPr) > 0
        membPr = MS_S2D_ReadProbabilities(options.membPr, 2);
        membPr = membPr(subWin(1,:),subWin(2,:));
    end

    mitoPr = [];
    if length(options.mitoPr) > 0
        mitoPr = MS_S2D_ReadProbabilities(options.mitoPr, 3);
        mitoPr = mitoPr(subWin(1,:),subWin(2,:));
    end

    % Estimate the threshold GBF and BIR (at thr_upp)
    Ibw = get_BW(Igr1, M_thr_upp, 0, options);
    [NT, NL, nl, NM, nm, TB, CSV, X,Y, IMF, TBV, GBF, BIR, oBox] = ...
          get_white_components(Igr1, Ibw, mitoPr, 0.1, minGBF, maxBIR, ...
          [1 im_size(2)], [1 im_size(1)], options);

    % Estimate the bacjkground mito probability
    mitoPrBg = 1;
    disp(' ');
    disp('Computing the background mito probability ...');

    if length(options.mitoPr) > 0
        mitoPrBg = compute_background_mitoPr(mitoPr, options);
    end

    % Use recursive approach to adaptive thresholding
    disp(' ');
    disp('=========================================================')
    disp(' ');
    disp(['Refinement of intensity thresholds (in ' num2str(options.numRef) ...
          ' cycles) ...']);

    dx = double(im_size(2))/nx;
    dy = double(im_size(1))/ny;
    
    for n = 1:uint8(options.numRef)
        [M_thr, M] = perform_refinement_cycle(n, Igr1, membPr, mitoPr, mitoPrBg, ...
                        minGBF, maxBIR, M_thr, M_thr_upp, M_thr_low, ...
                        nx, ny, dx, dy, options);
    end

    if options.showIter == 1 || options.showIter == 2
        Ibw  = get_BW(Igr1, M_thr, 1, options);
        MS_S2D_ShowImageAndBoxes(Ibw, ['BW after refinement cycle ' num2str(n) ], '', ...
            [], [], [], [], options);
        Ibw = produce_watershed_binary_map(Ibw, 0, options);
        MS_S2D_ShowImageAndBoxes(Ibw, ['BW after refinement cycle ' num2str(n) ], '', ...
            [], [], [], [], options);
    end

% -----------------------------------------------------------------------------

function M_thr_ref = eliminate_mito_boundaries(Igr1, M_thr_ref, mitoPr, mitoPrBg, ...
                     xb, yb, options)
    Igr1b = Igr1(yb(1):yb(2),xb(1):xb(2));
    M_thr_ref1b = M_thr_ref(yb(1):yb(2),xb(1):xb(2));
    mitoPr1b = mitoPr(yb(1):yb(2),xb(1):xb(2));
    Ibw1b = get_BW(Igr1b, M_thr_ref1b, 0, options);
    options.maxSize = options.minSize;
    [NT, NL, nl, NM, nm, TB, CSV, X,Y, IMF, TBV, GBF, BIR, oBox] = ...
        get_white_components(Igr1b, Ibw1b, mitoPr1b, mitoPrBg, 0, 0, ...
            [1 size(Igr1b,2)],[1 size(Igr1b,1)], options);

    Ibw2b = produce_watershed_binary_map(Ibw1b, 0, options);
    [NT, NL, nl, NM, nm, TB, CSV, X,Y, IMF, TBV, GBF, BIR, oBox] = ...
        get_white_components(Igr1b, Ibw2b, mitoPr1b, mitoPrBg, 0, 0, ...
            [1 size(Igr1b,2)],[1 size(Igr1b,1)], options);

    CC  = bwconncomp(Ibw2b);
    master_BW_map = uint8(zeros(size(Igr1b)));       
    num_mito_regions = 0;    
    for k = 1:CC.NumObjects
        if is_mito(Ibw1b, CC.PixelIdxList{k}, mitoPr1b, mitoPrBg) > options.minMitoFr
            num_mito_regions = num_mito_regions + 1;
            current_BW_map = uint8(zeros(size(Igr1b)));
            current_BW_map(CC.PixelIdxList{k}) = 255;
            current_BW_map = imdilate(current_BW_map,strel('disk', 1));                   
            M_thr_ref1b(master_BW_map & current_BW_map) = 0;
            master_BW_map = (master_BW_map | current_BW_map);
        end
    end
    M_thr_ref(yb(1):yb(2),xb(1):xb(2)) = M_thr_ref1b;

% -----------------------------------------------------------------------------

function visualize_refinement_results(ref_type, thr_orig, thr_ref, Igr1, ...
             M_thr1, M_thr, M_thr_ref, ...
             mitoPr, mitoPrBg, minGBF, maxBIR, xb, yb, options)

    % At upper threshold
    Ibw  = get_BW(Igr1, M_thr1, 0, options);
    thr_upp = mean(mean(M_thr1(yb(1):yb(2),xb(1):xb(2))));
    [NT, NL, nl, NM, nm, TB, CSV, X,Y, IMF, TBV, GBF, BIR, oBox] = ...
        get_white_components(Igr1, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);
    MS_S2D_ShowImageAndBoxes(Ibw, ['At upper thr =' num2str(thr_upp)], ...
      ref_type, [xb(1) yb(1) xb(2)-xb(1) yb(2)-yb(1)], ...
      [oBox(1) oBox(2) oBox(3)-oBox(1) oBox(4)-oBox(2)], X(1:nl), Y(1:nl), options);

    % Before refinement
    M_thr(yb(1):yb(2),xb(1):xb(2))= thr_orig;
    Ibw   = get_BW(Igr1, M_thr, 0, options);
    [NT, NL, nl, NM, nm, TB, CSV, X,Y, IMF, TBV, GBF, BIR, oBox] = ...
        get_white_components(Igr1, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);
    MS_S2D_ShowImageAndBoxes(Ibw, ['Before ' ref_type ' refinement ' ...
      ' thr=' num2str(thr_orig) ], ref_type, ...
      [xb(1)     yb(1)     xb(2)-xb(1)     yb(2)-yb(1)], ...
      [oBox(1),oBox(2),oBox(3)-oBox(1),oBox(4)-oBox(2)], X(1:nl), Y(1:nl), options);

    % After refinement
    M_thr(yb(1):yb(2),xb(1):xb(2))= thr_ref;
    disp(['In visualize_refinement_results thr_ref=' num2str(thr_ref)]);
    Ibw   = get_BW(Igr1, M_thr, 0, options);
    [NT, NL, nl, NM, nm, TB, CSV, X,Y, IMF, TBV, GBF, BIR, oBox] = ...
        get_white_components(Igr1, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);
    MS_S2D_ShowImageAndBoxes(Ibw, ['After ' ref_type ' refinement ' ...
      ' thr=' num2str(thr_ref) ], ref_type, ...
      [xb(1) yb(1) xb(2)-xb(1) yb(2)-yb(1)], ...
      [oBox(1),oBox(2),oBox(3)-oBox(1),oBox(4)-oBox(2)], X(1:nl), Y(1:nl), options);

    % After refinement, with watershed in the window
    Ibwb = produce_watershed_binary_map(Ibw(yb(1):yb(2),xb(1):xb(2)), 0, options);
    M_thr_ref1b =  M_thr_ref(yb(1):yb(2),xb(1):xb(2));
    Ibwb(M_thr_ref1b == 0) = 0;
    Ibw(yb(1):yb(2),xb(1):xb(2)) = Ibwb(:,:);
    [NT, NL, nl, NM, nm, TB, CSV, X,Y, IMF, TBV, GBF, BIR, oBox] = ...
        get_white_components(Igr1, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);
    MS_S2D_ShowImageAndBoxes(Ibw, ['After ' ref_type ' refinement ' ...
      ' thr=' num2str(thr_ref) ], ref_type, ...
      [xb(1) yb(1) xb(2)-xb(1) yb(2)-yb(1)], ...
      [oBox(1),oBox(2),oBox(3)-oBox(1),oBox(4)-oBox(2)], X(1:nl), Y(1:nl), options);

    % The enrire watershed so far
    Ibw   = get_BW(Igr1, M_thr_ref, 0, options);
    Ibw = produce_watershed_binary_map(Ibw, 0, options);
    Ibw = bwareaopen(Ibw, double(options.minSize));
    MS_S2D_ShowImageAndBoxes(Ibw, ['The entire watershed so far'...
      ' (after ' ref_type ')'], ref_type, [], [], [], [], options);

% -----------------------------------------------------------------------------

function [thr_ref, NL, NL3, NM, nm]  = ...
         refine_one_threshold(thr_max, thr_min, gal, ...
                              NT1, NL1, nl1, NM1, nm1, CSV1, X1, Y1, IMF1, M_thr1, ...
                              NL, NL3, NM, nm, M_thr, Igr1, xb, yb, ref_type, ...
                              mitoPr, mitoPrBg, minGBF, maxBIR, options)
    thr = thr_max;
    prev_NL = NL;
    prev_NL3 = NL3;
    prev_nm = nm;
    prev_NM1 = NM1;
    prev_thr_max = thr_max;
    prev_thr = thr_max;
    adjust2  = 0;
    if strcmp(ref_type, 'incrementation')
        prev_thr = thr_min;
        adjust2 = 1;
    end
    while  (thr_max - thr_min > 1) && ...
          ((strcmp(ref_type, 'incrementation') && ~gal || ...
           (strcmp(ref_type, 'decrementation') &&  gal && (NM1 < nm || NL <= NL3))))

        
        thr = get_threshold_using_bisection(Igr1(yb(1):yb(2),xb(1):xb(2)), thr_max, thr_min, adjust2, options);
%       if thr == thr_max
%           if strcmp(ref_type, 'decrementation')
%               thr = thr_min;
%           end
%           disp([' Breaking1 at thr_bs=' num2str(thr)]);
%           break;
%       end

        M_thr( yb(1):yb(2),xb(1):xb(2)) = thr;
        Ibw  = get_BW(Igr1, M_thr, 0, options);
        [NT, NL, nl, NM, nm, TB, CSV, X, Y, IMF, TBV, GBF, BIR, oBox] = ...
            get_white_components(Igr1, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);

        M_thr(yb(1):yb(2),xb(1):xb(2)) = thr - 1;
        Ibw = get_BW(Igr1, M_thr, 0, options);
        [NT3, NL3, nl3, NM3, nm3, TB3, CSV3, X3, Y3, IMF3, TBV3, GBF3, BIR3, oBox3] = ...
            get_white_components(Igr1, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);
        M_thr(yb(1):yb(2),xb(1):xb(2)) = thr;

        disp([' Current ' ref_type ' thr=' num2str(thr) ' NL=' num2str(NL) ' gal=' num2str(gal) ...
              ' NM1=' num2str(NM1) ' nm=' num2str(nm) ' NL3=' num2str(NL3)]);

        [gal, alignment1] = is_good_alignment(NT1, NL1, nl1, NM1, nm1, CSV1, ...
                                X1, Y1, IMF1, NT , NL , nl , NM , nm , CSV , ...
                                X , Y , IMF , Igr1, options);


        if strcmp(ref_type, 'incrementation')
            if ~gal                     
                thr_min = thr;
            elseif gal && NM1 < nm  
                if thr_max == thr
                   prev_thr_max = thr_max;
                   prev_NL = NL;
                   thr_max = round((thr_max + thr_min)/2);
                   disp(['   Case1: setting thr_max = prev_thr_max = ' num2str(prev_thr_max) ', thr_min = thr']);
                   continue;
                else
                    thr_max = thr;
                    disp(['   Case2: setting thr_max = thr = ' num2str(thr)]);
                end
            else
                disp([' Breaking2 at thr_bs=' num2str(thr)]);
                break;
            end
        elseif strcmp(ref_type, 'decrementation')
            if ~gal                  
                thr_max = prev_thr_max;
                thr_min = thr;
                thr = thr_max;
                NL = prev_NL;
                disp(['   Case3: setting thr_max = prev_thr_max = ' num2str(prev_thr_max) ', thr_min = thr=' num2str(thr)]);
                NL3 = prev_NL3;
                nm = prev_nm;
                NM1 = prev_NM1;
                continue;
            elseif NL <= NL3 || NM1 < nm
                if thr_max == thr
                   prev_thr_max = thr_max;
                   prev_NL = NL;
                   thr_max = thr_max - 1;                  
                   disp(['   Case4: setting thr_max = thr - 1 = ' num2str(thr_max)]);
                   prev_NM1 = NM1;
                   prev_nm = nm;
                   prev_NL3 = NL3;
                   continue;
                else
                   thr_max = thr;                
                   disp(['   Case5: setting thr_max = thr = ' num2str(thr)]);
                end
            else
                disp([' Breaking3 at thr_bs=' num2str(thr) ' since the goal is reached']);
                break;
            end
        end
        prev_NL = NL;
        prev_thr_max = thr_max;
        prev_thr = thr;
        prev_NM1 = NM1;
        prev_nm = nm;
        prev_NL3 = NL3;
    end
    thr_ref = thr;

% -----------------------------------------------------------------------------
function output_refinement_start_info(thr, thr_max, thr_min, ...
                                      ref_type, jx, iy, gal, NL1, ...
                                      NL, NL3, nl1, nl, nl3, NM1, NM, NM3, nm, ...
                                      TB1, CSV1, CSV, CSV3, alignment1, ...
                                      X1, Y1, X, Y, X3, Y3, IMF1, IMF, TBV1, GBF1, BIR1, options)
    disp(' ');
    disp('---------------------------------------------------');
    disp([' Start ' ref_type ' refinement for jx=' num2str(jx) ' iy=' num2str(iy) ]);  
    disp([' from thr=' num2str(thr) ' thr_max=' num2str(thr_max) ' thr_min=' num2str(thr_min)]);  
    disp([' options.incRef=' num2str(uint8(options.incRef)) ' options.numRef=' num2str(uint8(options.numRef))]); 
    disp([' gal=' num2str(gal) ' NL1=' num2str(NL1) ' NL=' num2str(NL) ...
          ' NL3=' num2str(NL3) ' nl=' num2str(nl) ' nl1=' num2str(nl1) ...
          ' NM1=' num2str(NM1) ' nm=' num2str(nm) ' TB1=' num2str(TB1) ]);
    disp([' CSV1=' num2str(CSV1(1:nl1)) ]);
    disp([' CSV =' num2str(CSV(1:nl)) ]);
    disp([' CSV3=' num2str(CSV3(1:nl3)) ]);
    disp([' alignment1=' num2str(alignment1)]);
    disp(['  X1 =' num2str(X1(1:nl1)) ]);
    disp(['  Y1 =' num2str(Y1(1:nl1)) ]);
    disp(['   X =' num2str(X(1:nl)) ]);
    disp(['   Y =' num2str(Y(1:nl)) ]);
    disp(['  X3 =' num2str(X3(1:nl3)) ]);
    disp(['  Y3 =' num2str(Y3(1:nl3)) ]);
    disp([' IMF1      =' num2str(IMF1) ]);
    disp([' IMF       =' num2str(IMF ) ]);
    disp([' TBV1      =' num2str(TBV1) ]);
    disp([' GBF1      =' num2str(GBF1) ]);
    disp([' BIR1      =' num2str(BIR1) ]);


% -----------------------------------------------------------------------------

function output_refinement_finish_info(thr, thr_ref, thr_orig, ...
                                       gal, NL1, nl1, NM1, ref_type, ...
                                       Igr1, M_thr, mitoPr, mitoPrBg, ...
                                       minGBF, maxBIR, xb,yb, options)
    M_thr(yb(1):yb(2),xb(1):xb(2)) = thr_ref;
    Ibw  = get_BW(Igr1, M_thr, 0, options);
    [NT, NL, nl, NM, nm, TB, CSV, X, Y, IMF, TBV, GBF, BIR, oBox] = ...
        get_white_components(Igr1, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);

    M_thr(yb(1):yb(2),xb(1):xb(2))= thr_orig;

    disp(' ');
    disp([' After ' ref_type ' refinement thr=' num2str(thr_ref) ...
          ' gal='  num2str(gal) ' NL =' num2str(NL) ...
          ' NL1 =' num2str(NL1) ' nl=' num2str(nl) ' nl1=' num2str(nl1) ...
          ' NM1 =' num2str(NM1) ' nm=' num2str(nm) '  TB=' num2str(TB)]);
    disp([' CSV =' num2str(CSV(1:nl))]);
    disp([' IMF =' num2str(IMF(1:nl))]);
    disp([' GBF =' num2str(GBF)]);
    disp([' BIR =' num2str(BIR)]);
    disp(' ');

% -----------------------------------------------------------------------------

function M = get_thresholds_matrix(nx, ny, dx, dy, M_thr)
    M = zeros([ny nx]);
    im_size = size(M_thr);
    for s = 1:(nx*ny)
        iy = int32(mod(s, ny));
        if iy == 0
            iy = ny;
        end
        if s > ny
            jx = int32((s - iy)/ny)+1;
        else
            jx = 1;
        end

        [yb(1) yb(2)] = get_bounds(iy, ny, dy, im_size(1));
        [xb(1) xb(2)] = get_bounds(jx, nx, dx, im_size(2));
        M(iy, jx) = mean(mean(M_thr(yb(1):yb(2), xb(1):xb(2))));
    end

% -----------------------------------------------------------------------------

function [M_thr_ref, M] = perform_refinement_cycle(cycle_id, Igr1, membPr, mitoPr,...
                          mitoPrBg, minGBF, maxBIR, ...
                          M_thr0, M_thr1, M_thr2, nx, ny, dx, dy, options)
    disp(' ');
    disp('=================================================');
    disp(' ');
    disp(['Start refinement cycle #' num2str(cycle_id)]);
    disp(' ');
    disp('=================================================');
    
    M = zeros([ny nx]);    % nx x ny matrix of thresholds 
    M_thr_ref = zeros(size(M_thr0));
    M_thr_ref(:,:) = M_thr0(:,:);
    M_thr     = zeros(size(M_thr0));
    im_size   = size(Igr1);
    
    for s = 1:(nx*ny) 
        iy = int32(mod(s, ny));
        if iy == 0
            iy = ny;
        end
        if s > ny
            jx = int32((s - iy)/ny)+1;
        else
            jx = 1;
        end

        [yb(1) yb(2)] = get_bounds(iy, ny, dy, im_size(1));
        [xb(1) xb(2)] = get_bounds(jx, nx, dx, im_size(2));

        % Before refinement
        Igr1b  = Igr1( yb(1):yb(2),xb(1):xb(2));

        thr     = round(mean(mean(M_thr0( yb(1):yb(2),xb(1):xb(2)))));
        thr_upp = round(mean(mean(M_thr1(yb(1):yb(2),xb(1):xb(2)))));
        thr_low = round(mean(mean(M_thr2(yb(1):yb(2),xb(1):xb(2)))));
        thr3    = thr-1;
        thr_ref = thr;
        thr_orig= thr;

        % Trial M_thr is M_thr1 everywhere, except the current inner box
        M_thr(:,:) = M_thr1(:,:);
        M_thr(yb(1):yb(2),xb(1):xb(2)) = thr;

        disp(' ');
        disp('=================================================');
        disp(' ');
        disp([' Refinement cycle= ' num2str(cycle_id) '/' num2str(options.numRef) ...
              ', section= ' num2str(s) '/' num2str(nx*ny) ...
              ', jx=' num2str(jx) ' iy=' num2str(iy) ]);
        disp(' ');
        disp([' thr_upp=' num2str(thr_upp) ' thr =' num2str(thr    ) ' thr_low=' num2str(thr_low)]);
        disp(' ');

        M_thr(yb(1):yb(2),xb(1):xb(2)) = thr_upp;
        Ibw1 = get_BW(Igr1, M_thr1, 0, options);
        [NT1, NL1, nl1, NM1, nm1, TB1, CSV1, X1, Y1, IMF1, TBV1, GBF1, BIR1, oBox1] = ...
            get_white_components(Igr1, Ibw1, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);

        disp([' NT1=' num2str(NT1) ' NL1=' num2str(NL1) ' nl1=' num2str(nl1) ...
              ' NM1=' num2str(NM1) ' nm1=' num2str(nm1) ' TB1=' num2str(TB1)]);
        disp([' CSV1=' num2str(CSV1) ]);
        disp([' IMF1=' num2str(IMF1) ]);
        disp(['   X1=' num2str(X1(1:nl1)) ]);
        disp(['   Y1=' num2str(Y1(1:nl1)) ]);
        disp([' GBF1=' num2str(GBF1) ]);
        disp([' BIR1=' num2str(BIR1) ]);

        disp(' ');
        disp(' ');

        M_thr(yb(1):yb(2),xb(1):xb(2)) = thr_orig;
        Ibw   = get_BW(Igr1, M_thr, 0, options);
        [NT, NL, nl, NM, nm, TB, CSV, X,Y, IMF, TBV, GBF, BIR, oBox] = ...
            get_white_components(Igr1, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);

        disp([' NT =' num2str(NT) ' NL =' num2str(NL)  ' nl=' num2str(nl) ...
              ' NM =' num2str(NM)  ' nm=' num2str(nm) ' TB=' num2str(TB)]);  
        disp([' CSV =' num2str(CSV) ]);
        disp([' IMF =' num2str(IMF) ]);
        disp(['   X =' num2str(X(1:nl)) ]);
        disp(['   Y =' num2str(Y(1:nl)) ]);
        disp([' GBF =' num2str(GBF) ]);
        disp([' BIR =' num2str(BIR)]);

        disp(' ');
        disp(' ');

        M_thr(yb(1):yb(2),xb(1):xb(2)) = thr3;
        Ibw3 = get_BW(Igr1, M_thr, 0, options);
        [NT3, NL3, nl3, NM3, nm3, TB3, CSV3, X3, Y3, IMF3, TBV3, GBF3, BIR3, oBox3] = ...
            get_white_components(Igr1, Ibw3, mitoPr, mitoPrBg, minGBF, maxBIR, xb,yb, options);

        disp([' NT3 =' num2str(NT3) ' NL3 =' num2str(NL3)  ' nl3=' num2str(nl3) ...
              ' NM3 =' num2str(NM3)  ' nm=' num2str(nm3) ' TB=' num2str(TB3)]);
        disp([' CSV3 =' num2str(CSV3) ]);
        disp([' IMF3 =' num2str(IMF3) ]);
        disp(['   X3 =' num2str(X3(1:nl3)) ]);
        disp(['   Y3 =' num2str(Y3(1:nl3)) ]);
        disp([' GBF3 =' num2str(GBF3) ]);
        disp([' BIR3 =' num2str(BIR3)]);

        M_thr(yb(1):yb(2),xb(1):xb(2)) = thr;

        disp(' ');
        disp(' ');

        [gal, alignment1] = is_good_alignment(NT1, NL1, nl1, NM1, nm1, CSV1, ...
                                X1, Y1, IMF1, NT , NL , nl , NM , nm , CSV , ...
                                X , Y , IMF , Igr1, options);

        disp(['gal=' num2str(gal) ' alignment1=' num2str(alignment1)]);
        disp(' ');

        if gal && (nm <= NM1 || NM1 == 0)  
            if options.verbose > 0
                disp('---------------------------------------------------');
                disp(['    Nothing to refine ' ...
                      ' for jx=' num2str(jx) ' iy=' num2str(iy) ]);   
            end
            M(iy, jx)  = thr_ref;
            M
        elseif ~gal               
            % Must perform incrementation refinement!
            num_iter = 0;
            while ~gal && num_iter < 1 
                thr_max = thr_upp;
%               if TB1 > 0
%                   thr_max = min(thr_upp, TB1);
%               end
                thr_min = thr;

                if options.verbose
                    output_refinement_start_info(thr, thr_max, thr_min, ...
                        'incrementation', jx, iy, ...
                        gal, NL1, NL, NL3, nl1, nl, nl3, NM1, NM, NM3, nm, TB1, ...
                        CSV1, CSV, CSV3, alignment1, ...
                        X1, Y1, X, Y, X3, Y3, IMF1, IMF, TBV1, GBF1, BIR1, options);

                end
                [thr_ref, NL, NL3, NM, nm] = ...
                    refine_one_threshold(thr_max, thr_min, gal, ...
                        NT1, NL1, nl1, NM1, nm1, CSV1, X1, Y1, IMF1, M_thr1, ...
                        NL, NL3, NM, nm, M_thr, Igr1, xb, yb, 'incrementation', ...
                        mitoPr, mitoPrBg, minGBF, maxBIR, options);

                M_thr_ref(yb(1):yb(2),xb(1):xb(2))= thr_ref;     
                    
%               M_thr_ref = eliminate_mito_boundaries(Igr1, M_thr_ref, mitoPr, mitoPrBg, xb, yb, options);

                disp([' Final thr_ref=' num2str(mean(mean(M_thr_ref(yb(1):yb(2),xb(1):xb(2)))))]);
                M(iy, jx)  = thr_ref;
                M

                if options.verbose
                    output_refinement_finish_info(thr, thr_ref, thr_orig, ...
                           gal, NL1, nl1, NM1, ...
                           'incrementation', Igr1, M_thr, mitoPr, mitoPrBg, ...
                           minGBF, maxBIR, xb,yb, options);
                end

                if options.showIter == 2 
                    visualize_refinement_results('incrementation', thr_orig, ...
                        thr_ref, Igr1, M_thr1, M_thr, M_thr_ref, ...
                        mitoPr, mitoPrBg, minGBF, maxBIR, xb, yb, options);
                end

                [gal, alignment1] = is_good_alignment(NT1, NL1, nl1, NM1, nm1, ...
                                        CSV1, X1, Y1, IMF1, NT , NL , nl , ...
                                        NM , nm , CSV , X , Y , IMF , Igr1, options);
                num_iter = num_iter + 1;
                thr = thr_ref;
            end
            M_thr(yb(1):yb(2),xb(1):xb(2))= thr_orig;

        elseif uint8(options.incRef) == 0 && gal && (nm > NM1 || NL == NL3)
            % Can perform decrementation refinement
            thr_max = thr;
%           thr_min = thr_low;
            thr_min = 0;

            if options.verbose
                output_refinement_start_info(thr, thr_max, thr_min, ...
                         'decrementation', jx, iy, ...
                         gal, NL1, NL, NL3, nl1, nl, nl3, NM1, NM, NM3, nm, TB1, ...
                         CSV1, CSV, CSV3, alignment1, ...
                         X1, Y1, X, Y, X3, Y3, IMF1, IMF, TBV1, GBF1, BIR1, options);
            end

            [thr_ref, NL, NL3, NM, nm] = ...
                refine_one_threshold(thr_max, thr_min, gal, ...
                    NT1, NL1, nl1, NM1, nm1, CSV1, X1, Y1, IMF1, M_thr1, ...
                    NL, NL3, NM, nm, M_thr, Igr1, xb, yb, 'decrementation', ...
                    mitoPr, mitoPrBg, minGBF, maxBIR, options);
        

            M_thr_ref(yb(1):yb(2),xb(1):xb(2))= thr_ref;               

%           M_thr_ref = eliminate_mito_boundaries(Igr1, M_thr_ref, mitoPr, mitoPrBg, xb, yb, options);

            disp([' Final thr_ref=' num2str(mean(mean(M_thr_ref(yb(1):yb(2),xb(1):xb(2)))))]);
            M(iy, jx)  = thr_ref;
            M

            if options.verbose
                output_refinement_finish_info(thr, thr_ref, thr_orig, ...
                       gal, NL1, nl1, NM1, ...
                       'decrementation', Igr1, M_thr, mitoPr, mitoPrBg, ...
                       minGBF, maxBIR, xb,yb, options);
            end

            if options.showIter == 2
                visualize_refinement_results('decrementation', thr_orig, ...
                    thr_ref, Igr1, M_thr1, M_thr, M_thr_ref, ...
                    mitoPr, mitoPrBg, minGBF, maxBIR, xb, yb, options);
            end

            M_thr(yb(1):yb(2),xb(1):xb(2))= thr_orig;
        end
        if jx == options.sx && iy == options.sy
            visualize_refinement_results('uncknown', thr_orig, thr_ref, ...
                Igr1, M_thr1, M_thr, M_thr_ref, ...
                mitoPr, mitoPrBg, minGBF, maxBIR, xb, yb, options);
        end
        % This updating makes the result more accurate, although 
        % it makes it dependent on the order in which windows were processed              
        M_thr(:,:) = M_thr_ref(:,:);
    end
    disp(['Final thresholds matrix=']);
    M_final = get_thresholds_matrix(nx, ny, dx, dy, M_thr_ref)

% -----------------------------------------------------------------------------

function oBox = update_outer_box(thisBB, oBox)
    if  oBox(1) > thisBB(1)
        oBox(1) = thisBB(1);
    end
    if  oBox(2) > thisBB(2)
        oBox(2) = thisBB(2);
    end
    if  oBox(3) < thisBB(1) + thisBB(3)
        oBox(3) = thisBB(1) + thisBB(3);
    end
    if  oBox(4) < thisBB(2) + thisBB(4)
        oBox(4) = thisBB(2) + thisBB(4);
    end

% -----------------------------------------------------------------------------

function bc = is_boundary_component(thisBB, sz)
    bc = 0;
    if thisBB(1) == 1 || ...
       thisBB(2) == 1 || ...
       thisBB(1) + thisBB(3) >= sz(2) || ... 
       thisBB(2) + thisBB(4) >= sz(1)
       bc = 1;
    end

% -----------------------------------------------------------------------------
% Returns the fraction of pixels where mitoPr > mitoPrBg 
function imp = is_mito(Ibw, pixelsList, mitoPr, mitoPrBg)
    sz = size(Ibw);
    imp = 0;
    if numel(mitoPr) > 0
        npix_mito = sum(mitoPr(ind2sub(sz, pixelsList)) > mitoPrBg);
        if npix_mito > 0
            imp = npix_mito/numel(pixelsList);
        end
    end

% -----------------------------------------------------------------------------

function [NT, CSV] = get_num_components(Ibw, xb, yb, options)
    sz = size(Ibw);
    Ibw = bwareaopen(Ibw, double(options.minSize-1));
    CC  = bwconncomp(Ibw, 8);
    BB = regionprops(CC,'Boundingbox');
    CSV_unsort = [];
    if numel(xb) > 0 && numel(yb) > 0
        for k = 1:CC.NumObjects
            thisBB = BB(k).BoundingBox;  % 4 elements: x, y, width, height
            if xb(1) <= thisBB(1) + thisBB(3) && xb(2) >= thisBB(1) && ...
               yb(1) <= thisBB(2) + thisBB(4) && yb(2) >= thisBB(2)

                CSV_unsort = [CSV_unsort numel(CC.PixelIdxList{k})];
            end
        end
    else
        CSV_unsort = cellfun(@numel, CC.PixelIdxList);            
    end
    NT  = numel(CSV_unsort);
    CSV = sort(CSV_unsort, 'descend');

% -----------------------------------------------------------------------------

function [TBV, GBF, BIR] =  good_boundary_fractions(Igr, Ibw, options)
    Igr    = MS_S2D_AddBoundaryPadding(Igr, 0);
    Ibw    = MS_S2D_AddBoundaryPadding(Ibw, 0);

    if 0
        MS_S2D_ShowImageAndBoxes(Igr, 'Grayscale mage for computing boundary fraction', ...
            '', [], [], [], [], options);
        Irgb = cat(3, Igr, Igr, Igr);
        Irgb1 = Irgb(:,:,1);
        Irgb1(Ibw == 0) = 255;
        Irgb(:,:,1) = Irgb1;
        Irgb2 = Irgb(:,:,2);
        Irgb2(Ibw == 0) = 0;
        Irgb(:,:,2) = Irgb2;
        Irgb(:,:,3) = Irgb2;
%       Irgb(:,:,1) = Irgb1;
        MS_S2D_ShowImageAndBoxes(Irgb, 'Image for computing boundary fraction', ...
            '', [], [], [], [], options);
    end

    [B, L] = bwboundaries(Ibw,8,'noholes');
    TBV = [];
    GBF = []; % fraction of region boundary occupied by 'dark' pixels
    BIR = []; % ratio of the mean intensities of the dark and light boundary pixels
%   disp(['numel(B)=' num2str(numel(B))]);
    for i=1:numel(B)  % i = component id
        b = B{i};
        I1= Igr(b(:,1),b(:,2))*255;
        I = round(reshape(I1,1,numel(I1))); % make it a vector
        try
            counts = get_counts(I, max(I));
            if options.a == 0 && options.q == 0
                [thr, mu]  = MS_S2D_GetOtsuThresholds(counts, 1);
            else
                thr = MS_S2D_GetEntropyBasedThresholds(counts, options.a, options.q)
            end
        catch
            thr = max(I);
        end
        GBF = [GBF sum(I <= thr)/numel(I)];
        median1 = median(I(I <= thr));
        median2 = median(I(I >= thr));
        BIR   = [BIR median1/median2];
        TBV   = [TBV thr];
    end

% -----------------------------------------------------------------------------

function NM = get_num_mito(IMF_unsorted, CSV, options)
    IMS_unsorted = IMF_unsorted .* CSV; % vector of the #'s of mito pixels
    [IMS, idx] = sort(IMS_unsorted, 'descend');
    IMF = IMF_unsorted(idx);  % sor
    NM = 0;
    if IMF(1) > options.minMitoFr
        NM = 1;
        % Compute the index preceding the largest drop in the # of mito pixels
        if numel(IMS) == 2
            if IMF(2) > options.minMitoFr && IMF(2) > IMF(1)/2
                NM = 2;
            end
        else % numel(IMS) > 2
            max_drop_ind = 1;
            max_drop_ratio = 0;
            for j=2:(numel(IMS)-1)
                 if IMF(j) < options.minMitoFr
                     break;
                 end
                 drop_ratio = (IMS(j-1) - IMS(j))/IMS(j);
                 if max_drop_ratio < drop_ratio
                     max_drop_ratio = drop_ratio;
                     max_drop_ind = j-1;
                 end
            end
            NM = max_drop_ind;
        end
    end

% -----------------------------------------------------------------------------

% Compute the total # of components and the # of good ones
function [NT, NL, nl, NM, nm, TB, CSV, X, Y, IMF, TBV, GBF, BIR, oBox] = ...
          get_white_components(Igr, Ibw, mitoPr, mitoPrBg, minGBF, maxBIR, xb, yb, options)
    % NT = total number of white components overlapping with a givem window
    % nl = number of large (both mito and non-mito) components
    % NL = number of large, but only non-mito components
    % NM = number of the 'true' mito components
    % nm = number of all 'fragments' with mito fraction above a threshold
    % TB = threshold intensity just below the bifurcation of 'bad' boundary
    sz = size(Igr);

    % Process components included in the outer box; label them
    cc  = bwconncomp(Ibw, 8);
    nc  = cc.NumObjects;
    CSV_unsort = []; % component sizes vector
    IMF_unsort = []; % vector of the fraction of pixels with high mito probability/get_white_components
    X_unsort   = [];
    Y_unsort   = [];
    BB = regionprops(cc,'Boundingbox');
    oBox = [xb(1) yb(1) xb(2) yb(2)];
    L = zeros(sz); % labels indicating the size of component
    NM = 0;
    if numel(xb) == 2 && numel(yb) == 2
        for k = 1:nc
            thisBB = BB(k).BoundingBox;  % 4 elements: x, y, width, height
            
            if xb(1) <= thisBB(1) + thisBB(3) && xb(2) >= thisBB(1) && ...
               yb(1) <= thisBB(2) + thisBB(4) && yb(2) >= thisBB(2)

                if numel(cc.PixelIdxList{k}) > options.minSize
                    % Label fragments of large components
                    % Each label == component size
                    L(cc.PixelIdxList{k}) = numel(cc.PixelIdxList{k});
                end

                % Update outerBox
                oBox = update_outer_box(thisBB, oBox);
            end
        end

        % Process only large components overlapping with inner box
        Ibw1b    =    Ibw(yb(1):yb(2),xb(1):xb(2));
        L1b      =      L(yb(1):yb(2),xb(1):xb(2));
        mitoPr1b = [];
        if numel(mitoPr) > 0
            mitoPr1b = mitoPr(yb(1):yb(2),xb(1):xb(2));
        end

%       Ibw1b  = MS_S2D_AddBoundaryPadding(Ibw1b, 0);
        CC     = bwconncomp(Ibw1b);
        NC     = CC.NumObjects;
        disp(['Num white components in the inner box=' num2str(NC)]);
        rpc    = regionprops(CC,'Centroid');
        coords = cat(1, rpc.Centroid);

        for K = 1:NC
            % Pick only the fragments of large components
            if mean(mean(L1b(CC.PixelIdxList{K}))) >= options.maxSize
                
                CSV_unsort = [CSV_unsort numel(CC.PixelIdxList{K})];
                IMF_unsort = [IMF_unsort is_mito(Ibw1b, CC.PixelIdxList{K}, mitoPr1b, mitoPrBg)];
                X_unsort   = [X_unsort (coords(K,1) + xb(1)-1)];
                Y_unsort   = [Y_unsort (coords(K,2) + yb(1)-1)];
            end
        end
    else
        CSV_unsort = cellfun(numel, CC.PixelIdxList);
%       IMF_unsort = [IMF_unsort is_mito(Ibw, CC.PixelIdxList{k}, mitoPr)];
    end
    clear L;

    NT = numel(CSV_unsort);
    [CSV, idx] = sort(CSV_unsort, 'descend');
    IMF = IMF_unsort(idx);  % sorted IMF
    X   =   X_unsort(idx);
    Y   =   Y_unsort(idx);

    nl = numel(CSV);                          % both mito and non-mito
    NL = sum(IMF <  options.minMitoFr); % # large non-mito regions
    nm = sum(IMF >= options.minMitoFr); % # of all mito 'fragments'
    NM = sum(IMF >= options.maxMitoFr); % number of the true mito regions

    [TBV, GBF, BIR] = good_boundary_fractions(Igr(yb(1):yb(2),xb(1):xb(2)), ...
                                              Ibw(yb(1):yb(2),xb(1):xb(2)), options);
    TB = 0;
    if numel(GBF) >= nl &&  numel(BIR) >= nl
        inds = find(GBF(1:nl) >= 0.35 & BIR(1:nl) <= 0.35);
        if numel(inds) >= 1
            TB = min(TBV(inds));
        end
    end

% -----------------------------------------------------------------------------

function [gal, alignment1] = is_good_alignment(NT1, NL1, nl1, NM1, nm1, CSV1, ...
                                 X1, Y1, IMF1, NT , NL , nl , NM , nm , CSV , ...
                                 X,  Y,  IMF , Igr1, options)
    alignment1 = AlignCSVs(CSV1(1:nl1), CSV(1:nl), X1, Y1, X, Y, nl1, nl, options);
    gal = 1;              
    for i = 1:nl
        % Determine if the alignment is good
        ids1 = find(alignment1 == i);
        if numel(ids1) > 1
            gal = 0;
            for j=1:numel(ids1)
                % Alignment is still good if 
                % - at least one of merging regions is mito
                % - at least one of merging regions is of small size
                if IMF1(ids1(j)) > options.minMitoFr || CSV1(ids1(j)) < options.maxSize
                    gal =1;
                    break
                end
            end
        end
    end

% -----------------------------------------------------------------------------

function MS_S2D_ShowImageAndBoxes(I, image_title, ref_type, iB, oB, X, Y, options)
    if strcmp(ref_type, 'decrementation')
        color = [0 0.5 0];
    else
        color = 'r';
    end

    figure
    imshow(I);
    hold on;
    % Inner box
    if numel(iB) > 0
        rectangle('Position',[iB(1),iB(2),iB(3),iB(4)],'EdgeColor',color,...
              'LineWidth',2, 'LineStyle', '-');
        hold on;
    end

    % Outer box
    if numel(oB) > 0
        rectangle('Position',[oB(1),oB(2),oB(3),oB(4)],'EdgeColor',color,...
                  'LineWidth',2, 'LineStyle', '--');
        hold on;
    end

    if numel(X) > 0 && numel(Y) == numel(X)
        if  strcmp(ref_type, 'decrementation')
            plot(X, Y, 'g+');
        else
            plot(X, Y, 'r+');
        end
        hold on;
    end

    title(image_title);
    impixelinfo;

    waitforbuttonpress;
    if options.closeAll
        close all;
    end

% -----------------------------------------------------------------------------

function dist = CSVDistance(CSV1, CSV2)
    NT = min(numel(CSV1), numel(CSV2));
    dist = 0.;
    for i=1:NT
        dist = dist + (CSV1(i) - CSV2(i))^2;
    end
    dist = sqrt(dist);

% -----------------------------------------------------------------------------

function [BWdist] = BWDistance(Ibw1, Ibw2, xb, yb, options)
    % Compute a distance between two binary images
%   disp(['size(Ibw1)=' num2str(size(Ibw1)) ' size(Ibw2)=' num2str(size(Ibw2))]);
    [NT1, CSV1] = get_num_components(Ibw1, xb, yb, options);
    [NT2, CSV2] = get_num_components(Ibw2, xb, yb, options);

    NT = min(NT1, NT2);                    
    BWdist = CSVDistance(CSV1, CSV2);

% -----------------------------------------------------------------------------

function alignment1 = AlignCSVs(CSV1, CSV, X1, Y1, X, Y, nl1, nl, options)
    % alignedCSV1 =      array of ids of aligned components for the upper threshold
    % alignedIDS1 = cell array of ids of aligned components for the curr  threshold
    disp(' ');
    NT1 = numel(CSV1);
    NT  = numel(CSV );
%   disp(['In AlignCSVs: nl1=' num2str(nl1) ' NT1=' num2str(NT1) ' CSV1=' num2str(CSV1)]);
%   disp(['In AlignCSVs: nl =' num2str(nl ) ' NT =' num2str(NT ) ' CSV =' num2str(CSV )]);
    N1 = nl1 + nl1*(2*NT1 - 1 - nl1)/2;
    alignedCSV1 = zeros(1, N1);
    alignedX1   = zeros(1, N1);
    alignedY1   = zeros(1, N1);
    alignedIDS1 =  cell(   N1, 1);
    alignment1  = zeros(1, nl1);
    alignment   = zeros(1, nl );
    if options.debug == 2
        for i=1:nl
            disp(['k=' num2str(i) ' X,Y=' num2str(round([X(i) Y(i)]))]);
        end
        disp(' ');
    end
    for i=1:nl1
        alignedCSV1(i) = CSV1(i);
        alignedX1(i)   = X1(i);
        alignedY1(i)   = Y1(i);
        alignedIDS1{i} = i;
        if options.debug == 2
            disp(['k=' num2str(i) ' X1,Y1=' num2str(round([X1(i) Y1(i)]))]);
        end
    end
%   disp(['   Initial alignedCSV1=' num2str(alignedCSV1)]);
    k = nl1;
    for i=1:nl1
        for j=(i+1):NT1
            k = k + 1;
            alignedCSV1(k) = CSV1(i) + CSV1(j);
            alignedX1(k)   = (X1(i)*CSV1(i) + X1(j)*CSV1(j))/alignedCSV1(k);
            alignedY1(k)   = (Y1(i)*CSV1(i) + Y1(j)*CSV1(j))/alignedCSV1(k);
            if options.debug == 2
                disp(['k=' num2str(k) ' X1,Y1=' num2str(round([alignedX1(k) alignedY1(k)])) '   i,j=' num2str([i j])]);
            end
            alignedIDS1{k} = [i j];
        end
    end

    % Populate the dsitance matrix
%   disp(['nl1=' num2str(nl1) ' alignedCSV1_orig=' num2str(alignedCSV1)]);
%   disp(['N1=' num2str(N1) ' alignedCSV1_sorted=' num2str(alignedCSV1)]);
%   disp(['nl =' num2str(nl ) '        CSV_sorted=' num2str(CSV)]);
%   disp('alignedIDS1=');
%   for i=1:N1
%       disp([' (i=' num2str(i) ') ' num2str(alignedIDS1{i})]);
%   end
    N1 = numel(alignedCSV1);
    N  = numel(CSV);
    distMatr = zeros(N, N1);
    for i=1:N
        for j=1:N1
            distMatr(i,j) = sqrt(double((X(i)-alignedX1(j)))^2 + double((Y(i)-alignedY1(j)))^2);
        end
    end
    if options.debug == 2
        distMatr
    end
    % Compute the alignment vectors
    min_ind = 0;
    for n=1:N
        min_value =  min(min(distMatr));    
        min_ind = find(distMatr == min_value);
        if numel(min_ind) > 1
            min_ind = min_ind(1);
        end
        row = mod(min_ind, N);
        if row == 0
            row = N;
        end
        col = (min_ind - row)/N + 1;
%       disp(['   i=' num2str(i) ' min_ind=' num2str(min_ind) ' row=' num2str(row) ...
%             ' col=' num2str(col) ' numel(alignedIDS1)=' num2str(numel(alignedIDS1))]);
        if ~isinf(min_value) && min_value <= CSV(row)
            alignment(row) = 1;
            ids = alignedIDS1{col};
%           disp(['   row=' num2str(row) ' col=' num2str(col) ' ids=' num2str(ids) ' component size=' num2str(CSV(row))]);
%           disp(['   i=' num2str(i) ' min_ind=' num2str(min_ind) ' aligned_ids=' num2str(ids)]);
            for k=1:numel(ids)
                if ids(k) <= nl1
                    alignment1(ids(k)) = row;
                end
            end
%           disp(['    Now alignment1=' num2str(alignment1)]);
            distMatr(row,:) = Inf;
            distMatr(:,col) = Inf;
            for j=1:N1
                ids1 = alignedIDS1{j};
                for k=1:numel(ids1)
                    if ids1(k) <= nl1
                        if ids1(k) == ids(1) || (numel(ids) == 2 && ids1(k) == ids(2))
                            distMatr(:,j) = Inf;
                        end
                    end
                end
            end
%           distMatr
        end
    end
%   disp(['    Final alignment1=' num2str(alignment1)]);
