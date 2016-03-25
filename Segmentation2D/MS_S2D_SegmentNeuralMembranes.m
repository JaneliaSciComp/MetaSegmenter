%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

% Identify regions corresponding to neural cells
function Ibwn = MS_S2D_SegmentNeuralMembranes(inputName,fracBlack,fracBlack2,varargin)   
    warning off;
    try
        p = MS_S2D_InputParser;              % call constructor
        p.parse(fracBlack, fracBlack2, varargin{:});
        options = MS_S2D_ExtractOptions(p, fracBlack, fracBlack2);
    catch
        output_usage_message();
        return
    end

    % Input image data  
    if strcmp(class(inputName),'char') && exist(inputName, 'file') == 2   % input is image file
        Igr = imread(inputName); 
    else
        Igr = inputName;        % input is image
    end

    disp(['class(Igr)=' class(Igr)]);
    if length(size(Igr) > 2)
        Igr = Igr(:,:,1);
    end

    subWin = [1:size(Igr,1); 1:size(Igr,2)];
    if numel(options.subWin) > 0
        limits = [ cellfun(@str2num,strsplit(options.subWin,','))];
        limits
        subWin = [limits(2):limits(4); limits(1):limits(3)];
        Igr = Igr(subWin(1,:),subWin(2,:)); 
    end

    disp(['In MS_S2D_SegmentNeuralMembranes: size(Igr)=' num2str(size(Igr))]);
    disp(['    options.fracBlack=' num2str(options.fracBlack)]);
    disp(['    max(Igr)=' num2str(max(max(Igr))) ' min(Igr)=' num2str(min(min(Igr)))]);
    if options.dispOn
        MS_S2D_ShowImage(Igr, 'Original image (Igr)', options);
    end

    if options.vesicles > 0
        radii_range = [10 10];
        disp(' ');
        disp(['Detecting vesicles in the range of radii [' num2str(radii_range) '] ...']);
        % Detect and remove vesicles
        if options.vesicles == 1
            [centers,radii] = imfindcircles(Igr, radii_range);
            disp(['numel(radii)='  num2str(numel(radii)) ' numel(centers)=' num2str(numel(centers))]);
            disp(['radii=' num2str(radii')]);
            for ic=1:numel(radii)
                Igr = rgb2gray(insertShape(Igr, ...
                               'circle', [centers(2*ic-1) centers(2*ic) radii(ic)], ...
                               'LineWidth', 5, 'Color', 'white'));
            end
            if options.dispOn
                MS_S2D_ShowImage(Igr, 'After deleting vesicles', options);
            end
        elseif options.vesicles == 2
            Ivesicles =  mat2gray(MS_UT_CircularHough(Igr, radii_range));
            if options.dispOn
                MS_S2D_ShowImage(Ivesicles, 'Detected vesicles ', options);
            end
        elseif options.vesicles == 3
            rthick = 3;
            quality = 1.0;
            flt2do = 0;
            for r=radii_range(1):radii_range(2)
                centers = MS_UT_HoughDetect(Igr, r, rthick, quality, flt2do); 
                if numel(centers) > 0
                    disp(['r=' num2str(r) ' num_centers=' num2str(size(centers,1))]);
                end
                for ic = 1:size(centers, 1)
                    Igr = rgb2gray(insertShape(Igr, ...
                               'circle', [centers(ic,1) centers(ic,2) r], ...
                               'LineWidth', 2, 'Color', 'white'));
                end
            end
            if options.dispOn
                MS_S2D_ShowImage(Igr, 'After deleting vesicles', options);
            end
        end
    end

    % If membrane probabilities are available,
    % use them to accent the membrane signals
    options.thr_min = 0;
    options.thr_max = 255;
    if length(options.membPr) > 0
        disp(' ');
        disp('Enhancing membrane signals ...');
        membPr = MS_S2D_ReadProbabilities(options.membPr, 2);
        membPr = membPr(subWin(1,:),subWin(2,:));
        disp(['min(membPr)=' num2str(min(min(membPr))) ' max(membPr)=' num2str(max(max(membPr)))]);
        if options.dispOn
            MS_S2D_ShowImage(mat2gray(membPr), 'Membrane probabilities', options);
        end
        % If only membrane probabilities are provided and mito probabilities are not,
        % use the membrane probabilities instead of original signals
        if length(options.mitoPr) > 0 && ~options.useMembPr
            Igr = enhance_membrane_signals_by_probabilities(Igr, membPr, options);
        else
            Igr = use_membrane_probabilities_instead_of_signals(Igr, membPr, options);
        end
        clear membPr;
        if options.dispOn && ~options.useMembPr
            MS_S2D_ShowImage(mat2gray(Igr), 'Enhanced original image (Igr)', options);
        end
    end

        % If mitochondria probabilities are available,
    % use them to whiten out mitochondria signals
    if length(options.mitoPr) > 0 
        disp(' ');
        disp('Whitening out mitochondria ...');
        mitoPr = MS_S2D_ReadProbabilities(options.mitoPr, 3);
        mitoPr = mitoPr(subWin(1,:),subWin(2,:));
        if options.dispOn || options.dispOn2
            MS_S2D_ShowImage(mat2gray(mitoPr), 'Mitochondria probabilities', options);
        end
        Igr = whiten_out_mitochondria(Igr, mitoPr, subWin, options);
%       Igr = whiten_out_mitochondria_new(Igr, mitoPr, options);
        disp(' ');
        if options.dispOn || options.dispOn2
            MS_S2D_ShowImage(Igr, 'Enhanced raw map, mitochondria whitened out (Igr)', options);
        end
        clear mitoPr;
        if length(options.membPr) > 0
            membPr = MS_S2D_ReadProbabilities(options.membPr, 2);
            membPr = membPr(subWin(1,:),subWin(2,:));
            Igrw = enhance_membrane_signals_by_probabilities(Igr, membPr, options);
        end
    end

    if length(options.mitoMembPr) > 0 &&  ~options.useMembPr
        mitoMembPr = MS_S2D_ReadProbabilities(options.mitoMembPr, 4);
        disp(['max(mitoMembPr)=' max(max(mitoMembPr))]);
        mitoMembPr = normalize_probabilities(mitoMembPr);
        mitoMembPr = mitoMembPr(subWin(1,:),subWin(2,:));
        if options.dispOn
            MS_S2D_ShowImage(mat2gray(mitoMembPr), 'Mitochondria probabilities', options);
        end
        if length(options.membPr) > 0
            membPr = MS_S2D_ReadProbabilities(options.membPr, 2);
            Igr(mitoMembPr > membPr) = 255;
        end
        if options.dispOn
            MS_S2D_ShowImage(Igr, 'Enhanced raw map, mito-membranes whitened out (Igr)', options);
        end
    end

    clear I;

    % Handle subsections
    im_size = size(Igr);
    Ibwn = zeros(im_size(1),im_size(2));

    [M_thr, M] = MS_S2D_GetThresholdIntensity(Igr, options);           

    disp(['max(M_thr)=' num2str(max(max(M_thr))) ' min(M_thr)=' num2str(min(min(M_thr))) ' max(Igr)=' num2str(max(max(Igr)))]);
    if max(max(Igr)) > 1
        Igr = Igr / 255;
    end
    disp(['size(Igr)=' num2str(size(Igr)) ' size(M_thr)=' num2str(size(M_thr))]);

    % If requested an output of thresholds, don't do anything else
    if length(options.outThr) > 0
        return;
    end

    % Output binary map for optimal threshold
    Ibwn = segment_neural_membranes(Igr, M_thr, options);
    if length(options.outBW) > 0
        imwrite(Ibwn, options.outBW);                               
    end
   
    % Show a heat map of the thresho;lds matrix
    disp(['options.heatmap=' num2str(options.heatmap)]);
    if options.heatmap 
        MS_S2D_ShowImage(M_thr, 'Heat map of the thresholds matrix', options);
    end
 
    % Output SEG and RGB files
    if (options.padding || options.RGB || ...
        options.dispOn  || options.dispOn2)
        Ibwn = MS_S2D_AddBoundaryPadding(Ibwn, 0);
    end
    if options.RGB || length(options.outSeg) > 0 
        L = MS_S2D_GenerateLabelsMatrix(Ibwn, options.verbose);
        max_label = max(max(L));
        disp(' ');
        disp(['num_regions=' num2str(max_label)]);
        if length(options.outSeg) > 0
            disp(['size(L)=' num2str(size(L))]);
            disp(['Saving the labels matrix in file ' options.outSeg ]);
            imwrite(L, options.outSeg);
        end

        % Color components
        if options.RGB
            Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
            MS_S2D_ShowImage(Lrgb, 'Colored labels of neurons (Lrgb)', options);
        end
    end

% -----------------------------------------------------------------------------

function Ibwn = segment_neural_membranes(Igr, M_thr, options)
    % Black-white image
    %
    % Determine intensity threshold for conversion to black-white image
    % NOTE: imhist assumes that I is a grayscale image, with range of values [0, 1]
    % NOTE2: edges are computed incorrectly by Matlab's imhist function,
    %        so a custom function has been implemented for compuring the edges
    %
    disp(' ');
    disp(['In segment_neural_membranes max(max(Igr))=' num2str(max(max(Igr)))]);
    Ibw = im2bw(Igr, 1);     % all values are == logical(0)
    imwrite(Igr, 'Igr1_MS_S2D_SegmentNeuralMembranes.png');
    imwrite(Ibw, 'Ibw_MS_S2D_SegmentNeuralMembranes.tiff');
    disp(['size(Igr)=' num2str(size(Igr)) ' size(M_thr)=' num2str(size(M_thr))]);
    if size(Igr) == size(M_thr)
        Ibw(Igr > M_thr/255.) = logical(1);
    end

    % removing mito boundaries % removing mito boundariesnd
    clear Igr;

    if options.dispOn
        MS_S2D_ShowImage(Ibw, 'Black-white image (Ibw)', options);
    end

     % Fill holes and open gaps
    Ibw  = MS_S2D_AddBoundaryPadding(Ibw, 0);
    Ibwf = imfill(Ibw, 'holes');
    clear Ibw;

    imwrite(Ibwf, 'Ibwf_MS_S2D_SegmentNeuralMembranes.png');
    Ibwn = bwareaopen(Ibwf, double(options.minSize));
    clear Ibwf;

    if options.dispOn | options.dispOn2
        MS_S2D_ShowImage(Ibwn, 'Cleaned BW segmentation (Ibwn)', options);
    end
    % Compute the distance transform
    D = bwdist(Ibwn);
%   D = bwdist(Ibwn, 'chessboard');
    clear Ibwn;
%   % Segmentation matrix
    L = watershed(D);
    Ibwn = ones(size(D));
    % BW image with thin boundaries: 
    Ibwn(L == 0)     = 0;
    Ibwn(M_thr == 0) = 1; % removing mito boundaries
    Ibwn = imfill(Ibwn, 'holes');

    if options.dispOn | options.dispOn2
        MS_S2D_ShowImage(Ibwn, 'Final segmentation of neural membranes (Ibwn)', options);
    end

% -----------------------------------------------------------------------------

function Iw = use_membrane_probabilities_instead_of_signals(I, membPr, options)
    Iw = round(membPr.*255.);
    Iw = imcomplement(mat2gray(Iw));
    if options.dispOn
        MS_S2D_ShowImage(mat2gray(Iw), 'Signals weighted by membrane probabilities (Iw)', options);
    end

% -----------------------------------------------------------------------------

function norm_probs = normalize_probabilities(probs)
    max_prob = max(max(probs));
    norm_probs = probs/max_prob;

% -----------------------------------------------------------------------------

function Iw = enhance_membrane_signals_by_probabilities(I, membPr, options)
    Iw   = mat2gray(round(double(I) - membPr .* double(I)));

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
        Ibw1 = bwareaopen(Ibw, 200);
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
            Ibw = bwareaopen(Ibw, 200);
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

function counts = get_counts2(Igr, max_intensity)
    counts = zeros(1, max_intensity + 1);
    for i=0:max_intensity
        counts(i+1) = sum(sum(round(Igr) == i));
    end

% -----------------------------------------------------------------------------

function Igrw = whiten_out_mitochondria(Igr1, mitoPr, subWin, options)
    % First, whiten out the mito regions in an 'old' way
    suppression_multiplier = 1.5; 
    suppression_factor = exp(-mitoPr .* suppression_multiplier);
    if length(options.membPr) > 0
        membPr = normalize_probabilities(MS_S2D_ReadProbabilities(options.membPr, 2));
        membPr = membPr(subWin(1,:),subWin(2,:));
        mitoPr(membPr > 0.5*mitoPr) = 0.;
        suppression_factor = exp(-mitoPr .* (1 - membPr) .* suppression_multiplier);
    end
    max_Igr = double(max(max(Igr1)));
    Igrw = mat2gray(round((double(Igr1).*suppression_factor + (1.-suppression_factor) * max_Igr)/max_Igr * 255.));
    if options.dispOn
        MS_S2D_ShowImage(Igrw, 'After the 1st step of whitening out', options);
    end

    % Second, detect an whiten out the entire mito regions 
    % where fraction of mito pixels >= maxMitoFr
    if 0
        Igr1 = Igrw;
        counts = get_counts2(round(Igr1*255), 255);
        [thr_Otsu,mu]  = MS_S2D_GetOtsuThresholds(counts, 1);
        Ibw = im2bw(Igr1, 1);
        Ibw(Igr1 > thr_Otsu/255.) = logical(1);
        Ibw = bwareaopen(Ibw, double(options.minSize));
        Ibw = imfill(Ibw, 'holes');
        if options.dispOn
            MS_S2D_ShowImage(Ibw, 'Binary map at Otsu threshold ', options);
        end
        bgSignal = mean(mean(Igr1(Ibw == 1)));
        % Detect mito components, dilate them and 
        % whiten out tyhe grayscale image at their location
        mitoPrBg = compute_background_mitoPr(mitoPr, options);
        CC = bwconncomp(Ibw);
        num_mito_regions = 0;
        for k = 1:CC.NumObjects
            if is_mito(Ibw, CC.PixelIdxList{k}, mitoPr, mitoPrBg) >= options.maxMitoFr
                num_mito_regions = num_mito_regions + 1;
                Ibw1 = im2bw(Igr1, 1);
                Ibw1(CC.PixelIdxList{k}) = logical(1);
%               for i = 1:9 % 9 pixels is the expected thickness of a mitochondrial membrane at Otsu threshold
                    Ibw1 = imdilate(Ibw1,strel('disk', 1));
%               end
                Igr1(Ibw1 > 0) = bgSignal;
            end
        end
        disp([' Number of mito regions with mito_frac > maxMitoFr (=' num2str(options.maxMitoFr) '): ', num2str(num_mito_regions)]);
        if options.dispOn
            MS_S2D_ShowImage(Igr1, 'After the 2nd step of whitening out', options);
        end

        Igrw = Igr1;
    end

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_SegmentNeuralMembranes(inputName,fracBlack,[,parameters])');
    disp('Required arguments:');
    disp('    inputName    - name of input data (file or image)');
    disp('    fracBlack    - desired fraction of black for black/wite conversion');
    disp('Optional parameters, specified as parameter-value pairs:');
    disp('    dispOn       - weather or not to display covariance matrix (default=0)');
    disp('    closeAll     - close previous image when displayong next (default=1)');
    disp('    nx, ny       - numbers of section to which subdivide the image (defaults=1)');
    disp('    ix, iy       - ids of the section to be processed/shown (defaults=1)');
    disp('    membPr       - path to HDF5 file containing membrane probabilities');
    disp('    mitoPr       - path to HDF5 file containing mitochondria probabilities');
    disp('    mitoMembPr   - path to HDF5 file containing mitochondrial membrane probabilities');
    disp('    useMembPr    - whether or not to use membrane probabilities instead of original signals');
    disp('    outBW        - name of output black/white image file (default='', no output)');
    disp('    outSeg       - name of output segmentation file (default='', no output)');
    disp('    RGB          - whether or not to show RGB labels (default = no)');
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return;

