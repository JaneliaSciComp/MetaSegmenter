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
    if length(options.membPr) > 0
        disp(' ');
        disp('Enhancing membrane signals ...');
        membPr = MS_S2D_ReadProbabilities(options.membPr, 2);
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
        if options.dispOn
            MS_S2D_ShowImage(mat2gray(mitoPr), 'Mitochondria probabilities', options);
        end
        Igr = whiten_out_mitochondria(Igr, mitoPr, options);
        disp(' ');
        disp(['distType=' num2str(options.distType)]);
        if options.dispOn
            MS_S2D_ShowImage(Igr, 'Enhanced raw map, mitochondria whitened out (Igr)', options);
        end
        clear mitoPr;
        if length(options.membPr) > 0
            membPr = MS_S2D_ReadProbabilities(options.membPr, 2);
            Igrw = enhance_membrane_signals_by_probabilities(Igr, membPr, options);
        end
    end

    if length(options.mitoMembPr) > 0 &&  ~options.useMembPr
        mitoMembPr = MS_S2D_ReadProbabilities(options.mitoMembPr, 4);
        disp(['max(mitoMembPr)=' max(max(mitoMembPr))]);
        mitoMembPr = normalize_probabilities(mitoMembPr);
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
    subsections = MS_S2D_DivideImageIntoSubsections(im_size, options);
    [M_thr1, M_thr2, M_thr] = MS_S2D_GetThresholdIntensity(Igr, subsections, options);           
    disp(['max(M_thr)=' num2str(max(max(M_thr))) ' max(Igr)=' num2str(max(max(Igr)))]);
    if max(max(Igr)) > 1
        Igr = Igr / 255;
    end
    disp(['size(Igr)=' num2str(size(Igr)) ' size(M_thr)=' num2str(size(M_thr))]);

    % If requested an output of thresholds, don't do anything else
    if length(options.outThr) > 0
        return;
    end

    % Display binary map for upper threshold
    if options.dispOn
        Ibw1 = im2bw(Igr, 1);     % all values are == logical(0)
        if size(Igr) == size(M_thr1)
            Ibw1(Igr > M_thr1/255.) = logical(1);
        end
        MS_S2D_ShowImage(Ibw1, 'Binary map at upper threshold ', options);
    end
    clear M_thr1;

    % Display binary map for lower threshold
    if options.dispOn
        Ibw2 = im2bw(Igr, 1);     % all values are == logical(0)
        if size(Igr) == size(M_thr2)
            Ibw2(Igr > M_thr2/255.) = logical(1);
        end
        MS_S2D_ShowImage(Ibw2, 'Binary map at lower threshold ', options);
    end
    clear M_thr2;

    % Display binary map for optimal threshold
    Ibwn = segment_neural_membranes(Igr, M_thr, options);
    if length(options.outBW) > 0
        imwrite(Ibwn, 'Ibwn_MS_S2D_SegmentNeuralMembranes.png');
    end
    if (options.padding || options.RGB || ...
        options.dispOn  || options.dispOn2)
        Ibwn = MS_S2D_AddBoundaryPadding(Ibwn, 0);
    end

    % Optionally output BW file
    if length(options.outBW) > 0
        imwrite(Ibwn, options.outBW);
    end

     % Display labels
    disp(['options.RGB=' num2str(options.RGB)]);
    if options.RGB || length(options.outSeg) > 0 
        L = MS_S2D_GenerateLabelsMatrix(Ibwn, options.verbose);
        max_label = max(max(L));
        disp(' ');
        disp(['num_regions=' num2str(max_label)]);
        if length(options.outSeg) > 0
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
    clear Igr;

    Ibw = bwareaopen(Ibw,20);
    imwrite(Ibw, 'Ibw1_MS_S2D_SegmentNeuralMembranes.tiff');
    if options.dispOn
        MS_S2D_ShowImage(Ibw, 'Black-white image (Ibw)', options);
    end

     % Fill holes and open gaps
    Ibw  = MS_S2D_AddBoundaryPadding(Ibw, 0);
    Ibwf = imfill(Ibw, 'holes');
    clear Ibw;

    imwrite(Ibwf, 'Ibwf_MS_S2D_SegmentNeuralMembranes.tiff');
    Ibwn = bwareaopen(Ibwf,200);
%   Ibwn = imerode(Ibwn, strel('disk', 1));
    clear Ibwf;

    if options.dispOn | options.dispOn2
        MS_S2D_ShowImage(Ibwn, 'Cleaned BW segmentation (Ibwn)', options);
    end
    % Compute the distance transform
    D = bwdist(Ibwn);
    clear Ibwn;
%   % Segmentation matrix
    L = watershed(D);
    Ibwn = ones(size(D));
    % BW image with thin boundaries: 
    Ibwn(L == 0) = 0;
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
    weight_multiplier = 4;
    weighting_factor = 1. - exp(-membPr.*weight_multiplier);
    Imin = min(min(I));
    Iw   = mat2gray(round(double(I) - membPr .* double(I)));

% -----------------------------------------------------------------------------

function Igrw = whiten_out_mitochondria(Igr, mitoPr, options)
    suppression_multiplier = 1.0;   
    suppression_factor = exp(-mitoPr .* suppression_multiplier);
    if length(options.membPr) > 0
        membPr = normalize_probabilities(MS_S2D_ReadProbabilities(options.membPr, 2));
        mitoPr(membPr > 0.3*mitoPr) = 0.;
        suppression_factor = exp(-mitoPr .* (1 - membPr) .* suppression_multiplier);
    end
    % Suppress motochondria signals
    disp(['max(suppression_factor)=' num2str(max(max(suppression_factor))) ...
         ' min(suppression_factor)=' num2str(min(min(suppression_factor)))]);
    max_Igr = double(max(max(Igr)));
    Igrw = mat2gray(round((double(Igr).*suppression_factor + (1.-suppression_factor) * max_Igr)/max_Igr * 255.));
    clear Igr;
%   if length(options.membPr) > 0
%       membPr = normalize_probabilities(MS_S2D_ReadProbabilities(options.membPr));
%       Igrw = enhance_membrane_signals_by_probabilities(Igrw, membPr, options);
%   end

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

