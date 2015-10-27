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
    disp(['options.fracBlack=' num2str(options.fracBlack)]);
    disp(['max(Igr)=' num2str(max(max(Igr))) ' min(Igr)=' num2str(min(min(Igr)))]);
    disp(['max(Igr)=' num2str(max(max(Igr))) ' min(Igr)=' num2str(min(min(Igr)))]);
    if options.dispOn
        MS_S2D_ShowImage(Igr, 'Original image (Igr)', options);
    end

    % If membrane probabilities are available,
    % use them to accent the membrane signals
    if length(options.membPr) > 0
        membPr = MS_S2D_ReadProbabilities(options.membPr);
        if options.dispOn
            MS_S2D_ShowImage(mat2gray(membPr), 'Membrane probabilities', options);
        end
        Igr = weight_membrane_signals_by_probabilities(Igr, membPr, options);
        clear membPr;
        if options.dispOn
            MS_S2D_ShowImage(mat2gray(Igr), 'Weighted original image (Igr)', options);
        end
    end

        % If mitochondria probabilities are available,
    % use them to whiten out mitochondria signals
    if length(options.mitoPr) > 0
        mitoPr = MS_S2D_ReadProbabilities(options.mitoPr);
        if options.dispOn
            MS_S2D_ShowImage(mat2gray(mitoPr), 'Mitochondria probabilities', options);
        end
        Igr = whiten_out_mitochondria(Igr, mitoPr, options);
        clear mitoPr;
        if options.dispOn
            MS_S2D_ShowImage(Igr, 'Enhanced raw map, mitochondria whitened out (Igr)', options);
        end
    end

    clear I;

    im_size = size(Igr);
    Ibwn = zeros(im_size(1),im_size(2));
    subsections = MS_S2D_DivideImageIntoSubsections(im_size, options);
    [M_thr, M_thr2] = MS_S2D_GetThresholdIntensity(Igr, 1, subsections, options);           
    clear M_thr2;
    disp(['M_thr=' num2str(max(max(M_thr)))]);
    Ibwn = segment_neural_membranes(Igr, M_thr, options);
    if length(options.outBW) > 0
        imwrite(Ibwn, 'Ibwn_MS_S2D_SegmentNeuralMembranes.png');
    end
    if (options.padding || length(options.outRGB) > 0 || ...
        options.dispOn  || options.dispOn2)
        Ibwn = MS_S2D_AddBoundaryPadding(Ibwn, 0);
    end

    % Optionally output BW file
    if length(options.outBW) > 0
        imwrite(Ibwn, options.outBW);
    end

     % Display labels
    if (options.dispOn || options.dispOn2 || length(options.outSeg) > 0 ...
                       || length(options.outRGB) > 0)
        L = MS_S2D_GenerateLabelsMatrix(Ibwn, options.verbose);
        if length(options.outSeg) > 0
            imwrite(L, options.outSeg);
        end

        % Color components
        if options.dispOn || options.dispOn2 || length(options.outRGB) > 0
            Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
            MS_S2D_ShowImage(Lrgb, 'Colored labels of neurons (Lrgb)', options);
            if length(options.outRGB) > 0
                imwrite(Lrgb, options.outRGB);
            end
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
    Ibw = im2bw(Igr, 1);     % all values are == logical(0)
    imwrite(Igr, 'Igr1_MS_S2D_SegmentNeuralMembranes.png');
    imwrite(Ibw, 'Ibw_MS_S2D_SegmentNeuralMembranes.tiff');
    disp(['size(Igr)=' num2str(size(Igr)) ' size(M_thr)=' num2str(size(M_thr))]);
    Ibw(Igr > M_thr) = logical(1);
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
    % Segmentation matrix
    L = watershed(D);
    Ibwn = ones(size(D));
    % BW image with thin boundaries: 
    Ibwn(L == 0) = 0;
    if options.dispOn | options.dispOn2
        MS_S2D_ShowImage(Ibwn, 'Final segmentation of neural membranes (Ibwn)', options);
    end


% -----------------------------------------------------------------------------

function Iw = use_probabilities_instead_of_signals(I, membPr, options)
    Iw = round(membPr.*255.);
    if options.dispOn
        MS_S2D_ShowImage(mat2gray(Iw), 'Signals weighted by membrane probabilities (Iw)', options);
    end
    Iw = imcomplement(mat2gray(Iw));

% -----------------------------------------------------------------------------

function Iw = weight_membrane_signals_by_probabilities(I, membPr, options)
    weight_multiplier = 4;
    weighting_factor = 1. - exp(-membPr.*weight_multiplier);
    Imin = min(min(I));
    Iw   = mat2gray(round(double(I) - membPr .* double(I)));

% -----------------------------------------------------------------------------

function Igrw = whiten_out_mitochondria(Igr, mitoPr, options)
    suppression_multiplier = 1.;  
    suppression_factor = exp(-mitoPr .* suppression_multiplier);
    if length(options.membPr) > 0
        membPr = MS_S2D_ReadProbabilities(options.membPr);
%       mitoPr(membPr > mitoPr) = 0.;
        suppression_factor = exp(-mitoPr .* (1 - membPr) .* suppression_multiplier);
    end
    % Suppress motochondria signals
    disp(['max(suppression_factor)=' num2str(max(max(suppression_factor))) ...
         ' min(suppression_factor)=' num2str(min(min(suppression_factor)))]);
    max_Igr = double(max(max(Igr)));
    Igrw = mat2gray(round((double(Igr).*suppression_factor + (1.-suppression_factor) * max_Igr)/max_Igr * 255.));
    clear Igr;

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
    disp('    outBW        - name of output black/white image file (default='', no output)');
    disp('    outSeg       - name of output segmentation file (default='', no output)');
    disp('    outRGB       - name of output colored labels file (default='', no output)');
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return;

