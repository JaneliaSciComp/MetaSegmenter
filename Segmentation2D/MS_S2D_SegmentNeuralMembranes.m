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
        I = imread(inputName); 
    else
        I = inputName;        % input is image
    end

    if length(size(I) > 2)
        I = I(:,:,1);
    end
    disp(['In MS_S2D_SegmentNeuralMembranes: size(I)=' num2str(size(I))]);
    disp(['options.fracBlack=' num2str(options.fracBlack)]);
    disp(['max(I)=' num2str(max(max(I))) ' min(I)=' num2str(min(min(I)))]);
    Igr = mat2gray(I);
    disp(['max(Igr)=' num2str(max(max(Igr))) ' min(Igr)=' num2str(min(min(Igr)))]);
    if options.dispOn
        MS_S2D_ShowImage(Igr, 'Original image (Igr)', options);
    end

    if length(options.membPr) > 0
        membPr = MS_S2D_ReadProbabilities(options.membPr);
        if options.dispOn
            MS_S2D_ShowImage(mat2gray(membPr), 'Membrane probabilities', options);
        end
        weight_signals_by_probabilities = 1;
        if weight_signals_by_probabilities
            Igr = weight_image_by_probabilities(I, membPr, options);
        else 
            Igr = use_probabilities_instead_of_signals(I, membPr, options);
        end
        clear membPr;
        if options.dispOn
            MS_S2D_ShowImage(mat2gray(Igr), 'Weighted original image (Igr)', options);
        end
    end
    clear I;

    im_size = size(Igr);
    Ibwn = zeros(im_size(1),im_size(2));
    subsections = MS_S2D_DivideImageIntoSubsections(im_size, options);
    [M_thr, M_thr2] = MS_S2D_GetThresholdIntensity(Igr, 1, subsections, options);           
    clear M_thr2;
    disp(['M_thr=' num2str(max(max(M_thr)))]);
    Ibwn = segment_neurons(Igr, M_thr, options);
    imwrite(Ibwn, 'Ibwn_MS_S2D_SegmentNeuralMembranes.tiff');
    Ibwn = MS_S2D_AddBoundaryPadding(Ibwn, 0);

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

function Ibwn = segment_neurons(Igr, M_thr, options)
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

    orig_imsize = size(Ibw);
    scale = options.resize; 
    Ibw = imresize(Ibw, scale);
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
%   Ibwf = imerode(Ibwf, strel('disk', 1));
    Ibwn = bwareaopen(Ibwf,20);
    Ibwn = imerode(Ibwn, strel('disk', 1));
    clear Ibwf;

%   Ibwn = MS_S2D_AddBoundaryPadding(Ibwn, 0);
    if options.dispOn | options.dispOn2
        MS_S2D_ShowImage(Ibwn, 'Black-white image (Ibwn)', options);
    end

    Ibwn = imresize(Ibwn, orig_imsize);
    % Optionally output results to file
    if length(options.outBW) > 0
        imwrite(Ibwn, options.outBW);
    end
    Ibwn = bwareaopen(Ibwn,200);
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

function Iw = weight_image_by_probabilities(I, membPr, options)
    Iw = imcomplement(mat2gray(double(imcomplement(I)).*membPr));


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

