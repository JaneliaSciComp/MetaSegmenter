%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

% Identify regions corresponding to neural cells
function Ibwn = MS_S2D_SegmentNeurons(inputName,fracBlack,varargin)   
    warning off;
    clc;
    try
        p = MS_S2D_InputParser;              % call constructor
        p.parse(fracBlack, varargin{:});
        options = MS_S2D_ExtractOptions(p, fracBlack);
    catch
        output_usage_message();
        return
    end

    disp(['In MS_S2D_SegmentNeurons: options.dispOn2=' num2str(options.dispOn2)]);
    
    % Input image data  
    if strcmp(class(inputName),'char') && exist(inputName, 'file') == 2   % input is image file
        I = imread(inputName); 
    else
        I = inputName;                 % input is image
    end

    size(I)
    if length(size(I) > 2)
        I = I(:,:,1);
    end
    Igr = mat2gray(I);
    % Original image
    if options.dispOn
        MS_S2D_ShowImage(Igr, 'Original image (Igr)', options);
    end

    im_size = size(I);
    Ibwn = zeros(im_size(1),im_size(2));
    subsections = MS_S2D_DivideImageIntoSubsections(im_size, options);
    M_thr = MS_S2D_GetThresholdIntensity(Igr, fracBlack, subsections, options);           
    Ibwn = segment_neurons(Igr, M_thr, options);
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
            MS_S2D_ShowImage(Lrgb, 'Colored labels of neurons (Lrgb)', options)
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
    Ibw = im2bw(Igr, 1);
    mean(mean(M_thr))
    Ibw(Igr > M_thr) = logical(1);
    Ibw = bwareaopen(Ibw,20);
    if options.dispOn
        MS_S2D_ShowImage(Ibw, 'Black-white image (Ibw)', options);
    end

    % Fill holes and open gaps
    Ibw  = MS_S2D_AddBoundaryPadding(Ibw, 0);
    Ibwf = imfill(Ibw, 'holes');
    Ibwn = bwareaopen(Ibwf,20);
%   Ibwn = MS_S2D_AddBoundaryPadding(Ibwn, 0);
    if options.dispOn
        MS_S2D_ShowImage(Ibwn, 'Black-white image (Ibwn)', options);
    end

    % Optionally output results to file
    if length(options.outBW) > 0
        imwrite(Ibwn, options.outBW);
    end

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_SegmentNeurons(inputName,fracBlack,[,parameters])');
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

