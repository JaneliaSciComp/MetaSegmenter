%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: detect the location of "dark" structures, 
% such as mitochondria, T-,E-bars, and glial cells 
% in a 2D intensity image  
% 

function Ibwd = MS_S2D_DetectDarkStructures2D(inputName,fracBlack,fracBlack2,varargin)   
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
        I = inputName;                 % input is image
    end

    imwrite(I, 'orig_image_MS_S2D_DetectDarkStructures2D.png')

    % Original image
    if options.dispOn
        MS_S2D_ShowImage(I, 'Original image (I)', options);
    end

    if length(size(I) > 2)
        I = I(:,:,1);
    end
    im_size = size(I);
    Ibwd = zeros(im_size(1),im_size(2));
    subsections = MS_S2D_DivideImageIntoSubsections(im_size, options);
    Igr = mat2gray(I);
    clear I;

    [M_thr, M_thr2] = MS_S2D_GetThresholdIntensity(Igr, 2, subsections, options);      
    Ibwd = segment_dark_structures(Igr, M_thr, M_thr2, options);
    Ibwd = MS_S2D_AddBoundaryPadding(Ibwd, 0);

    if options.dispOn
        MS_S2D_ShowImage(Ibwd, 'Final (dilated, inversed, filled and eroded) black-white image (Ibwd)', options);
    end

    % Display/outpu labels
    if (options.dispOn || options.dispOn2 || length(options.outSeg) > 0 ...
                       || length(options.outRGB) > 0)
        % Label components in the inverse BW image
        L = MS_S2D_GenerateLabelsMatrix(Ibwd, options.verbose);
        if length(options.outSeg) > 0
            imwrite(L, options.outSeg);
        end

        % Color components
        if options.dispOn2 || options.outRGB
            Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
            MS_S2D_ShowImage(Lrgb, 'Colored labels of dark structures (Lrgb)', options);
            if length(options.outRGB) > 0
                imwrite(Lrgb, options.outRGB);
            end
        end
    end

% -----------------------------------------------------------------------------

function Ibwd = segment_dark_structures(Igr, M_thr, M_thr2, options)
    % Black-white image
    %
    % Determine intensity threshold for conversion to black-white image
    % NOTE: imhist assumes that I is a grayscale image, with range of values [0, 1]
    % NOTE2: edges are computed incorrectly by Matlab's imhist function,
    %        so a custom function has been implemented for compuring the edges
    Ibw = im2bw(Igr, 1);
    Ibw(Igr > M_thr)  = logical(1);
    Ibw2 = im2bw(Igr, 1);
    Ibw2(Igr < M_thr2) = logical(1);
    clear Igr;
    Ibw(Ibw2) = logical(0); % subtract Ibw2 from Ibw
    clear Ibw2; 

    orig_imsize = size(Ibw);
    scale = options.resize;
    Ibw = imresize(Ibw, scale);
    Ibw = bwareaopen(Ibw,20);
    Ibw = MS_S2D_AddBoundaryPadding(Ibw, 1);
    % Ibw  = imfill(Ibw, 'holes');
    if options.dispOn
        MS_S2D_ShowImage(Ibw, 'Black-white image without holes (Ibw)', options);
    end

    % Handling complement
    Ibwc = imcomplement(Ibw);
    clear Ibw;
    Ibwc = imerode(Ibwc, strel('disk', 2));
    Ibwc = bwareaopen(Ibwc,round(300*scale));
    if options.dispOn | options.dispOn2
        MS_S2D_ShowImage(Ibwc, 'Complement of dilated black-white image (Ibwc)', options);
    end

    % Dilated BW image
    Ibwc = imcomplement(Ibwc);
    Ibwcd = imdilate(Ibwc, strel('disk', 2));
    clear Ibwc;
    if options.dispOn
        MS_S2D_ShowImage(Ibwcd, 'Dilated, inversed, filled and eroded black-white image (Ibwd)', options);
    end

    % Inversed 
    Ibwcdc = imcomplement(Ibwcd);
    clear Ibwcd;
    Ibwcdc = imerode(Ibwcdc, strel('disk', 4));
    Ibwcdc = bwareaopen(Ibwcdc,round(300*scale));
    Ibwcdc = imfill(Ibwcdc, 'holes');
    Ibwcdc = imdilate(Ibwcdc, strel('disk', 3));
    Ibwd  = imfill(Ibwcdc, 'holes');
    clear Ibwcdc;

    Ibwd = imresize(Ibwd, orig_imsize);
%   Ibwd   = MS_S2D_AddBoundaryPadding(Ibwd, 0);
    if options.dispOn | options.dispOn2
        MS_S2D_ShowImage(Ibwd, 'Final (dilated, inversed, filled and eroded) black-white image (Ibwd)', options);
    end

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_DetectDarkStructures2D(image_file,fraction_of_black[,parameters])');
    disp('Required arguments:');
    disp('    inputName    - name of input data (file or image)');
    disp('    fracBlack    - desired fraction of black for black/wite conversion');
    disp('Optional parameters, specified as parameter-value pairs:');
    disp('    dispOn       - weather or not to display covariance matrix (default=0)');
    disp('    closeAll     - close previous image when displayong next (default=1)');
    disp('    nx, ny       - numbers of section to which subdivide the image (defaults=1)');
    disp('    ix, iy       - ids of the section to be processed/shown (defaults=1)');
    disp('    resize       - scale to be used when resizing BW image (default = 1)');
    disp('    outBW        - name of output black/white image file (default='', no output)');
    disp('    outSeg       - name of output segmentation file (default='', no output)');
    disp('    outRGB       - name of output colored labels file (default='', no output)');
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return;

