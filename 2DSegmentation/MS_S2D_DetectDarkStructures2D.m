%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: detect the location of "dark" structures, 
% such as mitochondria, T-,E-bars, and glial cells 
% in a 2D intensity image  
% 

function Ibwd = MS_S2D_DetectDarkStructures2D(inputName,fracBlack,varargin)   
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

    disp(['In MS_S2D_DetectDarkStructures2D: options.dispOn2=' num2str(options.dispOn2)]); 
    frac_black = 0.6;
    if options.fracBlack > 0
        frac_black = options.fracBlack;
    end

    % Input image data
    if strcmp(class(inputName),'char') && exist(inputName, 'file') == 2   % input is image file
        I = imread(inputName);
    else
        I = inputName;                 % input is image
    end

    % Original image
    if options.dispOn
        MS_S2D_ShowImage(Igr, 'Original image (Igr)', options);
    end

    if length(size(I) > 2)
        I = I(:,:,1);
    end
    im_size = size(I);
    Ibwd = zeros(im_size(1),im_size(2));
    subsections = MS_S2D_DivideImageIntoSubsections(im_size, options);
    Igr = mat2gray(I);
    M_thr = MS_S2D_GetThresholdIntensity(Igr, fracBlack, subsections, options);      
    Ibwd = segment_dark_structures(Igr, M_thr, options);
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

function Ibwd = segment_dark_structures(Igr, M_thr, options)
    % Black-white image
    %
    % Determine intensity threshold for conversion to black-white image
    % NOTE: imhist assumes that I is a grayscale image, with range of values [0, 1]
    % NOTE2: edges are computed incorrectly by Matlab's imhist function,
    %        so a custom function has been implemented for compuring the edges
    Ibw = im2bw(Igr, 1);
    Ibw(Igr > M_thr) = logical(1);
    Ibw = bwareaopen(Ibw,20);
    Ibw = MS_S2D_AddBoundaryPadding(Ibw, 1);
    % Ibw  = imfill(Ibw, 'holes');
    if options.dispOn
        MS_S2D_ShowImage(Ibw, 'Black-white image without holes (Ibw)', options);
    end

    % Handling complement
    Ibwc = imcomplement(Ibw);
    Ibwc = bwareaopen(Ibwc,200);
    Ibwc = imcomplement(Ibwc);
    if options.dispOn
        MS_S2D_ShowImage(Ibwc, 'Complement of dilated black-white image (Ibwc)', options);
    end

    % Dilated BW image
    Ibwcd = imdilate(Ibwc, strel('disk', 2));
    if options.dispOn
        MS_S2D_ShowImage(Ibwcd, 'Dilated, inversed, filled and eroded black-white image (Ibwd)', options);
    end

    % Inversed 
    Ibwcdc = imcomplement(Ibwcd);
    Ibwcdc = imerode(Ibwcdc, strel('disk', 2));
    Ibwcdc = bwareaopen(Ibwcdc,200);
    Ibwcdc = imfill(Ibwcdc, 'holes');
    Ibwcdc = imdilate(Ibwcdc, strel('disk', 2));
    Ibwd  = imfill(Ibwcdc, 'holes');
%   Ibwd   = MS_S2D_AddBoundaryPadding(Ibwd, 0);
    if options.dispOn
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
    disp('    outBW        - name of output black/white image file (default='', no output)');
    disp('    outSeg       - name of output segmentation file (default='', no output)');
    disp('    outRGB       - name of output colored labels file (default='', no output)');
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return;

% -----------------------------------------------------------------------------

function [xpixels, ypixels] = get_pixels_to_show(im_size, options)
    xlen = round(im_size(1)/options.nx);
    ylen = round(im_size(2)/options.ny);
    xbeg = (options.ix -1)*xlen + 1;
    xend =  options.ix    *xlen;
    ybeg = (options.iy -1)*ylen + 1;
    yend =  options.iy    *ylen;
    if options.verbose
        disp(['xbeg=' num2str(xbeg) ' xend=' num2str(xend) ...
             ' ybeg=' num2str(ybeg) ' yend=' num2str(yend)]);
    end
    xpixels = [xbeg:xend];
    ypixels = [ybeg:yend];

% -----------------------------------------------------------------------------

function L = generate_labels_matrix(Ibwf, verbose)
    CC = bwconncomp(Ibwf); % NOTE: this function does not label holes
    Nc = CC.NumObjects;
    imsize = CC.ImageSize;

    if verbose
        disp(['imsize=' num2str(imsize) ' num_comp=' num2str(Nc) ]);
    end
    L = zeros(imsize);
    for k = 1:Nc
        for m=1:numel(CC.PixelIdxList{k})
            lin_ind = CC.PixelIdxList{k}(m);
            [r, c] = ind2sub(CC.ImageSize, lin_ind);
%           disp(['lin_ind=' num2str(lin_ind) ' r=' num2str(r) ' c=' num2str(c)]);
            L(r,c) = k;
        end
    end

