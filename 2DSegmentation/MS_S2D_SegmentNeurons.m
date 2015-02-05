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

        compiled_code = 0;
        if isnumeric(fracBlack)     % original Matlab code
            options.fracBlack = p.Results.fracBlack;
            options.verbose   = p.Results.verbose;
            options.dispOn    = p.Results.dispOn;
            options.dispOn2   = p.Results.dispOn2;
            options.closeAll  = p.Results.closeAll;
            options.nx        = p.Results.nx;        
            options.ny        = p.Results.ny;        
            options.ix        = p.Results.ix;        
            options.iy        = p.Results.iy;       
            options.outName   = p.Results.outName; 
        else
            compiled_code     = 1;    % compiled code
            options.fracBlack =       str2double(p.Results.fracBlack);
            options.verbose   = int32(str2double(p.Results.verbose));
            options.dispOn    = int32(str2double(p.Results.dispOn));
            options.dispOn2   = int32(str2double(p.Results.dispOn2));
            options.closeAll  = int32(str2double(p.Results.closeAll));
            options.nx        = int32(str2double(p.Results.nx));
            options.ny        = int32(str2double(p.Results.ny));
            options.ix        = int32(str2double(p.Results.ix));
            options.iy        = int32(str2double(p.Results.iy));    
            options.outName   =                  p.Results.outName;
        end
    catch
        output_usage_message();
        return
    end

    disp(['In MS_S2D_SegmentNeurons: options.dispOn2=' num2str(options.dispOn2)]);
    frac_black = 0.34;
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
    im_size = size(I);
    [xpixels, ypixels] = get_pixels_to_show(im_size, options);
    Igr = mat2gray(I(xpixels,ypixels,1));
    In  = round(double(Igr)*255.0/max(max(double(Igr))));
    if options.verbose
        disp(['min(I)=' num2str(min(min(Igr))) ' max(Igr)=' num2str(max(max(Igr))) ...
              ' mean2(Igr)=' num2str(mean2(Igr)) ' std2(Igr)=' num2str(std2(Igr))]);
        disp(['min(In)=' num2str(min(min(In))) ' max(In)=' num2str(max(max(In))) ...
              ' mean2(In)=' num2str(mean2(In)) ' std2(In)=' num2str(std2(In))]);
    end
    if options.dispOn
        imshow(Igr);
        title('Original image (Igr)');
        waitforbuttonpress;
        if options.closeAll
            close all;
        end
    end

    % Black-white image
    %
    % Determine intensity threshold for conversion to black-white image
    % NOTE: imhist assumes that I is a grayscale image, with range of values [0, 1]
    % NOTE2: edges are computed incorrectly by Matlab's imhist function,
    %        so a custom function has been implemented for compuring the edges
    threshold = MS_S2D_GetThresholdIntensity(Igr, options.verbose, frac_black, 1001);
    Ibw = im2bw(Igr, threshold);
    Ibw = bwareaopen(Ibw,20);
    if options.verbose
        disp(['final threshold_intensity=' num2str(threshold)]);
        disp(['size(Igr)=' num2str(size(Igr)) ' class(Igr)=' class(Igr(1,1))]);
        disp(['size(Ibw)=' num2str(size(Ibw)) ' class(Ibw)=' class(Ibw(1,1))]);
        disp(['my frac_black=' num2str(sum(sum(Ibw == 0))/numel(Ibw))]);
    end
    if options.dispOn
        figure
        imshow(Ibw), title('Black-white image (Ibw)')
        drawnow;
        waitforbuttonpress;
        if options.closeAll
            close all;
        end
    end

    % Fill holes and open gaps
    Ibw = adjust_edges(Ibw, 1);
    Ibwf = imfill(Ibw, 'holes');
    Ibwn = bwareaopen(Ibwf,20);
    if options.dispOn
        imshow(Ibwn)
        drawnow;
        waitforbuttonpress;
        if options.closeAll
            close all;
        end
    end

    % Optionally output results to file
    if length(options.outName) > 0
        imwrite(Ibwn, options.outName);
    end

    % Display labels
    disp(['In MS_S2D_SegmentNeurons: options.dispOn2=' num2str(options.dispOn2)]);
    if (options.dispOn || options.dispOn2) 
        % Label components in the inverse BW image
        L = generate_labels_matrix(Ibwn, options.verbose);

        % Color components
        if options.dispOn
            Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
            figure
            imshow(Lrgb)
            title('Colored watershed label matrix (Lrgb)')
            drawnow;
            waitforbuttonpress;
            if options.closeAll
                close all;
            end
        end
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
    disp('    outName      - name of output file (default='', no output)');
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return;

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

% -----------------------------------------------------------------------------

% Add zero padding at the adges, 
% to ensure thgat that entire image != one BW component
function Ibw = adjust_edges(Ibw, value)
    size1 = size(Ibw, 1);
    size2 = size(Ibw, 2);

    i = 1;
    while i <= size1 && sum(Ibw(i,:)) >= size2*0.9
        Ibw(i,:) = logical(value);
        i = i+1;
    end

    i = 1;
    while i <= size2 && sum(Ibw(:,i)) >= size1*0.9
        Ibw(:,i) = logical(value);
        i = i+1;
    end

    i = 0;
    while size1 - i > 1 && sum(Ibw(size1 - i,:)) >= size2*0.9
        Ibw(size1 - i,:) = logical(value);
        i = i+1;
    end

    i = 0;
    while size2 - i > 1 && sum(Ibw(:,size2 - i)) >= size1*0.9
        Ibw(:,size2 - i) = logical(0);
        i = i+1;
    end

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

