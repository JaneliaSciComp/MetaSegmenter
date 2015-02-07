%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

% Combine regions corresponding to neurons and dark structures
% in order to produce final segmentation
function Ibw = MS_S2D_Segmentation2D(inputName,fracBlack,fracBlack2,varargin)   
    debug = 0;
%   warning off;
    clc;
    try
        p = MS_S2D_InputParser;              % call constructor
        p.parse(fracBlack, varargin{:});

        compiled_code = 0;
        if isnumeric(fracBlack)     % original Matlab code
            options.fracBlack =         p.Results.fracBlack;
            options.verbose   = logical(p.Results.verbose);
            options.dispOn    = logical(p.Results.dispOn);
            options.dispOn2   = logical(p.Results.dispOn2);
            options.closeAll  = logical(p.Results.closeAll);
            options.nx        =   int32(p.Results.nx);
            options.ny        =   int32(p.Results.ny);
            options.ix        =   int32(p.Results.ix);
            options.iy        =   int32(p.Results.iy);
            options.outName   =         p.Results.outName; 
        else
            compiled_code     = 1;    % compiled code
            options.fracBlack =         str2double(p.Results.fracBlack);
            options.verbose   = logical(str2double(p.Results.verbose));
            options.dispOn    = logical(str2double(p.Results.dispOn));
            options.dispOn2   = logical(str2double(p.Results.dispOn2));
            options.closeAll  = logical(str2double(p.Results.closeAll));
            options.nx        =   int32(str2double(p.Results.nx));
            options.ny        =   int32(str2double(p.Results.ny));
            options.ix        =   int32(str2double(p.Results.ix));
            options.iy        =   int32(str2double(p.Results.iy));    
            options.outName   =                    p.Results.outName;
        end
    catch
        output_usage_message();
        return
    end

    % Input image data  
    if exist(inputName, 'file') == 2   % input is image file
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
    disp(['In MS_S2D_Segmentation2D: options.dispOn2=' num2str(options.dispOn2)]);
    
    if (options.dispOn || options.dispOn2 || debug)
        imshow(Igr);
        title('Original image (Igr)');
        waitforbuttonpress;
        if options.closeAll
            close all;
        end
    end

    % Black-white images
    %
    Ibwn = MS_S2D_SegmentNeurons(        I, fracBlack, ...
                                 'dispOn',  options.dispOn, ...
                                 'dispOn2', options.dispOn2, ...
                                 'closeAll',options.closeAll);
    Ibwd = MS_S2D_DetectDarkStructures2D(I, fracBlack2,...
                                 'dispOn',  options.dispOn, ...
                                 'dispOn2', options.dispOn2, ...
                                 'closeAll',options.closeAll);
    if debug
        figure
        imshow(Ibwn)
        title('Initially detectred neural cells (Ibwn)');
        drawnow;
        waitforbuttonpress;

        figure
        imshow(Ibwd)
        title('Initially detectred dark structures (Ibwd)');
        drawnow;
        waitforbuttonpress;
    end

    % #########
    % Step 1:
    % #########

    % Identify the dark strutures each of which is entirely contained 
    % within a single neural cell region.   
    Ibwnd = imdilate(Ibwn, strel('disk', 1));
    Ibwn1 = (Ibwnd | Ibwd) - Ibwd;       
    % What we are looking for, should be holes in the remaing BW components
    Ibwn2 = imfill(Ibwn1, 'holes');
    Ibwd1 = Ibwn2 - Ibwn1;
    Ibwd11 = bwareaopen(Ibwd1,200);

    % Identify the dark strutures each of which is ALMOST entirely contained
    % within a single neural cell region.
    Ibwde  = imerode(Ibwd - Ibwd11, strel('disk', 1));
    Ibwn1 = (Ibwnd | Ibwde) - Ibwde;
    Ibwn2 = imfill(Ibwn1, 'holes');
    Ibwd12 = Ibwn2 - Ibwn1;
    Ibwd12  = bwareaopen(Ibwd12,200);
    Ibwd1f = Ibwd11 | Ibwd12;
    if debug 
        figure
        imshow(Ibwd1f)
        title('Final dark structures which are entirely contained in neural cell regions (Ibwd1f)');
        drawnow;
        waitforbuttonpress;
    end

    % #########
    % Step 2:
    % #########

    % Identify the (parts of) dark strutures that are located
    % - either outside of (slightly dilated) large neural cell regions
    % - or entirely cover very small neural cell regions.
    Ibwnd = imdilate(Ibwn, strel('disk', 2));
    Ibwd1 = Ibwd - Ibwd1f;
    Ibwd2 = (Ibwnd | Ibwd1) - Ibwnd;
    Ibwd2 = bwareaopen(Ibwd2,400);
    Ibwd2f = imfill(Ibwd2, 'holes');
%   Ibwd1 = imclose(Ibwd1, strel('disk', 2));

    if debug
        figure
        imshow(Ibwd2f)
        title('Final dark structures which are outside of neural cell regions (Ibwd2f)');
        drawnow;
        waitforbuttonpress;
        close all;
    end

    % #########
    % Step 3:
    % #########

    % Add the resulting dark structures to the neural cell segmentation
    % - identify boundary pixels of dark structures (Ibwd1f and Ibwd2f)
    % - identify internal pixels of dark structures (Ibwd1f and Ibwd2f)
    % - set Ibwn to 0 at the boundary pixels
    % - set Ibwn to 1 at the internal pixels
    Ibwf = Ibwn + Ibwd1f + Ibwd2f;
    Bd1f = bwboundaries(Ibwd1f, 4);
    for i=1:numel(Bd1f)  % i = component id
        % Boundary pixels
        bp1 = Bd1f{i};
        Q1  = size(bp1,1);
        for j=1:Q1          
            Ibwf(bp1(j,1),bp1(j,2)) = logical(0);
        end
    end
    Bd2f = bwboundaries(Ibwd2f, 4);
    for i=1:numel(Bd2f)
        bp2 = Bd2f{i};
        Q2  = size(bp2,1);
        for j=1:Q2            
            Ibwf(bp2(j,1), bp2(j,2)) = logical(0);
        end
    end

    % Display final labels
    if (options.dispOn || options.dispOn2)
        % Label components in the inverse BW image
        L = generate_labels_matrix(Ibwf, options.verbose);

        % Color components
        Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
        figure
        imshow(Lrgb)
        title('Final colored watershed label matrix (Lrgb)')
        drawnow;
        waitforbuttonpress;
        if options.closeAll
            close all;
        end
    end

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_Segmentation2D(inputName,fracBlack,[,parameters])');
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
function Ibw = add_boundary_padding(Ibw, value)
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

