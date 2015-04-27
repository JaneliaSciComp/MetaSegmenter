%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%
% Purpose: combine regions corresponding to neurons and dark structures
% in order to produce a final segmentation
function [] = MS_S2D_Segmentation2D(inputName,fracBlack,fracBlack2,varargin)   
    debug = 0;
    if nargin < 3
        output_usage_message();
        return
    end

    try
        p = MS_S2D_InputParser;              % call constructor
        p.parse(fracBlack,fracBlack2,varargin{:});
        options = MS_S2D_ExtractOptions(p, fracBlack, fracBlack2);
    catch
        output_usage_message();
        return
    end

    fracBlack  = options.fracBlack;
    fracBlack2 = options.fracBlack2;
    disp(['In MS_S2D_Segmentation2D: options.hist=' num2str(options.hist)]);
    options.debug = debug;

    % Input image data  
    if exist(inputName, 'file') == 2   % input is image file
        I = imread(inputName); 
    else
        I = inputName;                 % input is image
    end
    
    if length(size(I)) > 2
        I = squeeze(I(:,:,1));
    end
    
    if (options.dispOn || options.dispOn2 || options.debug)
        MS_S2D_ShowImage(I(:,:), 'Original image (I)', options);
    end

    imwrite(I, 'orig_image_MS_S2D_Segmentation2D.png');

    Ibwf = segment_image(I, fracBlack, fracBlack2, options);
    if options.padding
        Ibwf = MS_S2D_AddBoundaryPadding(Ibwf, 0);
    end
    
     % Optionally output BW file
    if length(options.outBW) > 0
        imwrite(Ibwf, options.outBW);
    end

    L = MS_S2D_GenerateLabelsMatrix(Ibwf, options.verbose);
    
    if (options.dispOn || options.dispOn2 || options.debug ...
                       || length(options.outRGB) > 0);          
        Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
        if (options.dispOn || options.dispOn2 || options.debug)
            MS_S2D_ShowImage(Lrgb, 'Final segmentation colored labels (Lrgb)', options);
        end
        if length(options.outRGB) > 0
            imwrite(Lrgb, options.outRGB);
        end
    end

    % Optionally output segmentation file 
    if length(options.outSeg) > 0  
        imwrite(L, options.outSeg);
    end

% -------------------------------------------------------------------------------

function Ibwf = segment_image(I, fracBlack, fracBlack2, options)
    % Black-white images
    %
    if options.compiled_code
        Ibwn = MS_S2D_SegmentNeurons(I, ...
            fracBlack, fracBlack2,...
            'nx', options.nx, 'ny', options.ny,...
            'dx', options.dx, 'dy', options.dy,...
            'maxSize', options.maxSize, ...
            'verbose', options.verbose);
        Ibwd = MS_S2D_DetectDarkStructures2D(I, ...
            fracBlack, fracBlack2,...
            'nx', options.nx, 'ny', options.ny,...
            'dx', options.dx, 'dy', options.dy,...
            'maxSize', options.maxSize, ...
            'verbose', options.verbose);
    else 
        Ibwn = MS_S2D_SegmentNeurons(I, fracBlack, fracBlack2,...
            'nx', options.nx, 'ny', options.ny,...
            'dx', options.dx, 'dy', options.dy,...
            'sx', options.sx, 'sy', options.sy,...
            'maxSize', options.maxSize, ...
            'verbose', options.verbose,...
            'dispOn', options.dispOn,...
            'dispOn2', options.dispOn2,...
            'closeAll',options.closeAll, ...
            'hist', options.hist);
        Ibwd = MS_S2D_DetectDarkStructures2D(I, fracBlack, fracBlack2,...
            'nx', options.nx, 'ny', options.ny,...
            'dx', options.dx, 'dy', options.dy,...
            'sx', options.sx, 'sy', options.sy,...
            'maxSize', options.maxSize, ...
            'verbose', options.verbose,...
            'dispOn', options.dispOn,...
            'dispOn2', options.dispOn2,...
            'closeAll',options.closeAll, ...
            'hist', options.hist);
    end
    clear I;

%   Ibwn = imfill(Ibwn, 'holes'); % instead, better to make a continuous, piecewise linear threshold
    if options.padding
        Ibwn = MS_S2D_AddBoundaryPadding(Ibwn, 0);
        Ibwd = MS_S2D_AddBoundaryPadding(Ibwd, 0);
    end
    if options.debug
        MS_S2D_ShowImage(Ibwn, 'Initially detectred neural cells (Ibwn)', options);
        MS_S2D_ShowImage(Ibwd, 'Initially detectred dark structures (Ibwd)', options);
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
    clear Ibwn1;
    clear Ibwn2;
    Ibwd11 = bwareaopen(Ibwd1,200);

    % Identify the dark strutures each of which is ALMOST entirely contained
    % within a single neural cell region.
    Ibwde  = imerode(Ibwd - Ibwd11, strel('disk', 1));
    Ibwn1 = (Ibwnd | Ibwde) - Ibwde;
    Ibwn2 = imfill(Ibwn1, 'holes');
    Ibwd12 = Ibwn2 - Ibwn1;
    clear Ibwn1;
    clear Ibwn2;
    Ibwd12  = bwareaopen(Ibwd12,200);
    Ibwd1f = Ibwd11 | Ibwd12;
    if options.debug 
        MS_S2D_ShowImage(Ibwd1f, 'Final dark structures which are entirely contained in neural cell regions (Ibwd1f)', options);
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
    clear Ibwnd;
    Ibwd2 = bwareaopen(Ibwd2, 400);
    Ibwd2f = imfill(Ibwd2, 'holes');
    clear Ibwd2;
%   Ibwd1 = imclose(Ibwd1, strel('disk', 2));

    if options.debug
         MS_S2D_ShowImage(Ibwd2f,'Final dark structures which are outside of neural cell regions (Ibwd2f)', options);
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
    if options.padding
        Ibwf = MS_S2D_AddBoundaryPadding(Ibwf, 0);
    end

    Bd1f = bwboundaries(Ibwd1f, 4);
    clear Ibwd1f;
    for i=1:numel(Bd1f)  % i = component id
        % Boundary pixels
        bp1 = Bd1f{i};
        Q1  = size(bp1,1);
        for j=1:Q1          
            Ibwf(bp1(j,1),bp1(j,2)) = logical(0);
        end
    end
    clear Bd1f;
    Bd2f = bwboundaries(Ibwd2f, 4);
    clear Ibwd2f;
    for i=1:numel(Bd2f)
        bp2 = Bd2f{i};
        Q2  = size(bp2,1);
        for j=1:Q2            
            Ibwf(bp2(j,1), bp2(j,2)) = logical(0);
        end
    end
    clear Bd2f;
%   Ibwf = imfill(Ibwf, 'holes'); % is it possible to fill only those holes which do not contain childs?
    return

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_Segmentation2D(inputName,fracBlack,fracBlack2[,parameters])');
    disp('Required arguments:');
    disp('    inputName    - name of input data (file or image)');
    disp('    fracBlack    - desired fraction of black for black/wite conversion');
    disp('Optional parameters, specified as parameter-value pairs:');
    disp('    dispOn       - weather or not to display covariance matrix (default=0)');
    disp('    closeAll     - close previous image when displayong next (default=1)');
    disp('    nx, ny       - numbers of section to which subdivide the image (defaults=1)');
    disp('    dx, dy       - # pixels of fragment overlap (default = 50)');
    disp('    sx, sy       - ids of the particular section to be processed/shown (defaults=1)');
    disp('    maxSize      - maximum allowed size of image subsection to be processed (default=Inf)');
    disp('    outBW        - name of output black/white image file (default='', no output)');
    disp('    outSeg       - name of output segmentation file (default='', no output)');
    disp('    outRGB       - name of output colored labels file (default='', no output)');
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return;

