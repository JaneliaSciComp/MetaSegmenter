%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%
% Purpose: combine regions corresponding to neurons and dark structures
% in order to produce a final segmentation
function [] = MS_S2D_Segmentation2D(inputName,fracBlack,fracBlack2,varargin)   
    if nargin < 3
        output_usage_message();
        return
    end

%   try
        p = MS_S2D_InputParser;              % call constructor
        p.parse(fracBlack,fracBlack2,varargin{:});
        options = MS_S2D_ExtractOptions(p, fracBlack, fracBlack2);
%   catch
%       output_usage_message();
%       return
%   end

    fracBlack  = options.fracBlack;
    fracBlack2 = options.fracBlack2;
    disp(['In MS_S2D_Segmentation2D: options.hist=' num2str(options.hist)]);
%   options.heatmap = 0;

    % Input image data  
    if exist(inputName, 'file') == 2   % input is image file
        I = imread(inputName); 
    else
        I = inputName;                 % input is image
    end
    
    if length(size(I)) > 2
        I = squeeze(I(:,:,1));
    end
   
    max_I = max(max(I));
    if max_I > 255
        I = mat2gray(double(I)/max_I);
    end
    disp(['In MS_S2D_Segmentation2D max(I)=' num2str(max(max(I)))]);

%   if (options.dispOn || options.dispOn2 || options.debug)
%       MS_S2D_ShowImage(I(:,:), 'Original image (I)', options);
%   end
    
    Ibwf = segment_image(I, fracBlack, fracBlack2, options);
    if options.padding
        Ibwf = MS_S2D_AddBoundaryPadding(Ibwf, 0);
    end
    
     % Optionally output BW file
    if length(options.outBW) > 0 
        imwrite(Ibwf, options.outBW);
    end

    if options.RGB || length(options.outRGB) > 0 || length(options.outSeg) > 0
        Ibwf = imerode(Ibwf, strel('disk', 1));
        L = MS_S2D_GenerateLabelsMatrix(Ibwf, options.verbose);
        if length(options.outSeg) > 0
            imwrite(L, options.outSeg);
        end
        Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
        if options.RGB
            MS_S2D_ShowImage(Lrgb, 'Final segmentation colored labels (Lrgb)', options);
        end
        if length(options.outRGB) > 0
            disp(['Producing the RGB file ' options.outRGB ' ...']);
            imwrite(Lrgb, options.outRGB);
        end
    end

% -------------------------------------------------------------------------------

function Ibwf = segment_image(I, fracBlack, fracBlack2, options)
    % Black-white images
    %
    if options.compiled_code
        Ibwn = MS_S2D_SegmentNeuralMembranes(I, ...
            fracBlack, fracBlack2,...
            'nx', options.nx, 'ny', options.ny,...
            'dx', options.dx, 'dy', options.dy,...
            'maxSize', options.maxSize, 'minSize', options.minSize, ...
            'thr', options.thr, 'thrMembPr', options.thrMembPr,...
            'membPr', options.membPr, 'mitoPr', options.mitoPr, ...
            'mitoMembPr', options.mitoMembPr, 'useMembPr', options.useMembPr, ...
            'minWin', options.minWin, 'maxWin', options.maxWin, ...
            'outThr', options.outThr, 'showIter', options.showIter, ...
            'debug', options.debug, 'incRef', options.incRef, 'heatmap', options.heatmap, ...
            'minMitoFr', options.minMitoFr, 'maxMitoFr', options.maxMitoFr, ...
            'subWin', options.subWin, 'numRef', options.numRef, ...
            'verbose', options.verbose, 'vesicles', options.vesicles);
        if options.darkStr
            Ibwd = MS_S2D_SegmentDarkStructures(I, ...
                fracBlack, fracBlack2,...
                'nx', options.nx, 'ny', options.ny,...
                'dx', options.dx, 'dy', options.dy,...
                'membPr', options.membPr, 'mitoPr', options.mitoPr, ...
                'thr',   options.thr2,  ...
                'subWin', options.subWin, ...
                'verbose', options.verbose);
        else
            Ibwd = [];
        end
    else 
        Ibwn = MS_S2D_SegmentNeuralMembranes(I, fracBlack, fracBlack2,...
            'nx', options.nx, 'ny', options.ny,...
            'dx', options.dx, 'dy', options.dy,...
            'membPr', options.membPr, 'mitoPr', options.mitoPr, 'thrMembPr', options.thrMembPr,...
            'mitoMembPr', options.mitoMembPr, 'useMembPr', options.useMembPr, ...
            'sx', options.sx, 'sy', options.sy, ...
            'maxSize', options.maxSize, 'minSize', options.minSize, ...
            'thr',    options.thr,  'verbose', options.verbose,...
            'minWin', options.minWin, 'maxWin', options.maxWin, 'heatmap', options.heatmap,  ...
            'debug', options.debug,  'incRef', options.incRef, ...
            'minMitoFr', options.minMitoFr, 'maxMitoFr', options.maxMitoFr, ...
            'dispOn', options.dispOn, 'dispOn2', options.dispOn2,...
            'outThr', options.outThr, ...
            'closeAll',options.closeAll, 'showIter', options.showIter, ...
            'subWin', options.subWin, 'numRef', options.numRef, ...
            'hist', options.hist, 'vesicles', options.vesicles);
        if options.darkStr
            Ibwd = MS_S2D_SegmentDarkStructures(I, fracBlack, fracBlack2,...
                'nx', options.nx, 'ny', options.ny,...
                'dx', options.dx, 'dy', options.dy,...
                'membPr', options.membPr, 'mitoPr', options.mitoPr, ...
                'mitoMembPr', options.mitoMembPr, ...
                'sx', options.sx, 'sy', options.sy,...
                'thr',   options.thr2, ...
                'verbose', options.verbose, 'subWin', options.subWin, ...
                'dispOn', options.dispOn,...
                'dispOn2', options.dispOn2,...
                'closeAll',options.closeAll, ...
                'hist', options.hist);
        else
            Ibwd = [];
        end
    end
    clear I;

%   Ibwn = imfill(Ibwn, 'holes'); % instead, better to make a continuous, piecewise linear threshold
    if options.padding
        Ibwn = MS_S2D_AddBoundaryPadding(Ibwn, 0);
        if options.darkStr
            Ibwd = MS_S2D_AddBoundaryPadding(Ibwd, 0);
        end
    end

    if length(options.mitoPr) > 0
        % Whiten out the area corresponding to the detecetd dark structures  
        Ibwf = Ibwn;
        clear Ibwn;
        if options.darkStr
            Ibwd_int = imerode(Ibwd, strel('disk', 2));
            Ibwd_bound = Ibwd - Ibwd_int;
            Ibwf(Ibwd_int   > 0) = 1; % interior part of dark structures
            Ibwf(Ibwd_bound > 0) = 0; % boundary of dark structures 
        end
        Ibwf = bwareaopen(Ibwf,100);
    else
        if options.darkStr
            % #########
            % Step 1:
            % #########

            % Identify the dark strutures each of which is entirely contained 
            % within a single neural cell region.   
            disp(['Step 1: dark strutures entirely contained within a single neural cell ...']);
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
            disp(['Step 2: dark strutures  ALMOST entirely contained...']);
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
%           Ibwd1 = imclose(Ibwd1, strel('disk', 2));

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
%           Ibwf = imfill(Ibwf, 'holes'); % is it possible to fill only those holes which do not contain childs?

        else
            Ibwf = Ibwn;
        end

        if options.padding
            Ibwf = MS_S2D_AddBoundaryPadding(Ibwf, 0);
        end
    end
    return

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_Segmentation2D(inputName,fracBlack,fracBlack2[,parameters])');
    disp('Required arguments:');
    disp('    inputName    - name of input data (file or image)');
    disp('    fracBlack    - desired fraction of black for black/wite conversion');
    disp('    fracBlack2   - desired fraction of black for black/wite conversion');
    disp('Optional parameters, specified as parameter-value pairs:');
    disp('    dispOn       - weather or not to display covariance matrix (default=0)');
    disp('    closeAll     - close previous image when displayong next (default=1)');
    disp('    nx, ny       - numbers of section to which subdivide the image (defaults=1)');
    disp('    dx, dy       - # pixels of fragment overlap (default = 50)');
    disp('    sx, sy       - ids of the particular section to be processed/shown (defaults=1)');
    disp('    membPr       - path to HDF5 file containing membrane probabilities');
    disp('    mitoPr       - path to HDF5 file containing mitochondria probabilities');
    disp('    mitoMembPr   - path to HDF5 file containing mitochondrial membrane probabilities');
    disp('    useMembPr    - whether or not to use membrane probabilities instead of original signals');
    disp('    outBW        - name of output black/white image file (default='', no output)');
    disp('    outSeg       - name of output segmentation file (default='', no output)');
    disp('    RGB          - whether or not to show colored labels (default = no)');
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return;

