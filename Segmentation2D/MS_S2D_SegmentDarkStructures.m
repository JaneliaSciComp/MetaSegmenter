%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: detect the location of "dark" structures, 
% such as mitochondria, T-,E-bars, and glial cells 
% in a 2D intensity image  
% 

function Ibwd = MS_S2D_SegmentDarkStructures(inputName,fracBlack,fracBlack2,varargin)   
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

    if length(size(I) > 2)
        I = I(:,:,1);
    end
    disp(['In MS_S2D_SegmentDarkStructures: size(I)=' num2str(size(I))]);
 
    % Original image
    if options.dispOn
        MS_S2D_ShowImage(I, 'Original image (I)', options);
    end

    im_size = size(I);
    Igr = mat2gray(I);
    clear I;

    [M_thr, M] = MS_S2D_GetThresholdIntensity(Igr, options);      
    Ibwd = segment_dark_structures(Igr, M_thr, M_thr, options);
    Ibwd = MS_S2D_AddBoundaryPadding(Ibwd, 0);

    if options.dispOn
        MS_S2D_ShowImage(Ibwd, 'Final (dilated, inversed, filled and eroded) black-white image (Ibwd)', options);
    end

    % Display/outpu labels
    if (options.dispOn || options.dispOn2 || length(options.outSeg) > 0 ...
                       || options.RGB)
        % Label components in the inverse BW image
        L = MS_S2D_GenerateLabelsMatrix(Ibwd, options.verbose);
        if length(options.outSeg) > 0
            imwrite(L, options.outSeg);
        end

        % Color components
        if options.dispOn2 || options.RGB
            Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
            MS_S2D_ShowImage(Lrgb, 'Colored labels of dark structures (Lrgb)', options);
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

    Ibwd = [];
    if length(options.mitoPr) > 0
        % Use probabilities to detect dartk structures
        norm_mitoPr = normalize_probabilities(MS_S2D_ReadProbabilities(options.mitoPr));
        weight_signals_by_probabilities = 1;
        if weight_signals_by_probabilities
            Igr = weight_image_by_probabilities(Igr, norm_mitoPr, options);
            Igrbw = im2bw(Igr, 0.2);
        else
            Igr = mat2gray(round(norm_mitoPr*255));
            Igrbw = im2bw(Igr, 0.2);
        end
        if options.dispOn | options.dispOn2
            MS_S2D_ShowImage(Igr, 'Normalized mito probability ', options);
        end 
%       if options.hist
%           figure
%           imhist(Ipr);
%           hold on;
%           title('Histogram of mito probabilities');
%           waitforbuttonpress;
%       end
        clear Igr;
        if options.dispOn | options.dispOn2
            MS_S2D_ShowImage(Igrbw, 'Thresholded mitoProb ', options);
        end
        Igrbw = imdilate(Igrbw,strel('disk', 8));
        Igrbw  = imfill(Igrbw, 'holes');
        Igrbw = imerode(Igrbw, strel('disk', 12));
        Igrbw = bwareaopen(Igrbw,round(400));
        if options.dispOn | options.dispOn2
            MS_S2D_ShowImage(Igrbw, 'After dilation/erosion ', options);
        end
        Ibwd = Igrbw;
        clear Igrbw;
    else
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
%       if options.dispOn
%           MS_S2D_ShowImage(Ibw, 'Black-white image without holes (Ibw)', options);
%       end

    % Handling complement
        Ibwc = imcomplement(Ibw);
        clear Ibw;
        Ibwc = imerode(Ibwc, strel('disk', 2));
        Ibwc = bwareaopen(Ibwc,round(300*scale));
%       if options.dispOn | options.dispOn2
%           MS_S2D_ShowImage(Ibwc, 'Complement of dilated black-white image (Ibwc)', options);
%       end

        % Dilated BW image
        Ibwc = imfill(Ibwc, 'holes');
        Ibwc = imcomplement(Ibwc);
        Ibwcd = imdilate(Ibwc, strel('disk', 2));
        clear Ibwc;
        if options.dispOn
            MS_S2D_ShowImage(Ibwcd, 'Dilated, inversed, filled and eroded black-white image (Ibwd)', options);
        end

        % Inversed 
        Ibwcd  = imerode(Ibwcd, strel('disk', 4));
        Ibwcdc = bwareaopen(Ibwcd,round(300*scale));
        Ibwcdc = imfill(Ibwcdc, 'holes');
        Ibwcdc = imdilate(Ibwcdc, strel('disk', 3));
        Ibwd  = imfill(Ibwcdc, 'holes');
        clear Ibwcdc;
    
        Ibwd = imresize(Ibwd, orig_imsize);
%       Ibwd   = MS_S2D_AddBoundaryPadding(Ibwd, 0);
    end
    if options.dispOn | options.dispOn2
        MS_S2D_ShowImage(Ibwd, 'Final black-white image (Ibwd)', options);
    end

% -----------------------------------------------------------------------------

function norm_probs = normalize_probabilities(probs)
    max_prob = max(max(probs));
    norm_probs = probs/max_prob;

% -----------------------------------------------------------------------------

function Iw = weight_image_by_probabilities(I, membPr, options)
    Iw = mat2gray(double(imcomplement(I)).*membPr);

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_SegmentDarkStructures(image_file,fraction_of_black[,parameters])');
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
    disp('    RGB          - whether or not to show colored labels (default = no)');                  
    disp('    verbose      - weather or not to increase verbosity of output (default=0)');
    return;

