%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

function MS_S2D_SegmentNeurons2D(varargin)   
warning off;
%%
clc;
try
    image_file = varargin{1};
catch
    output_usage_message();
    return
end

frac_black = 0.34;
if nargin > 1
    frac_black = double(varargin{2});
end 

% Original image    
I = imread(image_file); 
Igr = mat2gray(I(:,:,1));
disp(['min(I)=' num2str(min(min(Igr))) ' max(Igr)=' num2str(max(max(Igr))) ...
      ' mean2(Igr)=' num2str(mean2(Igr)) ' std2(Igr)=' num2str(std2(Igr))]);
In = round(double(Igr)*255.0/max(max(double(Igr))));
disp(['min(In)=' num2str(min(min(In))) ' max(In)=' num2str(max(max(In))) ...
      ' mean2(In)=' num2str(mean2(In)) ' std2(In)=' num2str(std2(In))]);
imshow(Igr);
title('Original image (Igr)');
waitforbuttonpress;
close all;

% Black-white image
%
% Determine intensity threshold for conversion to black-white image
% NOTE: imhist assumes that I is a grayscale image, with range of values [0, 1]
% NOTE2: edges are computed incorrectly by Matlab's imhist function,
%        so a custom function has been implemented for compuring the edges
threshold = MS_S2D_GetThresholdIntensity(Igr, frac_black, 1001);
disp(['final threshold_intensity=' num2str(threshold)]);
disp(['size(Igr)=' num2str(size(Igr)) ' class(Igr)=' class(Igr(1,1))]);
Ibw = im2bw(Igr, threshold);  
disp(['size(Ibw)=' num2str(size(Ibw)) ' class(Ibw)=' class(Ibw(1,1))]);
Ibw = bwareaopen(Ibw,20); 
disp(['my frac_black=' num2str(sum(sum(Ibw == 0))/numel(Ibw))]);
figure
imshow(Ibw), title('Black-white image (Ibw)')
drawnow;
waitforbuttonpress;
close all;

% Fill holes and open gaps
Ibw = adjust_edges(Ibw);
Ibwf = imfill(Ibw, 'holes');
Ibwf = bwareaopen(Ibwf,20);
imshow(Ibwf)
drawnow;
waitforbuttonpress;
close all;

% Label components in the inverse BW image
L = generate_labels_matrix(Ibwf);

% Color components
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure
imshow(Lrgb)
title('Colored watershed label matrix (Lrgb)')
drawnow;
waitforbuttonpress;
close all;

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_SegmentNeurons2D(image_file [ ,fraction_of_black ])');
    return;

% -----------------------------------------------------------------------------

function L = generate_labels_matrix(Ibwf)
    CC = bwconncomp(Ibwf); % NOTE: this function does not label holes
    Nc = CC.NumObjects;
    imsize = CC.ImageSize;

    disp(['imsize=' num2str(imsize) ' num_comp=' num2str(Nc) ]);
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
function Ibw = adjust_edges(Ibw)
    size1 = size(Ibw, 1);
    size2 = size(Ibw, 2);

    i = 1;
    while sum(Ibw(i,:)) >= size2*0.9
        Ibw(i,:) = logical(0);
        i = i+1;
    end

    i = 1;
    while sum(Ibw(:,i)) >= size1*0.9
        Ibw(:,i) = logical(0);
        i = i+1;
    end

    i = 0;
    while sum(Ibw(size1 - i,:)) >= size2*0.9
        Ibw(size1 - i,:) = logical(0);
        i = i+1;
    end

    i = 0;
    while sum(Ibw(:,size2 - i)) >= size1*0.9
        Ibw(:,size2 - i) = logical(0);
        i = i+1;
    end

