%
% Copyright (C) 2015 by Howard Hughes Medical Institute.
%

%
% Purpose: detect the location of "dark" structures, 
% such as mitochondria, T-,E-bars, and glial cells 
% in a 2D intensity image  
% 

function MS_S2D_DetectDarkStructures2D(varargin)   
warning off;
%%
clc;
try
    image_file = varargin{1};
catch
    output_usage_message();
    return
end

default_frac_black = 0.5;
if nargin > 1
    frac_black = double(varargin{2});
else
    frac_black = default_frac_black;
end 

% Original image
I = imread(image_file);
Igr = mat2gray(I(:,:,1));
disp(['min(Igr)=' num2str(min(min(Igr))) ' max(I)=' num2str(max(max(Igr))) ...
      ' mean2(I)=' num2str(mean2(Igr)) ' std2(Igr)=' num2str(std2(Igr))]);
In = round(double(Igr)*255.0/max(max(double(Igr))));
disp(['min(In)=' num2str(min(min(In))) ' max(In)=' num2str(max(max(In))) ...
      ' mean2(In)=' num2str(mean2(In)) ' std2(In)=' num2str(std2(In))]);
figure(1)
imshow(Igr);
title('Original image (Igr)');
waitforbuttonpress;
%close all;

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
Ibw  = adjust_edges(Ibw, 1);
%Ibw  = imfill(Ibw, 'holes');
disp(['my frac_black=' num2str(sum(sum(Ibw == 0))/numel(Ibw))]);
figure(2)
imshow(Ibw), title('Black-white image without holes (Ibw)')
waitforbuttonpress;
%close all;

% open  complement
Ibwc = imcomplement(Ibw);
Ibwc = bwareaopen(Ibwc,200);
%Ibwc = adjust_edges(Ibwc, 0);
Ibwc = imcomplement(Ibwc);
%Ibwc = imfill(Ibwc, 'holes');
figure(3)
imshow(Ibwc), title('Complement of dilated black-white image (Ibwc)')
waitforbuttonpress;
%close all;

% Dilated BW image
%bwboundaries();
Ibwcd = imdilate(Ibwc, strel('disk', 2));
figure(4)
imshow(Ibwcd), title('Dilated, inversed, filled and eroded black-white image (Ibwd)')
waitforbuttonpress;
%close all;

% Inversed 
Ibwcdc = imcomplement(Ibwcd);
Ibwcdc = imerode(Ibwcdc, strel('disk', 2));
Ibwcdc = bwareaopen(Ibwcdc,200);
Ibwcdc = imfill(Ibwcdc, 'holes');
Ibwcdc = imdilate(Ibwcdc, strel('disk', 2));
Ibwcdc = imfill(Ibwcdc, 'holes');
Ibwcdc = adjust_edges(Ibwcdc, 0);
figure(5)
imshow(Ibwcdc), title('Dilated, inversed, filled and eroded black-white image (Ibwd)')
waitforbuttonpress;
%close all;

% Label components in the inverse BW image
CC = bwconncomp(Ibwcdc); % NOTE: this fucntion does not label holes
Nc = CC.NumObjects;
imsize = CC.ImageSize;

disp(['imsize=' num2str(imsize) ' num_comp=' num2str(Nc) ]);
L = zeros(imsize);
for k = 1:Nc
    for m=1:numel(CC.PixelIdxList{k})
        lin_ind = CC.PixelIdxList{k}(m);
        [r, c] = ind2sub(CC.ImageSize, lin_ind);
%       disp(['lin_ind=' num2str(lin_ind) ' r=' num2str(r) ' c=' num2str(c)]);
        L(r,c) = k;
    end
end

Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure(6)
imshow(Lrgb)
title('Colored watershed label matrix (Lrgb)')
waitforbuttonpress;
close all;

% -----------------------------------------------------------------------------

% Add zero padding at the adges,
% to ensure thgat that entire image != one BW component
function Ibw = adjust_edges(Ibw, value)
    size1 = size(Ibw, 1);
    size2 = size(Ibw, 2);

    i = 1;
    while sum(Ibw(i,:)) >= size2*0.9
        Ibw(i,:) = logical(value);
        i = i+1;
    end

    i = 1;
    while sum(Ibw(:,i)) >= size1*0.9
        Ibw(:,i) = logical(value);
        i = i+1;
    end

    i = 0;
    while sum(Ibw(size1 - i,:)) >= size2*0.9
        Ibw(size1 - i,:) = logical(value);
        i = i+1;
    end

    i = 0;
    while sum(Ibw(:,size2 - i)) >= size1*0.9
        Ibw(:,size2 - i) = logical(0);
        i = i+1;
    end

% -----------------------------------------------------------------------------

function output_usage_message()
    disp('Usage: MS_S2D_DetectDarkStructures2D(image_file [ ,fraction_of_black ])');
    return;
