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

frac_black = 0.42;
if nargin > 1
    frac_black = double(varargin{2});
end 

% Original image    
I = imread(image_file); I = I(:,:,1);
imshow(I);
title('Original image (I)');
waitforbuttonpress;
close all;

% Black-white image
num_bins = 100;
[N,edges] = imhist(I,num_bins);
Nsum = sum(N);
Ncurr = 0;
threshold_frac = 0;
i = 1;
while i<= num_bins && threshold_frac < frac_black
    Ncurr = Ncurr + N(i);
    threshold_frac = Ncurr/Nsum; 
    i = i+1;   
end
disp(['threshold_frac=' num2str(threshold_frac)]);
Igr = mat2gray(I);
Ibw = im2bw(Igr, threshold_frac);  
Ibw = bwareaopen(Ibw,200); 
disp(['black/white=' num2str(sum(sum(Ibw == 0))/sum(sum(Ibw == 1)))]);
figure
imshow(Ibw), title('Black-white image (Ibw)')
drawnow;
waitforbuttonpress;
close all;

% Dilated BW image
Ibwd = imdilate(Ibw, strel('disk', 3));
imshow(Ibwd), title('Dilated black-white image (Ibwd)')
drawnow;
waitforbuttonpress;
close all;

% open  complement
Ibwdc = imcomplement(Ibwd);
Ibwdc = bwareaopen(Ibwdc,200);
imshow(Ibwdc), title('Complement of dilated black-white image (Ibwdc)')
drawnow;
waitforbuttonpress;
close all;

Ibwdco = imopen(Ibwdc, strel('disk', 6));
%Ibwdo = imcomplement(Ibwdco);
Ibwdco = bwareaopen(Ibwdco,400);
imshow(Ibwdco), title('Eroded black-white image (Ibwdo)')
drawnow;
waitforbuttonpress;
close all;

%Ibwc = imcomplement(Ibw);
%Ibwc = bwareaopen(Ibwc,120);
%figure
%imshow(Ibwc), title('Complement of black-white image (Ibwc)')
%drawnow;
%waitforbuttonpress;
%close all;

% Label components in the inverse BW image
CC = bwconncomp(Ibw); % NOTE: this fucntion does not label holes
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
figure
imshow(Lrgb)
title('Colored watershed label matrix (Lrgb)')
drawnow;
waitforbuttonpress;
close all;


function output_usage_message()
    disp('Usage: MS_S2D_DetectDarkStructures2D(image_file [ ,fraction_of_black ])');
    return;
